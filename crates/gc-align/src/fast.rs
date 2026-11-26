//! Schnellweg für leichte Reads.
//!
//! Die meisten Reads eines menschlichen Genoms passen in voller Länge ohne Lücke an genau eine
//! Stelle, mit höchstens einer Abweichung (HG002, GIAB-BAM von bwa auf chr20: 74 % ohne,
//! 16 % mit einer oder zwei Abweichungen, 92,5 % mit MAPQ 60). Für sie liefert das volle
//! BWA-MEM-Verfahren (Seed-Suche mit Nachsuche, Ketten, Verlängern) genau einen Bereich über den
//! ganzen Read. Der Schnellweg sucht diesen Bereich direkt:
//!
//! 1. Den Read in Stücke von mindestens 21 Basen teilen (7 bei 150 Basen) und jedes Stück im
//!    FM-Index suchen (beide Stränge).
//! 2. Hat ein Stück mehrere Fundstellen, gibt es eine zweite Kopie, und der Read geht den vollen
//!    Weg. Höchstens ein Stück darf gar keine haben.
//! 3. Ein Stück mit genau einer Fundstelle gibt den Anfang. Der ganze Read wird dort Base für
//!    Base mit der Referenz verglichen: höchstens eine Abweichung, kein N, innerhalb eines Contigs.
//! 4. Jedes Stück ohne die Abweichung muss genau eine Fundstelle haben (die bei diesem Anfang),
//!    das Stück mit der Abweichung keine.
//!
//! Mit höchstens einer Abweichung ist das ungelückte Alignment über den ganzen Read auch das
//! beste: Jede Lücke kostet mindestens 7 Punkte, die Abweichung nur 5, und das Abschneiden eines
//! Endes lohnt nie (bwa verlängert bis zum Ende, wenn das höchstens 5 Punkte kostet).
//!
//! Warum die MAPQ gleich bleibt: bwa senkt sie bei 150 Basen erst unter 60, wenn eine zweite
//! Kopie höchstens etwa 16 Punkte schlechter ist, also höchstens vier Ereignisse (Abweichung oder
//! Lücke) vom Read entfernt. Bei 7 Stücken passt eine solche Kopie in mindestens drei Stücken
//! genau; eines davon ist nicht das Stück mit der Abweichung und hat dann zwei Fundstellen.
//! Dasselbe gilt für Kopien bis sechs Ereignisse, die bwa als `XA` meldet.
//!
//! Was sich ändert: `XS` bei Reads, für die das volle Verfahren nur zufällige kurze Treffer
//! anderswo findet (Score unter 0,8 des besten). Auf chr20 (HG002, 843.128 Reads) war das die
//! einzige Abweichung, bei 9 % der Zeilen; MAPQ, Flags, Orte und CIGAR waren alle gleich.

use crate::region::Region;
use crate::sam::AlignOptions;
use gc_index::FmIndex;
use std::sync::atomic::{AtomicU64, Ordering};

/// Zähler des Schnellwegs, für das Protokoll.
pub struct FastCounts {
    pub tried: AtomicU64,
    pub taken: AtomicU64,
    pub repeat: AtomicU64,
    pub too_few_pieces: AtomicU64,
    pub disagree: AtomicU64,
    pub verify_failed: AtomicU64,
    /// Weitere Kopien geprüft und schwach genug; zu gut oder zu nah.
    pub copies_checked: AtomicU64,
    pub copy_too_good: AtomicU64,
    /// Nur im Prüfmodus: der volle Weg lieferte genau diesen einen Bereich.
    pub same_as_full: AtomicU64,
    pub differs_from_full: AtomicU64,
}

pub static COUNTS: FastCounts = FastCounts {
    tried: AtomicU64::new(0),
    taken: AtomicU64::new(0),
    repeat: AtomicU64::new(0),
    too_few_pieces: AtomicU64::new(0),
    disagree: AtomicU64::new(0),
    verify_failed: AtomicU64::new(0),
    copies_checked: AtomicU64::new(0),
    copy_too_good: AtomicU64::new(0),
    same_as_full: AtomicU64::new(0),
    differs_from_full: AtomicU64::new(0),
};

fn bump(c: &AtomicU64) {
    c.fetch_add(1, Ordering::Relaxed);
}

/// Mindestlänge eines Stücks (etwas über der kleinsten Seed-Länge 19 von bwa).
const PIECE_MIN: usize = 21;

/// Wie der Schnellweg eingesetzt wird.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FastMode {
    /// Aus: jeder Read geht den vollen Weg (Ausgabe wie bisher).
    #[default]
    Off,
    /// An: leichte Reads nehmen den Schnellweg, auch wenn einzelne Stücke an wenigen weiteren
    /// Stellen passen und diese Kopien nachweislich schwach genug sind.
    On,
    /// Einfach: wie `On`, aber jedes Stück mit mehr als einer Fundstelle schickt den Read den
    /// vollen Weg (Stand vor den Kopien).
    Simple,
    /// Prüfen: beide Wege rechnen, das Ergebnis des vollen Wegs ausgeben und zählen, wie oft der
    /// Schnellweg genau denselben Bereich geliefert hätte.
    Check,
}

/// Längen der Seeds, die das volle Verfahren an einer Abweichung bei `mm` hätte: die exakten
/// Stücke links und rechts davon, wenn sie mindestens `min_len` lang sind.
fn seed_coverage(len: i32, mismatch: Option<i32>, min_len: i32) -> (i32, i32) {
    match mismatch {
        None => (len, len),
        Some(m) => {
            let (a, b) = (m, len - m - 1);
            let cov = [a, b].into_iter().filter(|&x| x >= min_len).sum();
            (cov, a.max(b))
        }
    }
}

/// Höchstens so viele Fundstellen je Stück, damit der Schnellweg die übrigen Kopien noch prüft.
const MAX_PIECE_HITS: u64 = 8;
/// Höchstens so viele weitere Kopien (Anfänge, die nicht der gewählte sind).
const MAX_COPIES: usize = 4;
/// So weit über den Read hinaus wird bei einer weiteren Kopie ausgerichtet (Lücken bis zur
/// Bandbreite 100 von bwa passen hinein).
const COPY_MARGIN: u64 = 128;

/// Sucht den einen Bereich eines leichten Reads. `None` heißt: voller Weg.
///
/// Mit `copies` dürfen Stücke bis zu [`MAX_PIECE_HITS`] Fundstellen haben. Jede weitere Kopie
/// wird dann mit dem ganzen Read lokal ausgerichtet; der Read nimmt den Schnellweg nur, wenn
/// keine davon so gut ist, dass bwa sie in `XA` meldete, die MAPQ senkte oder von ihr aus den
/// Partner suchte (siehe [`copy_limit`]).
/// Obergrenze für den lokalen Score des Reads (Codes 0..3, ohne N) gegen `target`, ohne
/// Ausrichtung: die Zahl der Read-Positionen, deren 5-mer irgendwo in `target` vorkommt, plus 4.
///
/// Gilt für die Vorgabe-Bewertung (Treffer 1, Fehltreffer 4, Lücke mindestens 6 + 1): Ein
/// lokales Alignment besteht aus exakten Läufen, getrennt durch Fehltreffer oder Lücken. Ein
/// Lauf der Länge `L` bringt `L` Punkte und enthält `max(L - 4, 0)` 5-mere, die alle im Fenster
/// vorkommen; er bringt also höchstens 4 Punkte mehr als seine 5-mer-Starts. Jede Trennung
/// kostet mindestens 4 Punkte und hebt den Überschuss eines Laufs auf. Bleibt höchstens der
/// Überschuss eines Laufs: Score ≤ Starts + 4. Der Read darf kein N enthalten (der Schnellweg
/// nimmt solche Reads nicht).
/// `None`, wenn `target` ein N enthält: Gegen N kostet eine Base nur `ambiguous` (1), die
/// Herleitung gilt dann nicht.
fn kmer_bound(read: &[u8], target: &[u8]) -> Option<i32> {
    let mut mask = [0u64; 16];
    let (mut code, mut valid) = (0usize, 0);
    for &b in target {
        if b > 3 {
            return None;
        }
        code = ((code << 2) | usize::from(b)) & 1023;
        valid += 1;
        if valid >= 5 {
            mask[code >> 6] |= 1 << (code & 63);
        }
    }
    let mut starts = 0;
    let mut code = 0usize;
    for (i, &b) in read.iter().enumerate() {
        code = ((code << 2) | usize::from(b)) & 1023;
        if i >= 4 && mask[code >> 6] >> (code & 63) & 1 == 1 {
            starts += 1;
        }
    }
    Some(starts + 4)
}

/// Gilt [`kmer_bound`] für diese Bewertung?
fn kmer_bound_holds(sc: &crate::extend::Scoring) -> bool {
    sc.matching == 1
        && sc.mismatch >= 4
        && sc.o_del + sc.e_del >= 4
        && sc.o_ins + sc.e_ins >= 4
        && sc.e_del >= 0
        && sc.e_ins >= 0
}

pub fn fast_region(
    index: &FmIndex,
    read: &[u8],
    opt: &AlignOptions,
    copies: bool,
) -> Option<Region> {
    bump(&COUNTS.tried);
    let len = read.len();
    if len < 60 || read.iter().any(|&c| c > 3) {
        bump(&COUNTS.too_few_pieces);
        return None;
    }
    // Stücke von mindestens PIECE_MIN Basen. Eine zweite Kopie, die sich in höchstens
    // `n - 1` Ereignissen (Abweichung oder Lücke) vom Read unterscheidet, passt in mindestens
    // einem Stück genau und fällt dort als weitere Fundstelle auf.
    let n = (len / PIECE_MIN).max(3);
    let cut = |k: usize| k * len / n;
    let mut iv = [(0u64, 0u64); 16];
    let n = n.min(iv.len());
    let max_hits = if copies { MAX_PIECE_HITS } else { 1 };
    let mut many = false;
    if opt.seed.fast_pieces_together {
        let mut pieces: [&[u8]; 16] = [&[]; 16];
        for (k, p) in pieces.iter_mut().enumerate().take(n) {
            *p = &read[cut(k)..cut(k + 1)];
        }
        index.intervals_many(&pieces[..n], &mut iv[..n]);
        for &(lo, hi) in &iv[..n] {
            if hi - lo > max_hits {
                bump(&COUNTS.repeat);
                return None;
            }
            many |= hi - lo > 1;
        }
    } else {
        for (k, slot) in iv.iter_mut().enumerate().take(n) {
            let (lo, hi) = index.interval(&read[cut(k)..cut(k + 1)]);
            if hi - lo > max_hits {
                bump(&COUNTS.repeat);
                return None;
            }
            many |= hi - lo > 1;
            *slot = (lo, hi);
        }
    }
    let count = |k: usize| iv[k].1 - iv[k].0;
    // Anfänge je Stück. Ohne Stück mit mehreren Fundstellen genügt ein Nachschlagen.
    let mut starts: Vec<(u64, usize)> = Vec::new();
    let s = if many {
        for (k, &(lo, hi)) in iv.iter().enumerate().take(n) {
            for row in lo..hi {
                if let Some(s) = index.locate(row).checked_sub(cut(k) as u64) {
                    starts.push((s, k));
                }
            }
        }
        starts.sort_unstable();
        // Der Anfang, den alle Stücke bis auf höchstens eines teilen.
        let mut best: Option<(u64, usize)> = None;
        let mut i = 0;
        while i < starts.len() {
            let j = starts[i..]
                .iter()
                .take_while(|x| x.0 == starts[i].0)
                .count()
                + i;
            if j - i + 1 >= n {
                if best.is_some() {
                    bump(&COUNTS.disagree);
                    return None;
                }
                best = Some((starts[i].0, j - i));
            }
            i = j;
        }
        let Some((s, _)) = best else {
            bump(&COUNTS.too_few_pieces);
            return None;
        };
        s
    } else {
        let zero = (0..n).filter(|&k| count(k) == 0).count();
        let Some(k0) = (0..n).find(|&k| count(k) == 1).filter(|_| zero <= 1) else {
            bump(&COUNTS.too_few_pieces);
            return None;
        };
        let Some(s) = index.locate(iv[k0].0).checked_sub(cut(k0) as u64) else {
            bump(&COUNTS.disagree);
            return None;
        };
        s
    };
    let hit = index.hit_at(s, len as u64);
    if hit.is_none() {
        bump(&COUNTS.verify_failed);
        return None;
    }
    let mut mismatch = None;
    for (i, &c) in read.iter().enumerate() {
        let r = index.text_base(s + i as u64);
        if r != c {
            if mismatch.is_some() {
                bump(&COUNTS.verify_failed);
                return None;
            }
            mismatch = Some(i as i32);
        }
    }
    let bad = mismatch.map(|m| (0..n).find(|&k| (m as usize) < cut(k + 1)).unwrap_or(n - 1));
    if many {
        // Jedes Stück ohne die Abweichung muss bei `s` liegen.
        for k in 0..n {
            let at_s = starts.iter().any(|&(x, kk)| kk == k && x == s);
            if (Some(k) == bad) == at_s {
                bump(&COUNTS.disagree);
                return None;
            }
        }
    } else {
        // Jedes Stück ohne die Abweichung liegt genau bei `s` und hat dort seine eine
        // Fundstelle. Das Stück mit der Abweichung darf keine haben, sonst passt es anderswo.
        for k in 0..n {
            if (Some(k) == bad) != (count(k) == 0) {
                bump(&COUNTS.disagree);
                return None;
            }
        }
    }
    let hit = hit?;
    let sc = &opt.region.scoring;
    let l = len as i32;
    // `truesc` ist der Score über den ganzen Read. `score` ist der beste lokale Score, wie ihn
    // das Verlängern vom Seed aus findet: Liegt die Abweichung `k` Basen vor einem Ende, lohnt es
    // sich lokal erst ab k > 4 Basen (Gewinn k, Abzug 4), über sie hinwegzugehen.
    let (score, truesc) = match mismatch {
        None => (l * sc.matching, l * sc.matching),
        Some(m) => {
            let short = m.min(l - 1 - m);
            let across = (l - 1) * sc.matching - sc.mismatch;
            let stop = (l - 1 - short) * sc.matching;
            (across.max(stop), across)
        }
    };
    if many {
        // Die übrigen Anfänge sind weitere Kopien. Liegen mehrere dicht beieinander, ist es
        // dieselbe Kopie (mit Lücke); sie wird einmal ausgerichtet.
        let mut others: Vec<u64> = starts.iter().map(|x| x.0).filter(|&x| x != s).collect();
        others.dedup();
        let mut windows: Vec<(u64, u64)> = Vec::new();
        let total = index.reference().len();
        for o in others {
            if o.abs_diff(s) < len as u64 + 2 * COPY_MARGIN {
                bump(&COUNTS.copy_too_good);
                return None;
            }
            let (half_beg, half_end) = if o < total {
                (0, total)
            } else {
                (total, 2 * total)
            };
            let beg = o.saturating_sub(COPY_MARGIN).max(half_beg);
            let end = (o + len as u64 + COPY_MARGIN).min(half_end);
            match windows.last_mut() {
                Some(w) if beg <= w.1 => w.1 = w.1.max(end),
                _ => windows.push((beg, end)),
            }
        }
        if windows.len() > MAX_COPIES {
            bump(&COUNTS.repeat);
            return None;
        }
        let limit = copy_limit(score, opt);
        let byte = len as i32 * sc.matching < 250;
        let filter = opt.seed.fast_copy_filter && kmer_bound_holds(sc);
        for (beg, end) in windows {
            let target = index.text_range(beg, end);
            // Kann die Kopie die Grenze beweisbar nicht überschreiten, entfällt die Ausrichtung.
            if filter && kmer_bound(read, &target).is_some_and(|b| b <= limit) {
                continue;
            }
            let r = crate::local::local_align(read, &target, sc, byte, 0x10000, false);
            if r.score > limit {
                bump(&COUNTS.copy_too_good);
                return None;
            }
        }
        bump(&COUNTS.copies_checked);
    }
    let (seedcov, seedlen0) = seed_coverage(l, mismatch, opt.seed.min_seed_len as i32);
    bump(&COUNTS.taken);
    Some(Region {
        rb: s,
        re: s + len as u64,
        qb: 0,
        qe: l,
        contig: hit.contig,
        score,
        truesc,
        w: sc.band_width,
        seedcov,
        seedlen0,
        sub: 0,
        sub_n: 0,
        secondary: -1,
        csub: 0,
        rescued: false,
    })
}

/// Gleich in allem, was die Ausgabe bestimmt. `seedcov` und `seedlen0` zählen nicht: Sie
/// steuern nur das Verlängern innerhalb einer Kette und das Zusammenlegen von Bereichen, und
/// bei einem einzigen Bereich über den ganzen Read gehen sie in keine Zeile der Ausgabe ein
/// (die MAPQ rechnet ohne sie, siehe `finalize::mapq`).
/// Höchster Score einer weiteren Kopie, bei dem die Ausgabe gleich bleibt. Eine Kopie mit
/// Score `x` neben dem besten Score `b`:
/// - erscheint in `XA` erst ab `x >= 0,8 b`,
/// - senkt die MAPQ (60 bei 150 Basen) erst ab `b - x` unter etwa 17,
/// - wird erst ab `b - x <= 17` (Abzug für ungepaarte Reads) zum Ausgangspunkt der Partnersuche
///   und kann erst dann ein Paar bilden, das das beste Paar schlägt oder seine MAPQ senkt.
///
/// Der Abstand 31 lässt für alle drei Spielraum; der Score der Kopie im vollen Weg ist höchstens
/// der lokale Score über das ganze Fenster.
pub fn copy_limit(score: i32, opt: &AlignOptions) -> i32 {
    let xa = (f64::from(score) * f64::from(opt.sam.xa_drop_ratio)).ceil() as i32 - 1;
    xa.min(score - 31)
}

pub fn same_output(a: &Region, b: &Region) -> bool {
    let norm = |r: &Region| Region {
        seedcov: 0,
        seedlen0: 0,
        ..*r
    };
    norm(a) == norm(b)
}

/// Vergleicht im Prüfmodus den Bereich des Schnellwegs mit dem Ergebnis des vollen Wegs.
pub fn check_against_full(fast: &Region, full: &[Region], opt: &AlignOptions) {
    let limit = copy_limit(fast.score, opt);
    let same = !full.is_empty()
        && same_output(&full[0], fast)
        && full[1..].iter().all(|r| r.score <= limit);
    bump(if same {
        &COUNTS.same_as_full
    } else {
        &COUNTS.differs_from_full
    });
}

/// Stand der Zähler als (Name, Wert).
pub fn counts() -> Vec<(&'static str, u64)> {
    let c = &COUNTS;
    let g = |x: &AtomicU64| x.load(Ordering::Relaxed);
    vec![
        ("versucht", g(&c.tried)),
        ("genommen", g(&c.taken)),
        ("mehrere_fundstellen", g(&c.repeat)),
        ("zu_wenig_stuecke", g(&c.too_few_pieces)),
        ("stuecke_uneinig", g(&c.disagree)),
        ("pruefung_nicht_bestanden", g(&c.verify_failed)),
        ("kopien_geprueft", g(&c.copies_checked)),
        ("kopie_zu_gut", g(&c.copy_too_good)),
        ("gleich_wie_voller_weg", g(&c.same_as_full)),
        ("anders_als_voller_weg", g(&c.differs_from_full)),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::find_regions;

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn setup() -> (FmIndex, Vec<u8>) {
        let mut rng = 0x5eed_1234_u64;
        let s: Vec<u8> = (0..20_000)
            .map(|_| b"ACGT"[(xorshift(&mut rng) & 3) as usize])
            .collect();
        let mut r = gc_index::Reference::new();
        r.push_contig("a", &s);
        let index = FmIndex::build(r, gc_index::Variant::Kompakt, 0).unwrap().0;
        let codes = s
            .iter()
            .map(|&b| gc_index::reference::encode(b).unwrap())
            .collect();
        (index, codes)
    }

    fn revcomp(x: &[u8]) -> Vec<u8> {
        x.iter().rev().map(|&c| 3 - c).collect()
    }

    #[test]
    fn easy_reads_give_the_same_region_as_the_full_way() {
        let (index, g) = setup();
        let opt = AlignOptions::default();
        let mut rng = 99u64;
        let mut taken = 0;
        for _ in 0..300 {
            let at = (xorshift(&mut rng) % (g.len() as u64 - 200)) as usize;
            let mut read = g[at..at + 150].to_vec();
            // Ein Drittel ohne, ein Drittel mit einer Abweichung, ein Drittel mit zweien.
            let k = xorshift(&mut rng) % 3;
            for _ in 0..k {
                let p = (xorshift(&mut rng) % 150) as usize;
                read[p] = (read[p] + 1 + (xorshift(&mut rng) % 3) as u8) & 3;
            }
            if xorshift(&mut rng) & 1 == 1 {
                read = revcomp(&read);
            }
            let (full, _) = find_regions(&index, &read, &opt);
            if let Some(fast) = fast_region(&index, &read, &opt, true) {
                taken += 1;
                assert!(
                    full.len() == 1 && same_output(&full[0], &fast),
                    "Read bei {at}, {k} Abweichungen: {full:?} gegen {fast:?}"
                );
            }
        }
        assert!(taken > 150, "nur {taken} Reads auf dem Schnellweg");
    }

    #[test]
    fn a_second_copy_sends_the_read_the_full_way() {
        let (_, g) = setup();
        let mut s: Vec<u8> = g.iter().map(|&c| b"ACGT"[c as usize]).collect();
        // Die Stelle 1000..1300 noch einmal anhängen: jedes Stück hat dann zwei Fundstellen.
        let copy = s[1000..1300].to_vec();
        s.extend_from_slice(&copy);
        let mut r = gc_index::Reference::new();
        r.push_contig("a", &s);
        let index = FmIndex::build(r, gc_index::Variant::Kompakt, 0).unwrap().0;
        let read = g[1050..1200].to_vec();
        assert!(fast_region(&index, &read, &AlignOptions::default(), true).is_none());
    }

    #[test]
    fn a_weak_copy_keeps_the_read_on_the_fast_way() {
        let (_, g) = setup();
        let mut s: Vec<u8> = g.iter().map(|&c| b"ACGT"[c as usize]).collect();
        // Die Stelle 1000..1300 noch einmal anhängen, aber mit einer Abweichung alle 9 Basen:
        // Stücke von 21 Basen passen dort nirgends genau. Dazu 60 Basen von 5100 genau: Dort
        // passen zwei Stücke eines Reads ab 5050 genau, der Score bleibt aber bei etwa 60.
        let mut weak = s[1000..1300].to_vec();
        for (i, b) in weak.iter_mut().enumerate() {
            if i % 9 == 4 {
                *b = if *b == b'A' { b'C' } else { b'A' };
            }
        }
        let half = s[5100..5160].to_vec();
        s.extend_from_slice(&weak);
        s.extend_from_slice(&half);
        let mut r = gc_index::Reference::new();
        r.push_contig("a", &s);
        let index = FmIndex::build(r, gc_index::Variant::Kompakt, 0).unwrap().0;
        let opt = AlignOptions::default();
        for at in [1050, 5050] {
            let read = g[at..at + 150].to_vec();
            let fast = fast_region(&index, &read, &opt, true).expect("Schnellweg bei {at}");
            let (full, _) = find_regions(&index, &read, &opt);
            assert!(same_output(&full[0], &fast), "{full:?} gegen {fast:?}");
            let limit = copy_limit(fast.score, &opt);
            assert!(full[1..].iter().all(|r| r.score <= limit), "{full:?}");
            // Ohne Prüfung der Kopien geht der Read bei 5050 den vollen Weg.
            assert_eq!(
                fast_region(&index, &read, &opt, false).is_some(),
                at == 1050
            );
        }
    }

    #[test]
    fn reads_with_n_or_two_mismatches_go_the_full_way() {
        let (index, g) = setup();
        let opt = AlignOptions::default();
        let mut read = g[500..650].to_vec();
        read[10] = 4;
        assert!(fast_region(&index, &read, &opt, true).is_none());
        let mut read = g[500..650].to_vec();
        read[60] = (read[60] + 1) & 3;
        read[140] = (read[140] + 1) & 3;
        assert!(fast_region(&index, &read, &opt, true).is_none());
    }
    #[test]
    fn kmer_bound_is_never_below_the_local_score() {
        let mut rng = 0x5151_2626_u64;
        let sc = crate::extend::Scoring::default();
        assert!(kmer_bound_holds(&sc));
        let mut tight = 0;
        for round in 0..3000 {
            let n = 30 + (xorshift(&mut rng) % 150) as usize;
            let read: Vec<u8> = (0..n).map(|_| (xorshift(&mut rng) & 3) as u8).collect();
            // Fenster: Teile des Reads mit Fehlern, Lücken, N und Zufall dazwischen.
            let mut target: Vec<u8> = (0..xorshift(&mut rng) % 60)
                .map(|_| (xorshift(&mut rng) & 3) as u8)
                .collect();
            let mut i = (xorshift(&mut rng) as usize) % n;
            while i < n {
                let b = read[i];
                match xorshift(&mut rng) % [8, 20, 50][round % 3] {
                    0 => target.push((b + 1) % 4),
                    1 => {}
                    2 => {
                        target.push(b);
                        target.push((xorshift(&mut rng) & 3) as u8);
                    }
                    3 => target.push(4),
                    _ => target.push(b),
                }
                i += 1;
            }
            target.extend((0..xorshift(&mut rng) % 60).map(|_| (xorshift(&mut rng) & 3) as u8));
            let Some(bound) = kmer_bound(&read, &target) else {
                assert!(target.contains(&4));
                continue;
            };
            for byte in [false, true] {
                let byte = byte && n < 250;
                let r = crate::local::local_align(&read, &target, &sc, byte, 0x10000, false);
                assert!(
                    bound >= r.score,
                    "Runde {round}: Grenze {bound} < {}",
                    r.score
                );
                if bound == r.score {
                    tight += 1;
                }
            }
        }
        assert!(tight > 100, "{tight}");
    }
}
