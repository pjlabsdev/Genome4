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
    /// An: leichte Reads nehmen den Schnellweg.
    On,
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

/// Sucht den einen Bereich eines leichten Reads. `None` heißt: voller Weg.
pub fn fast_region(index: &FmIndex, read: &[u8], opt: &AlignOptions) -> Option<Region> {
    bump(&COUNTS.tried);
    let len = read.len();
    if len < 60 || read.iter().any(|&c| c > 3) {
        bump(&COUNTS.too_few_pieces);
        return None;
    }
    // Stücke von mindestens PIECE_MIN Basen. Eine zweite Kopie, die sich in höchstens
    // `n - 1` Ereignissen (Abweichung oder Lücke) vom Read unterscheidet, passt in mindestens
    // einem Stück genau und fällt dort als zweite Fundstelle auf.
    let n = (len / PIECE_MIN).max(3);
    let cut = |k: usize| k * len / n;
    let mut counts = [0u64; 16];
    let n = n.min(counts.len());
    let mut first = None;
    for (k, c) in counts.iter_mut().enumerate().take(n) {
        let (lo, hi) = index.interval(&read[cut(k)..cut(k + 1)]);
        *c = hi - lo;
        if *c > 1 {
            bump(&COUNTS.repeat);
            return None;
        }
        if *c == 1 && first.is_none() {
            first = Some((k, lo));
        }
    }
    let zero = counts[..n].iter().filter(|&&c| c == 0).count();
    let Some((k0, row)) = first.filter(|_| zero <= 1) else {
        bump(&COUNTS.too_few_pieces);
        return None;
    };
    let Some(s) = index.locate(row).checked_sub(cut(k0) as u64) else {
        bump(&COUNTS.disagree);
        return None;
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
    // Jedes Stück ohne die Abweichung liegt genau bei `s` und hat dort seine eine Fundstelle.
    // Das Stück mit der Abweichung darf keine haben, sonst passt es anderswo genau.
    let bad = mismatch.map(|m| (0..n).find(|&k| (m as usize) < cut(k + 1)).unwrap_or(n - 1));
    for (k, &c) in counts.iter().enumerate().take(n) {
        if (Some(k) == bad) != (c == 0) {
            bump(&COUNTS.disagree);
            return None;
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
pub fn same_output(a: &Region, b: &Region) -> bool {
    let norm = |r: &Region| Region {
        seedcov: 0,
        seedlen0: 0,
        ..*r
    };
    norm(a) == norm(b)
}

/// Vergleicht im Prüfmodus den Bereich des Schnellwegs mit dem Ergebnis des vollen Wegs.
pub fn check_against_full(fast: &Region, full: &[Region]) {
    let same = full.len() == 1 && same_output(&full[0], fast);
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
            if let Some(fast) = fast_region(&index, &read, &opt) {
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
        assert!(fast_region(&index, &read, &AlignOptions::default()).is_none());
    }

    #[test]
    fn reads_with_n_or_two_mismatches_go_the_full_way() {
        let (index, g) = setup();
        let opt = AlignOptions::default();
        let mut read = g[500..650].to_vec();
        read[10] = 4;
        assert!(fast_region(&index, &read, &opt).is_none());
        let mut read = g[500..650].to_vec();
        read[60] = (read[60] + 1) & 3;
        read[140] = (read[140] + 1) & 3;
        assert!(fast_region(&index, &read, &opt).is_none());
    }
}
