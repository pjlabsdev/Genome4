//! Strobemer-Weg: Reads über Randstrobes verorten statt über Seed-Suche, Ketten und Verlängern.
//!
//! Für jeden Read, den der Schnellweg nicht nimmt: Randstrobes beider Orientierungen im
//! [`StrobeIndex`] nachschlagen, Treffer nach Diagonale gruppieren und die größten Gruppen lokal
//! gegen die Referenz ausrichten. Genommen wird der Read nur, wenn
//! - kein Randstrobe wegen zu vieler Orte gesperrt war (sonst könnte der wahre Ort in einer
//!   Wiederholung verborgen sein),
//! - der beste Ort den ganzen Read ohne Abschneiden ausrichtet, mit höchstens `max_loss` Punkten
//!   unter dem vollen Score,
//! - der zweitbeste Ort mindestens `min_gap` Punkte schlechter ist.
//!
//! Dann entsteht genau ein Bereich wie beim vollen Weg; Partnersuche, Paarbewertung, MAPQ und
//! Zeilen laufen unverändert. Alle anderen Reads gehen den vollen Weg. Die Ausgabe ist damit
//! nicht mehr bwa-gleich und wird gegen die Wahrheit geprüft (Simulation, GIAB).

use crate::region::Region;
use crate::sam::AlignOptions;
use gc_index::strobe::{randstrobes, syncmers};
use gc_index::{FmIndex, StrobeIndex};
use std::sync::atomic::{AtomicU64, Ordering};

/// Einstellungen.
#[derive(Debug, Clone, PartialEq)]
pub struct StrobeOptions {
    /// Randstrobes mit mehr Orten gelten als gesperrt.
    pub max_occ: usize,
    /// So viele der größten Gruppen werden ausgerichtet.
    pub candidates: usize,
    /// Höchstens so viele Punkte unter dem vollen Score (Readlänge mal Punkte je Base).
    pub max_loss: i32,
    /// Der zweitbeste Ort muss mindestens so viele Punkte schlechter sein.
    pub min_gap: i32,
    /// Rand des Referenzfensters links und rechts (für Lücken).
    pub margin: u64,
    /// Höchstens so viel Prozent der Randstrobes beider Orientierungen dürfen gesperrt sein.
    pub max_blocked_percent: usize,
    /// Erst ohne Lücke entlang der Diagonale vergleichen, lokale Ausrichtung nur bei
    /// Lückenverdacht (Schalter zurück: `GC_STROBE_HAMMING=0`).
    pub ungapped_first: bool,
    /// Sobald ein Zweiter zu nah ist, abbrechen (Schalter zurück: `GC_STROBE_FRUEH=0`).
    pub stop_early: bool,
    /// Ohne Lücke nur, wenn keine kurze Lücke besser passt (Schalter zurück:
    /// `GC_STROBE_LUECKENPRUEFUNG=0`).
    pub gap_check: bool,
    /// Versuch: Reads mit zu nahem Zweitem selbst nehmen (Zweiter als `sub` und `csub`)
    /// statt den vollen Weg (`GC_STROBE_NAH_NEHMEN=1`).
    pub take_close: bool,
    /// Dazu alle fast gleich guten Orte zurückgeben (`GC_STROBE_NAH_NEHMEN=2`).
    pub take_close_multi: bool,
}

impl Default for StrobeOptions {
    fn default() -> Self {
        Self {
            max_occ: 1000,
            candidates: 16,
            max_loss: 30,
            min_gap: std::env::var("GC_STROBE_ABSTAND")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(15),
            margin: 12,
            max_blocked_percent: std::env::var("GC_STROBE_GESPERRT")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(25),
            ungapped_first: std::env::var("GC_STROBE_HAMMING").map_or(true, |v| v != "0"),
            stop_early: std::env::var("GC_STROBE_FRUEH").map_or(true, |v| v != "0"),
            gap_check: std::env::var("GC_STROBE_LUECKENPRUEFUNG").map_or(true, |v| v != "0"),
            take_close: std::env::var("GC_STROBE_NAH_NEHMEN").is_ok_and(|v| v == "1" || v == "2"),
            take_close_multi: std::env::var("GC_STROBE_NAH_NEHMEN").is_ok_and(|v| v == "2"),
        }
    }
}

pub struct StrobeCounts {
    pub tried: AtomicU64,
    pub taken: AtomicU64,
    pub blocked: AtomicU64,
    pub no_hit: AtomicU64,
    pub clipped: AtomicU64,
    pub too_weak: AtomicU64,
    pub second_close: AtomicU64,
    pub ungapped: AtomicU64,
    pub local: AtomicU64,
}

pub static COUNTS: StrobeCounts = StrobeCounts {
    tried: AtomicU64::new(0),
    taken: AtomicU64::new(0),
    blocked: AtomicU64::new(0),
    no_hit: AtomicU64::new(0),
    clipped: AtomicU64::new(0),
    too_weak: AtomicU64::new(0),
    second_close: AtomicU64::new(0),
    ungapped: AtomicU64::new(0),
    local: AtomicU64::new(0),
};

fn bump(c: &AtomicU64) {
    c.fetch_add(1, Ordering::Relaxed);
}

/// Zähler für das Protokoll.
pub fn counts() -> Vec<(&'static str, u64)> {
    let c = &COUNTS;
    vec![
        ("versucht", c.tried.load(Ordering::Relaxed)),
        ("genommen", c.taken.load(Ordering::Relaxed)),
        ("gesperrt", c.blocked.load(Ordering::Relaxed)),
        ("ohne_treffer", c.no_hit.load(Ordering::Relaxed)),
        ("abgeschnitten", c.clipped.load(Ordering::Relaxed)),
        ("zu_schwach", c.too_weak.load(Ordering::Relaxed)),
        ("zweiter_zu_nah", c.second_close.load(Ordering::Relaxed)),
        ("ohne_luecke", c.ungapped.load(Ordering::Relaxed)),
        ("lokal", c.local.load(Ordering::Relaxed)),
        (
            "paar_versucht",
            crate::paarweg::TRIED.load(Ordering::Relaxed),
        ),
        (
            "paar_genommen",
            crate::paarweg::TAKEN.load(Ordering::Relaxed),
        ),
        (
            "paar_unsicher",
            crate::paarweg::UNSURE.load(Ordering::Relaxed),
        ),
    ]
}

thread_local! {
    /// Der letzte Aufruf endete mit zu nahem Zweitem (für den Versuch Strobe-Saat).
    pub static LAST_CLOSE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static SCRATCH: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::default());
}

#[derive(Default)]
struct Scratch {
    sy: Vec<(u32, u64)>,
    sc: Vec<u64>,
    st: Vec<(u64, u32)>,
    rc: Vec<u8>,
    hits: Vec<u64>,
    groups: Vec<(u32, u64)>,
    gap: Vec<i32>,
}

/// Größte geprüfte Lücke beim Vergleich ohne Lücke (Schalter zurück: `GC_STROBE_LUECKE_BIS=3`).
fn max_gap() -> usize {
    static G: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *G.get_or_init(|| {
        std::env::var("GC_STROBE_LUECKE_BIS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(15)
    })
}

/// Größte Lückenlänge `d <= max`, mit der `n` Basen überhaupt mehr als `full` Punkte erreichen
/// können (Deletion höchstens `n` Treffer minus Abzug, Insertion `n - d`); 0, wenn keine. Längere
/// Lücken brauchen nicht geprüft zu werden, die Antwort von [`best_one_gap`] bleibt dieselbe.
fn gap_reach(n: usize, full: i32, sc: &crate::extend::Scoring, max: usize) -> usize {
    let all = n as i32 * sc.matching;
    (1..=max.min(n.saturating_sub(1)))
        .take_while(|&d| {
            let d32 = d as i32;
            let del = all - (sc.o_del + sc.e_del * d32);
            let ins = all - d32 * sc.matching - (sc.o_ins + sc.e_ins * d32);
            del.max(ins) > full
        })
        .last()
        .unwrap_or(0)
}

/// Welche Seite eines Stücks auf der gegebenen Diagonale bleiben muss.
#[derive(Clone, Copy, PartialEq)]
enum Fixed {
    /// Anfang (rechtes Ende eines Reads hinter einer Ausrichtung)
    Start,
    /// Ende (linkes Ende eines Reads vor einer Ausrichtung)
    End,
    /// eine von beiden (ganzer Read)
    Either,
}

/// Bester Score von `q` ganz, mit genau einer Lücke von 1 bis `g` Basen, wobei die
/// feste Seite ohne Lücke auf `t[off + i]` liegt. `None`, wenn das Fenster dafür nicht reicht.
fn best_one_gap(
    q: &[u8],
    t: &[u8],
    off: usize,
    g: usize,
    sc: &crate::extend::Scoring,
    fixed: Fixed,
    buf: &mut Vec<i32>,
) -> Option<i32> {
    let n = q.len();
    if off < g || off + n + g > t.len() || n < 2 {
        return None;
    }
    // Für jede Diagonale x = off - g ..= off + g: Vorsummen P[k] und Nachsummen S[k]
    let w = n + 1;
    let nd = 2 * g + 1;
    buf.clear();
    buf.resize(2 * nd * w, 0);
    let (pre, suf) = buf.split_at_mut(nd * w);
    for j in 0..nd {
        let x = off + j - g;
        let (p, s) = (&mut pre[j * w..(j + 1) * w], &mut suf[j * w..(j + 1) * w]);
        for i in 0..n {
            let v = if q[i] == t[x + i] {
                sc.matching
            } else {
                -sc.mismatch
            };
            p[i + 1] = p[i] + v;
        }
        for i in (0..n).rev() {
            let v = if q[i] == t[x + i] {
                sc.matching
            } else {
                -sc.mismatch
            };
            s[i] = s[i + 1] + v;
        }
    }
    let p = |j: usize, k: usize| pre[j * w + k];
    let s = |j: usize, k: usize| suf[j * w + k];
    let mut best = i32::MIN;
    for d in 1..=g {
        let del = sc.o_del + sc.e_del * d as i32;
        let ins = sc.o_ins + sc.e_ins * d as i32;
        // (linke Diagonale, rechte Diagonale) als Index in 0..nd; Mitte g = `off`
        let all = [
            (g, g + d, true, fixed != Fixed::End), // Deletion: rechts um d weiter
            (g, g - d, false, fixed != Fixed::End), // Insertion: rechts um d zurück
            (g - d, g, true, fixed != Fixed::Start),
            (g + d, g, false, fixed != Fixed::Start),
        ];
        for (lj, rj, deletion, _) in all.into_iter().filter(|a| a.3) {
            let (lo, hi) = match fixed {
                Fixed::Start => (0, n),
                Fixed::End => (0, n),
                Fixed::Either => (1, n - 1),
            };
            if deletion {
                for k in lo..=hi {
                    best = best.max(p(lj, k) + s(rj, k) - del);
                }
            } else {
                for k in lo..=hi.saturating_sub(d) {
                    if k + d <= n {
                        best = best.max(p(lj, k) + s(rj, k + d) - ins);
                    }
                }
            }
        }
    }
    Some(best)
}

/// Ein Bereich für den Read oder `None` (dann der volle Weg).
pub fn strobe_region(
    index: &FmIndex,
    strobes: &StrobeIndex,
    read: &[u8],
    opt: &AlignOptions,
) -> Option<Region> {
    strobe_regions(index, strobes, read, opt).map(|mut v| v.swap_remove(0))
}

/// Wie [`strobe_region`], im Versuch `GC_STROBE_NAH_NEHMEN=2` bei nahem Zweitem aber alle
/// fast gleich guten Orte (der beste zuerst), damit die Paarung unter ihnen wählen kann.
pub fn strobe_regions(
    index: &FmIndex,
    strobes: &StrobeIndex,
    read: &[u8],
    opt: &AlignOptions,
) -> Option<Vec<Region>> {
    bump(&COUNTS.tried);
    LAST_CLOSE.with(|c| c.set(false));
    let so = &opt.strobe;
    let len = read.len();
    if len < 60 || read.iter().any(|&c| c > 3) {
        bump(&COUNTS.no_hit);
        return None;
    }
    SCRATCH.with(|s| {
        let s = &mut *s.borrow_mut();
        s.rc.clear();
        s.rc.extend(read.iter().rev().map(|&c| 3 - c));
        s.hits.clear();
        let (mut blocked, mut usable) = (0usize, 0usize);
        // Treffer: Strang (Bit 63) | Diagonale + 2^32 (Ort im Vorwärtsstrang minus Ort im Read)
        for strand in 0..2u64 {
            let q: &[u8] = if strand == 0 { read } else { &s.rc };
            syncmers(q, &mut s.sy, &mut s.sc);
            randstrobes(&s.sy, &mut s.st);
            for &(h, qp) in &s.st {
                match strobes.lookup(h, so.max_occ) {
                    None => {
                        blocked += 1;
                    }
                    Some(it) => {
                        usable += 1;
                        for p in it {
                            let d = (p + (1 << 32)).wrapping_sub(u64::from(qp));
                            s.hits.push(strand << 63 | d);
                        }
                    }
                }
            }
        }
        // Gesperrte Randstrobes (Wiederholungen) könnten einen gleich guten Ort verbergen. Erlaubt
        // ist nur ein kleiner Anteil, damit der Rest des Reads jede nahe Kopie mitfindet.
        if blocked > 0 && blocked * 100 > (blocked + usable) * so.max_blocked_percent {
            bump(&COUNTS.blocked);
            return None;
        }
        if s.hits.is_empty() {
            bump(&COUNTS.no_hit);
            return None;
        }
        s.hits.sort_unstable();
        s.groups.clear();
        let mut i = 0;
        while i < s.hits.len() {
            let mut j = i + 1;
            while j < s.hits.len() && s.hits[j] - s.hits[j - 1] <= 10 {
                j += 1;
            }
            s.groups.push(((j - i) as u32, s.hits[(i + j) / 2]));
            i = j;
        }
        s.groups
            .sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        let sc = &opt.region.scoring;
        let l = index.reference().len();
        let byte = len as i32 * sc.matching < 250;
        let bound_ok = crate::fast::kmer_bound_holds(sc);
        // (Score, Textanfang, Textende, Read-Anfang, Read-Ende)
        let mut best: Option<(i32, u64, u64, i32, i32)> = None;
        let mut all: Vec<(i32, u64, u64, i32, i32)> = Vec::new();
        let mut second = 0i32;
        for &(_, g) in s.groups.iter().take(so.candidates) {
            let strand = g >> 63;
            let diag = (g & !(1 << 63)) as i64 - (1 << 32);
            // Fenster im Vorwärtsstrang, dann im Indextext (Gegenstrang: gespiegelt)
            let fb = (diag - so.margin as i64).max(0) as u64;
            let fe = ((diag + len as i64 + so.margin as i64).max(0) as u64).min(l);
            if fe <= fb + len as u64 / 2 {
                continue;
            }
            let (tb, te) = if strand == 0 {
                (fb, fe)
            } else {
                (2 * l - fe, 2 * l - fb)
            };
            let target = index.text_range(tb, te);
            if target.iter().any(|&c| c > 3) {
                continue;
            }
            // Beweisbar zu schwache Kandidaten nicht ausrichten: ihre Schranke zählt als Zweiter
            // (vorsichtig, die MAPQ wird dadurch höchstens kleiner).
            if let Some((b, ..)) = best
                && bound_ok
                && let Some(bound) = crate::fast::kmer_bound(read, &target)
                && bound <= b - so.min_gap
            {
                second = second.max(bound);
                continue;
            }
            // Ohne Lücke entlang der Diagonale: bester lokaler Abschnitt (Kadane). Reicht er für
            // die Annahme und liegt die Schranke nicht darüber, entfällt die lokale Ausrichtung.
            let off = if strand == 0 {
                diag - fb as i64
            } else {
                fe as i64 - diag - len as i64
            };
            let ungapped = if so.ungapped_first && off >= 0 && off as usize + len <= target.len() {
                let t = &target[off as usize..off as usize + len];
                let (mut run, mut run_b) = (0i32, 0usize);
                let mut top = (0i32, 0usize, 0usize);
                for (i, (a, b)) in read.iter().zip(t).enumerate() {
                    run += if a == b { sc.matching } else { -sc.mismatch };
                    if run <= 0 {
                        run = 0;
                        run_b = i + 1;
                    } else if run > top.0 {
                        top = (run, run_b, i + 1);
                    }
                }
                // Die Enden außerhalb des Abschnitts müssen sich ohne Lücke anhängen lassen
                // (wie unten beim Verlängern), sonst droht dort eine Lücke.
                let tail = |q: &[u8], t: &[u8]| -> i32 {
                    q.iter()
                        .zip(t)
                        .map(|(a, b)| if a == b { sc.matching } else { -sc.mismatch })
                        .sum()
                };
                let ends_ok = tail(&read[..top.1], &t[..top.1]) >= -opt.region.pen_clip5
                    && tail(&read[top.2..], &t[top.2..]) >= -opt.region.pen_clip3;
                // Eine kurze Lücke darf nicht besser passen als ohne Lücke (Indel nahe am Ende),
                // sonst richtet die lokale Ausrichtung sie aus.
                let full = tail(read, t);
                let g = gap_reach(len, full, sc, max_gap());
                let no_gap_better = !so.gap_check
                    || g == 0
                    || best_one_gap(
                        read,
                        &target,
                        off as usize,
                        g,
                        sc,
                        Fixed::Either,
                        &mut s.gap,
                    )
                    .is_some_and(|b| b <= full);
                let tight = ends_ok
                    && no_gap_better
                    && (!bound_ok
                        || crate::fast::kmer_bound(read, &target).is_none_or(|b| b <= top.0)
                        || best.is_none_or(|(b, ..)| top.0 > b));
                (top.0 >= len as i32 * sc.matching - so.max_loss && tight).then(|| {
                    let o = tb + off as u64;
                    (
                        top.0,
                        o + top.1 as u64,
                        o + top.2 as u64,
                        top.1 as i32,
                        top.2 as i32,
                    )
                })
            } else {
                None
            };
            let cand = if let Some(c) = ungapped {
                bump(&COUNTS.ungapped);
                c
            } else {
                bump(&COUNTS.local);
                let r = crate::local::local_align(read, &target, sc, byte, 0, true);
                (
                    r.score,
                    tb + r.tb as u64,
                    tb + r.te as u64 + 1,
                    r.qb,
                    r.qe + 1,
                )
            };
            if so.take_close_multi {
                all.push(cand);
            }
            match best {
                Some(b) if cand.0 <= b.0 => {
                    // dieselbe Stelle (überlappend) zählt nicht als zweiter Ort
                    if cand.1 >= b.2 || cand.2 <= b.1 {
                        second = second.max(cand.0);
                    }
                }
                _ => {
                    if let Some(b) = best
                        && (b.1 >= cand.2 || b.2 <= cand.1)
                    {
                        second = second.max(b.0);
                    }
                    best = Some(cand);
                }
            }
            // Ein zu naher Zweiter schickt den Read ohnehin den vollen Weg: dann gleich aufhören
            // statt weitere Kandidaten auszurichten (vorsichtig, ein später klar Bester ginge
            // dadurch nur den vollen Weg).
            if so.stop_early
                && !so.take_close
                && let Some((b, ..)) = best
                && b - second < so.min_gap
            {
                bump(&COUNTS.second_close);
                LAST_CLOSE.with(|c| c.set(blocked == 0));
                return None;
            }
        }
        let Some((score, mut rb, mut re, qb, qe)) = best else {
            bump(&COUNTS.no_hit);
            return None;
        };
        // Abgeschnittene Enden ohne Lücke bis zum Read-Ende verlängern, wenn das je Seite höchstens
        // den Abzug fürs Abschneiden kostet (wie bwa mit `-L 5,5`); sonst der volle Weg.
        let mut truesc = score;
        let side = |q: &[u8], t: &[u8]| -> i32 {
            q.iter()
                .zip(t)
                .map(|(a, b)| if a == b { sc.matching } else { -sc.mismatch })
                .sum()
        };
        if qb > 0 {
            let n = qb as u64;
            if rb < n || (rb < l) != (rb - n < l) {
                bump(&COUNTS.clipped);
                return None;
            }
            let d = side(&read[..qb as usize], &index.text_range(rb - n, rb));
            let gu = gap_reach(n as usize, d, sc, max_gap());
            let g = gu as u64;
            let gapped = if !so.gap_check || g == 0 {
                Some(i32::MIN)
            } else if rb < n + g || (rb < l) != (rb - n - g < l) || (rb + g <= l) != (rb <= l) {
                None
            } else {
                let w = index.text_range(rb - n - g, rb + g);
                best_one_gap(&read[..qb as usize], &w, gu, gu, sc, Fixed::End, &mut s.gap)
            };
            if d < -opt.region.pen_clip5 || gapped.is_none_or(|b| b > d) {
                bump(&COUNTS.clipped);
                return None;
            }
            truesc += d;
            rb -= n;
        }
        if qe < len as i32 {
            let n = (len as i32 - qe) as u64;
            if re + n > 2 * l || (re <= l) != (re + n <= l) {
                bump(&COUNTS.clipped);
                return None;
            }
            let d = side(&read[qe as usize..], &index.text_range(re, re + n));
            let gu = gap_reach(n as usize, d, sc, max_gap());
            let g = gu as u64;
            let gapped = if !so.gap_check || g == 0 {
                Some(i32::MIN)
            } else if re < g
                || re + n + g > 2 * l
                || (re - g < l) != (re < l)
                || (re <= l) != (re + n + g <= l)
            {
                None
            } else {
                let w = index.text_range(re - g, re + n + g);
                best_one_gap(
                    &read[qe as usize..],
                    &w,
                    gu,
                    gu,
                    sc,
                    Fixed::Start,
                    &mut s.gap,
                )
            };
            if d < -opt.region.pen_clip3 || gapped.is_none_or(|b| b > d) {
                bump(&COUNTS.clipped);
                return None;
            }
            truesc += d;
            re += n;
        }
        if truesc < len as i32 * sc.matching - so.max_loss {
            bump(&COUNTS.too_weak);
            return None;
        }
        let close = score - second < so.min_gap;
        // Versuch mehrere Orte: Nur wenn der Strobeweg alle Kopien sehen konnte (nichts
        // gesperrt, nicht mehr Gruppen als Kandidaten), sonst der volle Weg.
        let unseen = so.take_close_multi
            && match std::env::var("GC_STROBE_NAH_VORSICHT").as_deref() {
                Ok("0") => false,
                Ok("gesperrt") => blocked > 0,
                Ok("gruppen") => s.groups.get(so.candidates).is_some_and(|g| g.0 >= 2),
                Ok("gruppen3") => s.groups.get(so.candidates).is_some_and(|g| g.0 >= 3),
                _ => blocked > 0 || s.groups.len() > so.candidates,
            };
        // Variante: Reads mit gesperrten Randstrobes trotzdem nehmen, aber ohne Sicherheit
        // (`csub` = Score: MAPQ 0, und ein Partner darauf stützt das Paar nicht).
        let zero = so.take_close_multi
            && close
            && blocked > 0
            && std::env::var("GC_STROBE_NAH_GESPERRT_NULL").is_ok_and(|v| v == "1");
        let unseen = unseen && !zero;
        if close && (!so.take_close || unseen) {
            bump(&COUNTS.second_close);
            LAST_CLOSE.with(|c| c.set(blocked == 0));
            return None;
        }
        // Contig des Bereichs; über eine Contig-Grenze geht der volle Weg.
        let hit = index.hit_at(rb, re - rb)?;
        bump(&COUNTS.taken);
        let min_sub = opt.finalize.min_seed_len * sc.matching;
        let multi = close && so.take_close_multi;
        let mut out = vec![Region {
            rb,
            re,
            qb: 0,
            qe: len as i32,
            contig: hit.contig,
            score,
            truesc,
            w: sc.band_width,
            seedcov: len as i32,
            seedlen0: len as i32,
            sub: if multi {
                0
            } else if second > min_sub {
                second
            } else {
                0
            },
            sub_n: 0,
            secondary: -1,
            csub: if zero {
                score
            } else if close && !multi {
                if std::env::var("GC_STROBE_NAH_CSUB").is_ok_and(|v| v == "voll") {
                    score
                } else {
                    second
                }
            } else {
                0
            },
            rescued: false,
        }];
        if multi {
            all.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
            for &(csc, cb, ce, cqb, cqe) in &all {
                if csc < score - opt.pair.pen_unpaired || csc < min_sub {
                    break;
                }
                if out.iter().any(|r| cb < r.re && ce > r.rb) {
                    continue;
                }
                let Some(h) = index.hit_at(cb, ce - cb) else {
                    continue;
                };
                out.push(Region {
                    rb: cb,
                    re: ce,
                    qb: cqb,
                    qe: cqe,
                    contig: h.contig,
                    score: csc,
                    truesc: csc,
                    w: sc.band_width,
                    seedcov: cqe - cqb,
                    seedlen0: cqe - cqb,
                    sub: 0,
                    sub_n: 0,
                    secondary: -1,
                    csub: if zero { csc } else { 0 },
                    rescued: false,
                });
            }
        }
        Some(out)
    })
}

/// Versuch Strobe-Saat (`GC_STROBE_SAAT=1`): Seeds für den vollen Weg aus den Strobe-Treffern
/// statt aus der SMEM-Suche. Je Diagonale (höchstens `SAAT_MAX` Gruppen) werden Read und Text
/// Base für Base verglichen; jede Strecke gleicher Basen ab `min_seed_len` wird ein Seed.
pub fn strobe_seeds(
    index: &FmIndex,
    strobes: &StrobeIndex,
    read: &[u8],
    opt: &AlignOptions,
) -> Option<Vec<crate::seed::Seed>> {
    let so = &opt.strobe;
    let len = read.len();
    let min_len = opt.seed.min_seed_len.max(1) as usize;
    SCRATCH.with(|s| {
        let s = &mut *s.borrow_mut();
        s.rc.clear();
        s.rc.extend(read.iter().rev().map(|&c| 3 - c));
        s.hits.clear();
        for strand in 0..2u64 {
            let q: &[u8] = if strand == 0 { read } else { &s.rc };
            syncmers(q, &mut s.sy, &mut s.sc);
            randstrobes(&s.sy, &mut s.st);
            for &(h, qp) in &s.st {
                let it = strobes.lookup(h, so.max_occ)?;
                for p in it {
                    let d = (p + (1 << 32)).wrapping_sub(u64::from(qp));
                    s.hits.push(strand << 63 | d);
                }
            }
        }
        s.hits.sort_unstable();
        s.groups.clear();
        let mut i = 0;
        while i < s.hits.len() {
            let mut j = i + 1;
            while j < s.hits.len() && s.hits[j] - s.hits[j - 1] <= 10 {
                j += 1;
            }
            s.groups.push(((j - i) as u32, s.hits[(i + j) / 2]));
            i = j;
        }
        s.groups
            .sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        let l = index.reference().len();
        let n_groups = s.groups.len() as u64;
        let mut seeds = Vec::new();
        for &(_, g) in s.groups.iter().take(SAAT_MAX) {
            let strand = g >> 63;
            let diag = (g & !(1 << 63)) as i64 - (1 << 32);
            let fb = diag.max(0) as u64;
            let fe = ((diag + len as i64).max(0) as u64).min(l);
            if fe <= fb + min_len as u64 {
                continue;
            }
            let (tb, te) = if strand == 0 {
                (fb, fe)
            } else {
                (2 * l - fe, 2 * l - fb)
            };
            let target = index.text_range(tb, te);
            // Read-Position des Textanfangs
            let q0 = if strand == 0 {
                fb as i64 - diag
            } else {
                diag + len as i64 - fe as i64
            };
            let mut run = 0usize;
            for (k, &t) in target.iter().enumerate() {
                let qi = q0 + k as i64;
                let same = qi >= 0 && (qi as usize) < len && read[qi as usize] == t && t < 4;
                if same {
                    run += 1;
                }
                if (!same || k + 1 == target.len()) && run >= min_len {
                    let end = if same { k + 1 } else { k };
                    let rb = tb + (end - run) as u64;
                    if let Some(h) = index.hit_at(rb, run as u64) {
                        seeds.push(crate::seed::Seed {
                            rbeg: rb,
                            qbeg: (q0 + (end - run) as i64) as u32,
                            len: run as u32,
                            contig: h.contig,
                            occurrences: n_groups,
                        });
                    }
                }
                if !same {
                    run = 0;
                }
            }
        }
        seeds.sort_unstable_by_key(|s| (s.qbeg, s.rbeg, s.len));
        seeds.dedup_by(|a, b| a.qbeg == b.qbeg && a.rbeg == b.rbeg && a.len == b.len);
        (!seeds.is_empty()).then_some(seeds)
    })
}

/// Alle Anker eines Reads (Zahl der Treffer, Strang in Bit 63 | Diagonale + 2^32), die größten
/// zuerst, dazu die Zahl gesperrter und nutzbarer Randstrobes. Für den Strobe-Paarweg.
pub fn strobe_groups(
    strobes: &StrobeIndex,
    read: &[u8],
    so: &StrobeOptions,
) -> (Vec<(u32, u64)>, usize, usize) {
    SCRATCH.with(|s| {
        let s = &mut *s.borrow_mut();
        s.rc.clear();
        s.rc.extend(read.iter().rev().map(|&c| 3 - c));
        s.hits.clear();
        let (mut blocked, mut usable) = (0usize, 0usize);
        for strand in 0..2u64 {
            let q: &[u8] = if strand == 0 { read } else { &s.rc };
            syncmers(q, &mut s.sy, &mut s.sc);
            randstrobes(&s.sy, &mut s.st);
            for &(h, qp) in &s.st {
                match strobes.lookup(h, so.max_occ) {
                    None => blocked += 1,
                    Some(it) => {
                        usable += 1;
                        for p in it {
                            let d = (p + (1 << 32)).wrapping_sub(u64::from(qp));
                            s.hits.push(strand << 63 | d);
                        }
                    }
                }
            }
        }
        s.hits.sort_unstable();
        let mut groups = Vec::new();
        let mut i = 0;
        while i < s.hits.len() {
            let mut j = i + 1;
            while j < s.hits.len() && s.hits[j] - s.hits[j - 1] <= 10 {
                j += 1;
            }
            groups.push(((j - i) as u32, s.hits[(i + j) / 2]));
            i = j;
        }
        groups.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        (groups, blocked, usable)
    })
}

/// Höchstzahl Diagonalen für die Strobe-Saat.
const SAAT_MAX: usize = 100;

#[cfg(test)]
mod tests {
    use super::*;

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn a_short_gap_near_the_end_beats_the_ungapped_comparison() {
        let sc = crate::extend::Scoring::default();
        let mut state = 0x5eed_u64;
        let t: Vec<u8> = (0..200).map(|_| (xorshift(&mut state) & 3) as u8).collect();
        let mut buf = Vec::new();
        let off = 20;
        let n = 150;
        // ohne Änderung: keine Lücke ist besser als ohne Lücke
        let q: Vec<u8> = t[off..off + n].to_vec();
        let full: i32 = n as i32 * sc.matching;
        let b = best_one_gap(&q, &t, off, 3, &sc, Fixed::Either, &mut buf).unwrap();
        assert!(b < full);
        // eine Base 8 vor dem Ende fehlt im Read (Deletion): mit Lücke besser als ohne
        let mut q: Vec<u8> = t[off..off + n - 8].to_vec();
        q.extend_from_slice(&t[off + n - 7..off + n + 1]);
        let ungapped: i32 = q
            .iter()
            .zip(&t[off..])
            .map(|(a, b)| if a == b { sc.matching } else { -sc.mismatch })
            .sum();
        let b = best_one_gap(&q, &t, off, 3, &sc, Fixed::Either, &mut buf).unwrap();
        assert_eq!(b, n as i32 * sc.matching - sc.o_del - sc.e_del);
        assert!(b > ungapped);
        // eine eingefügte Base 6 vor dem Ende, als Ende hinter einer Ausrichtung geprüft
        let tail_at = off + 100;
        let mut tail: Vec<u8> = t[tail_at..tail_at + 44].to_vec();
        tail.insert(44, 1 ^ t[tail_at + 44]);
        tail.extend_from_slice(&t[tail_at + 44..tail_at + 49]);
        let b = best_one_gap(&tail, &t, tail_at, 3, &sc, Fixed::Start, &mut buf).unwrap();
        assert_eq!(b, 49 * sc.matching - sc.o_ins - sc.e_ins);
        // zu nah am Rand des Fensters: keine Aussage
        assert!(best_one_gap(&q, &t, 1, 3, &sc, Fixed::Either, &mut buf).is_none());
    }

    #[test]
    fn longer_gaps_in_repeats_are_found_and_the_reach_bound_changes_no_answer() {
        let sc = crate::extend::Scoring::default();
        let mut buf = Vec::new();
        // 4 Basen fehlen 21 vor dem Ende, in einer Wiederholung (ACGT)n, die kurz vor dem Ende
        // aufhört: ohne Lücke und mit 1 bis 3 Basen Lücke Fehler am Ende, mit 4 Basen keine
        let mut state = 0x77_u64;
        let mut t: Vec<u8> = (0..120).map(|_| (xorshift(&mut state) & 3) as u8).collect();
        t.extend((0..48).map(|i| (i % 4) as u8));
        t.extend((0..40).map(|_| (xorshift(&mut state) & 3) as u8));
        let off = 20;
        let n = 151;
        let cut = off + 130;
        let mut q: Vec<u8> = t[off..cut].to_vec();
        q.extend_from_slice(&t[cut + 4..off + n + 4]);
        let full: i32 = q
            .iter()
            .zip(&t[off..])
            .map(|(a, b)| if a == b { sc.matching } else { -sc.mismatch })
            .sum();
        let b3 = best_one_gap(&q, &t, off, 3, &sc, Fixed::Either, &mut buf).unwrap();
        let g = gap_reach(n, full, &sc, 15);
        assert!(g >= 4);
        let b = best_one_gap(&q, &t, off, g, &sc, Fixed::Either, &mut buf).unwrap();
        assert!(b > full);
        assert_eq!(b, n as i32 * sc.matching - sc.o_del - 4 * sc.e_del);
        assert!(b > b3);
        // Die Schranke kürzt nur Lücken, die ohnehin nicht besser als ohne Lücke sein können
        for round in 0..300u64 {
            let mut q: Vec<u8> = t[off..off + n].to_vec();
            for _ in 0..(round % 7) {
                let i = (xorshift(&mut state) % n as u64) as usize;
                q[i] = (xorshift(&mut state) & 3) as u8;
            }
            if round % 3 == 0 {
                let at = (xorshift(&mut state) % 140) as usize + 5;
                let d = (xorshift(&mut state) % 12) as usize + 1;
                q.drain(at..(at + d).min(n));
                q.extend_from_slice(&t[off + n..off + n + (n - q.len())]);
            }
            let full: i32 = q
                .iter()
                .zip(&t[off..])
                .map(|(a, b)| if a == b { sc.matching } else { -sc.mismatch })
                .sum();
            let all = best_one_gap(&q, &t, off, 15, &sc, Fixed::Either, &mut buf).unwrap();
            let g = gap_reach(n, full, &sc, 15);
            let cut = if g == 0 {
                i32::MIN
            } else {
                best_one_gap(&q, &t, off, g, &sc, Fixed::Either, &mut buf).unwrap()
            };
            assert_eq!(all > full, cut > full, "Runde {round}");
        }
    }

    #[test]
    fn reads_with_few_changes_land_at_their_origin_and_copies_go_the_full_way() {
        let mut rng = 0x5eed_4321_u64;
        let mut s: Vec<u8> = (0..60_000)
            .map(|_| b"ACGT"[(xorshift(&mut rng) & 3) as usize])
            .collect();
        // Eine Kopie von 2 kb: Reads daraus sind nicht eindeutig.
        let copy = s[10_000..12_000].to_vec();
        s[40_000..42_000].copy_from_slice(&copy);
        let mut r = gc_index::Reference::new();
        r.push_contig("a", &s);
        let index = FmIndex::build(r, gc_index::Variant::Kompakt, 0).unwrap().0;
        let strobes = StrobeIndex::build(index.reference(), 2);
        let opt = AlignOptions::default();
        let codes: Vec<u8> = s
            .iter()
            .map(|&b| gc_index::reference::encode(b).unwrap())
            .collect();
        let l = index.reference().len();
        let (mut taken, mut tried) = (0, 0);
        for k in 0..400 {
            let at = 20_000 + (xorshift(&mut rng) % 15_000) as usize;
            let mut read = codes[at..at + 150].to_vec();
            // bis zu drei Abweichungen
            for _ in 0..k % 4 {
                let p = (xorshift(&mut rng) % 150) as usize;
                read[p] = (read[p] + 1) % 4;
            }
            let reverse = k % 2 == 1;
            if reverse {
                read = read.iter().rev().map(|&c| 3 - c).collect();
            }
            tried += 1;
            if let Some(reg) = strobe_region(&index, &strobes, &read, &opt) {
                taken += 1;
                let want = if reverse {
                    2 * l - (at as u64 + 150)
                } else {
                    at as u64
                };
                assert_eq!(reg.rb, want, "Read {k}");
                assert_eq!((reg.qb, reg.qe, reg.re - reg.rb), (0, 150, 150));
                assert_eq!(reg.sub, 0);
            }
        }
        assert!(taken * 10 >= tried * 9, "{taken} von {tried}");
        // Reads aus der Kopie: zweiter Ort gleich gut, also voller Weg
        for k in 0..50 {
            let at = 10_200 + k * 30;
            let read = codes[at..at + 150].to_vec();
            assert!(strobe_region(&index, &strobes, &read, &opt).is_none());
        }
    }
}
