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
}

pub static COUNTS: StrobeCounts = StrobeCounts {
    tried: AtomicU64::new(0),
    taken: AtomicU64::new(0),
    blocked: AtomicU64::new(0),
    no_hit: AtomicU64::new(0),
    clipped: AtomicU64::new(0),
    too_weak: AtomicU64::new(0),
    second_close: AtomicU64::new(0),
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
    ]
}

thread_local! {
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
}

/// Ein Bereich für den Read oder `None` (dann der volle Weg).
pub fn strobe_region(
    index: &FmIndex,
    strobes: &StrobeIndex,
    read: &[u8],
    opt: &AlignOptions,
) -> Option<Region> {
    bump(&COUNTS.tried);
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
            let r = crate::local::local_align(read, &target, sc, byte, 0, true);
            let cand = (
                r.score,
                tb + r.tb as u64,
                tb + r.te as u64 + 1,
                r.qb,
                r.qe + 1,
            );
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
            if d < -opt.region.pen_clip5 {
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
            if d < -opt.region.pen_clip3 {
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
        if score - second < so.min_gap {
            bump(&COUNTS.second_close);
            return None;
        }
        // Contig des Bereichs; über eine Contig-Grenze geht der volle Weg.
        let hit = index.hit_at(rb, re - rb)?;
        bump(&COUNTS.taken);
        let min_sub = opt.finalize.min_seed_len * sc.matching;
        Some(Region {
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
            sub: if second > min_sub { second } else { 0 },
            sub_n: 0,
            secondary: -1,
            csub: 0,
            rescued: false,
        })
    })
}

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
