//! Versuch Schubfach-Paarprüfung (`GC_STROBE_PAAR_TAUBE=1`) für die unsicheren Paare des
//! Strobe-Paarwegs: Paare, bei denen beide Reads in Wiederholungen liegen.
//!
//! Für Ort und MAPQ eines Paares zählen nur der beste und der zweitbeste Paar-Score. Jeder Read
//! wird in getrennte Stücke zu [`K`] Basen zerlegt (bei 150 Basen sieben). Eine Kopie mit
//! höchstens sechs Abweichungen enthält nach dem Schubfachprinzip mindestens eines der Stücke
//! exakt. Die Stücke des Reads mit den wenigeren Fundstellen werden im Index nachgeschlagen
//! (höchstens `GC_STROBE_PAAR_TAUBE_MAX` Fundstellen, sonst der bisherige Weg); jede Diagonale
//! wird ausgerichtet. Der Partner wird im Paarfenster über seine eigenen Stücke gesucht statt mit
//! einer Ausrichtung über das ganze Fenster. Die Orte werden nach ihrem Score abgearbeitet; es
//! wird aufgehört, sobald Score plus der höchstmögliche Partner-Score das bisher zweitbeste Paar
//! nicht mehr erreicht (verzweigen und begrenzen).

use crate::region::Region;
use crate::sam::AlignOptions;
use gc_index::FmIndex;

/// Länge der Stücke.
pub const K: usize = 21;

pub fn mode() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_PAAR_TAUBE").is_ok_and(|v| v == "1"))
}

fn cap() -> u64 {
    static V: std::sync::OnceLock<u64> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("GC_STROBE_PAAR_TAUBE_MAX")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(500)
    })
}

pub static TRIED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static TAKEN: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static CROWDED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

fn pack(s: &[u8]) -> u64 {
    s.iter().fold(0u64, |h, &c| h << 2 | u64::from(c))
}

/// Strobe-Kodierung (Strang in Bit 63, Diagonale + 2^32) des Reads, der im Indextext bei `d`
/// beginnt.
fn anchor_of(d: u64, len: u64, l: u64) -> u64 {
    if d < l {
        d + (1 << 32)
    } else {
        1 << 63 | ((2 * l - d - len) + (1 << 32))
    }
}

fn same_place(o: &Region, r: &Region) -> bool {
    let ov = o.re.min(r.re) as i64 - o.rb.max(r.rb) as i64;
    ov * 100 >= (o.re - o.rb).min(r.re - r.rb) as i64 * 95
}

fn add(v: &mut Vec<Region>, r: Region) {
    match v.iter_mut().find(|o| same_place(o, &r)) {
        Some(o) if o.score >= r.score => {}
        Some(o) => *o = r,
        None => v.push(r),
    }
}

/// Die Bereiche beider Reads; `None`: bisheriger Weg. `known` sind die schon gefundenen Orte
/// (aus den Strobe-Ankern), sie kommen dazu.
pub fn pigeonhole_pair(
    index: &FmIndex,
    reads: [&[u8]; 2],
    known: &[Vec<Region>; 2],
    max_ins: i64,
    opt: &AlignOptions,
) -> Option<[Vec<Region>; 2]> {
    use std::sync::atomic::Ordering::Relaxed;
    TRIED.fetch_add(1, Relaxed);
    let l = index.reference().len();
    let sc = &opt.region.scoring;
    // Stücke und ihre Zeilenbereiche je Read
    let mut ivs: [Vec<(usize, (u64, u64))>; 2] = [Vec::new(), Vec::new()];
    let mut total = [0u64; 2];
    for r in 0..2 {
        let read = reads[r];
        let pieces: Vec<&[u8]> = (0..read.len() / K)
            .map(|i| &read[i * K..(i + 1) * K])
            .collect();
        let mut out = vec![(0, 0); pieces.len()];
        index.intervals_many(&pieces, &mut out);
        for (i, iv) in out.into_iter().enumerate() {
            total[r] += iv.1 - iv.0;
            ivs[r].push((i * K, iv));
        }
    }
    let a = if total[0] <= total[1] { 0 } else { 1 };
    let b = 1 - a;
    if total[a] > cap() || total[a] == 0 {
        CROWDED.fetch_add(1, Relaxed);
        return None;
    }
    let (la, lb) = (reads[a].len() as u64, reads[b].len() as u64);
    // Diagonalen des selteneren Reads
    let mut rows = Vec::new();
    let mut offs = Vec::new();
    for &(off, (lo, hi)) in &ivs[a] {
        for row in lo..hi {
            rows.push(row);
            offs.push(off as u64);
        }
    }
    let mut pos = Vec::with_capacity(rows.len());
    index.locate_many(&rows, &mut pos);
    let mut diags: Vec<u64> = pos
        .iter()
        .zip(&offs)
        .filter(|&(&p, &o)| p >= o)
        .map(|(&p, &o)| p - o)
        .collect();
    diags.sort_unstable();
    diags.dedup_by(|x, y| *x - *y <= 8);
    let mut cands: Vec<Region> = Vec::new();
    for &d in &diags {
        if d + la > 2 * l || (d < l) != (d + la <= l) {
            continue;
        }
        if let Some(r) = crate::paarweg::align_anchor(index, reads[a], anchor_of(d, la, l), opt) {
            add(&mut cands, r);
        }
    }
    for r in &known[a] {
        add(&mut cands, *r);
    }
    if cands.is_empty() {
        return None;
    }
    cands.sort_unstable_by(|x, y| y.score.cmp(&x.score).then(x.rb.cmp(&y.rb)));
    // Stücke des Partners
    let bpieces: Vec<(u64, u64)> = (0..reads[b].len() / K)
        .map(|i| (pack(&reads[b][i * K..(i + 1) * K]), (i * K) as u64))
        .collect();
    let mask = (1u64 << (2 * K)) - 1;
    let b_max = lb as i32 * sc.matching;
    let (mut p1, mut p2) = (i32::MIN, i32::MIN);
    let mut out_a: Vec<Region> = Vec::new();
    let mut out_b: Vec<Region> = known[b].clone();
    for c in &cands {
        if c.score.saturating_add(b_max) < p2 {
            break;
        }
        out_a.push(*c);
        // Fenster des Partners im Indextext: auf dem anderen Strang, dahinter (vorwärts gesehen)
        let (wb, we) = if c.rb < l {
            let f = c.rb as i64;
            let (fb, fe) = (f, f + max_ins);
            (
                (2 * l as i64 - fe).max(l as i64),
                (2 * l as i64 - fb).min(2 * l as i64),
            )
        } else {
            let fend = 2 * l as i64 - c.rb as i64;
            ((fend - max_ins).max(0), fend.min(l as i64))
        };
        if we - wb < lb as i64 {
            continue;
        }
        let text = index.text_range(wb as u64, we as u64);
        let mut h = 0u64;
        let mut run = 0usize;
        let mut found: Vec<u64> = Vec::new();
        for (t, &x) in text.iter().enumerate() {
            if x > 3 {
                run = 0;
                continue;
            }
            h = (h << 2 | u64::from(x)) & mask;
            run += 1;
            if run >= K {
                let start = (wb as u64) + (t + 1 - K) as u64;
                for &(key, off) in &bpieces {
                    if key == h && start >= off {
                        found.push(start - off);
                    }
                }
            }
        }
        found.sort_unstable();
        found.dedup_by(|x, y| *x - *y <= 8);
        for &d in &found {
            if d + lb > 2 * l || (d < l) != (d + lb <= l) {
                continue;
            }
            if let Some(r) = crate::paarweg::align_anchor(index, reads[b], anchor_of(d, lb, l), opt)
            {
                let s = c.score + r.score;
                if s > p1 {
                    p2 = p1;
                    p1 = s;
                } else if s > p2 {
                    p2 = s;
                }
                add(&mut out_b, r);
            }
        }
    }
    if out_a.is_empty() || out_b.is_empty() {
        return None;
    }
    for v in [&mut out_a, &mut out_b] {
        v.sort_unstable_by(|x, y| y.score.cmp(&x.score).then(x.rb.cmp(&y.rb)));
    }
    TAKEN.fetch_add(1, Relaxed);
    Some(if a == 0 {
        [out_a, out_b]
    } else {
        [out_b, out_a]
    })
}
