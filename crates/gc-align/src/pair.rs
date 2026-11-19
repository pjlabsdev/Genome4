//! Read-Paare (BWA-MEM: `bwamem_pair.c`): Insertgröße schätzen (`mem_pestat`), den Partner in
//! der Nähe suchen (`mem_matesw`), Paare bewerten (`mem_pair`) und beide Reads als SAM schreiben
//! (`mem_sam_pe`).
//!
//! Die beiden Reads eines Paares stammen von den Enden desselben DNA-Stücks. Liegen ihre Orte
//! im erwarteten Abstand und in der erwarteten Richtung, stützen sie einander: Ein Read, der
//! allein an mehrere Stellen passt, bekommt durch seinen Partner einen eindeutigen Ort.

use crate::finalize::{hash_64, mapq, mark_primary, sort_dedup};
use crate::local::local_align;
use crate::region::Region;
use crate::sam::{AlignOptions, Line, Rec, gen_xa, make_line, regions_to_recs};
use gc_index::FmIndex;

/// Einstellungen für Paare. Die Vorgaben sind die von BWA-MEM.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairOptions {
    /// Größerer Abstand zählt nicht in die Schätzung der Insertgröße.
    pub max_ins: i64,
    /// Abzug, wenn die beiden Reads nicht als Paar liegen (`-U`).
    pub pen_unpaired: i32,
    /// Für so viele Bereiche eines Reads wird der Partner gesucht (`-m`).
    pub max_matesw: usize,
}

impl Default for PairOptions {
    fn default() -> Self {
        Self {
            max_ins: 10_000,
            pen_unpaired: 17,
            max_matesw: 50,
        }
    }
}

/// Verteilung des Abstands für eine der vier Richtungen.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct InsertStats {
    /// Zahl der Paare, aus denen geschätzt wurde.
    pub n: usize,
    /// Die Richtung wird nicht benutzt (zu wenige Paare).
    pub failed: bool,
    /// Grenzen, in denen ein Paar als "richtig" gilt.
    pub low: i32,
    pub high: i32,
    pub avg: f64,
    pub std: f64,
    /// 25., 50. und 75. Perzentil.
    pub percentiles: [i32; 3],
}

/// Namen der vier Richtungen: Strang von Read 1 und Read 2 relativ zueinander und wer vorn liegt.
pub const DIRECTIONS: [&str; 4] = ["FF", "FR", "RF", "RR"];

const MIN_RATIO: f64 = 0.8;
const MIN_DIR_CNT: usize = 10;
const MIN_DIR_RATIO: f64 = 0.05;
const OUTLIER_BOUND: f64 = 2.0;
const MAPPING_BOUND: f64 = 3.0;
const MAX_STDDEV: f64 = 4.0;

/// Richtung (0..4) und Abstand zweier Orte im Indextext der Länge `2 * l` (bwa
/// `mem_infer_dir`). Der Abstand zählt zwischen den Anfängen, den zweiten Ort auf den Strang
/// des ersten gespiegelt.
fn infer_dir(l: i64, b1: i64, b2: i64) -> (usize, i64) {
    let (r1, r2) = (b1 >= l, b2 >= l);
    let p2 = if r1 == r2 { b2 } else { (l << 1) - 1 - b2 };
    let dist = (p2 - b1).abs();
    let dir = usize::from(r1 != r2) ^ if p2 > b1 { 0 } else { 3 };
    (dir, dist)
}

/// Score des besten Bereichs, der den ersten im Read stark überlappt (bwa `cal_sub`).
fn cal_sub(r: &[Region], opt: &AlignOptions) -> i32 {
    for x in &r[1..] {
        let b_max = x.qb.max(r[0].qb);
        let e_min = x.qe.min(r[0].qe);
        if e_min > b_max {
            let min_l = (x.qe - x.qb).min(r[0].qe - r[0].qb);
            if (e_min - b_max) as f32 >= min_l as f32 * opt.finalize.mask_level {
                return x.score;
            }
        }
    }
    opt.finalize.min_seed_len * opt.region.scoring.matching
}

/// Schätzt die Insertgröße aus einem Block von Paaren (bwa `mem_pestat`). `regs` enthält die
/// Bereiche aller Reads des Blocks in der Folge Read 1, Read 2, Read 1, ... aus
/// [`crate::sam::find_regions`]. Es zählen nur Paare, bei denen beide Reads einen deutlich
/// besten Ort auf demselben Contig haben.
///
/// Die Schätzung gilt nur für diesen Block; das Ergebnis eines Paares hängt deshalb davon ab,
/// mit welchen anderen Paaren es im Block liegt.
pub fn estimate_insert(
    index: &FmIndex,
    regs: &[Vec<Region>],
    opt: &AlignOptions,
) -> [InsertStats; 4] {
    let l = index.reference().len() as i64;
    let mut isize: [Vec<u64>; 4] = Default::default();
    for pair in regs.as_chunks::<2>().0 {
        let (r0, r1) = (&pair[0], &pair[1]);
        if r0.is_empty() || r1.is_empty() {
            continue;
        }
        if f64::from(cal_sub(r0, opt)) > MIN_RATIO * f64::from(r0[0].score)
            || f64::from(cal_sub(r1, opt)) > MIN_RATIO * f64::from(r1[0].score)
        {
            continue;
        }
        if r0[0].contig != r1[0].contig {
            continue;
        }
        let (dir, is) = infer_dir(l, r0[0].rb as i64, r1[0].rb as i64);
        if is != 0 && is <= opt.pair.max_ins {
            isize[dir].push(is as u64);
        }
    }
    let mut pes = [InsertStats::default(); 4];
    for d in 0..4 {
        let q = &mut isize[d];
        let r = &mut pes[d];
        r.n = q.len();
        if q.len() < MIN_DIR_CNT {
            r.failed = true;
            continue;
        }
        q.sort_unstable();
        let n = q.len() as f64;
        let at = |f: f64| q[(f * n + 0.499) as usize] as i32;
        let (p25, p50, p75) = (at(0.25), at(0.50), at(0.75));
        r.percentiles = [p25, p50, p75];
        let iqr = f64::from(p75 - p25);
        // Mittel und Streuung ohne Ausreißer.
        let low = ((f64::from(p25) - OUTLIER_BOUND * iqr + 0.499) as i32).max(1);
        let high = (f64::from(p75) + OUTLIER_BOUND * iqr + 0.499) as i32;
        let inside = |v: &&u64| **v >= low as u64 && high >= 0 && **v <= high as u64;
        let x = q.iter().filter(inside).count() as f64;
        r.avg = q.iter().filter(inside).map(|&v| v as f64).sum::<f64>() / x;
        r.std = (q
            .iter()
            .filter(inside)
            .map(|&v| (v as f64 - r.avg) * (v as f64 - r.avg))
            .sum::<f64>()
            / x)
            .sqrt();
        // Grenzen für "richtiges Paar": weiter als für die Schätzung.
        r.low = (f64::from(p25) - MAPPING_BOUND * iqr + 0.499) as i32;
        r.high = (f64::from(p75) + MAPPING_BOUND * iqr + 0.499) as i32;
        if f64::from(r.low) > r.avg - MAX_STDDEV * r.std {
            r.low = (r.avg - MAX_STDDEV * r.std + 0.499) as i32;
        }
        if f64::from(r.high) < r.avg + MAX_STDDEV * r.std {
            r.high = (r.avg + MAX_STDDEV * r.std + 0.499) as i32;
        }
        r.low = r.low.max(1);
    }
    let max = pes.iter().map(|p| p.n).max().unwrap_or(0);
    for p in &mut pes {
        if !p.failed && (p.n as f64) < max as f64 * MIN_DIR_RATIO {
            p.failed = true;
        }
    }
    pes
}

/// Sucht den Partner eines Reads dort, wo er nach der Insertgröße liegen müsste (bwa
/// `mem_matesw`). `a` ist ein Bereich des einen Reads, `mate` der andere Read, `ma` dessen
/// bisherige Bereiche. Für jede Richtung, in der `ma` noch keinen passenden Bereich hat, wird
/// der Partner im erwarteten Fenster lokal ausgerichtet; ein Fund kommt zu `ma` dazu.
/// Liefert die Zahl der durchsuchten Fenster.
pub fn mate_rescue(
    index: &FmIndex,
    pes: &[InsertStats; 4],
    a: &Region,
    mate: &[u8],
    ma: &mut Vec<Region>,
    opt: &AlignOptions,
) -> usize {
    let sc = &opt.region.scoring;
    let l = index.reference().len() as i64;
    let l_ms = mate.len() as i64;
    let a_rb = a.rb as i64;
    let mut skip = [false; 4];
    for r in 0..4 {
        skip[r] = pes[r].failed;
    }
    for m in ma.iter() {
        let (r, dist) = infer_dir(l, a_rb, m.rb as i64);
        if dist >= i64::from(pes[r].low) && dist <= i64::from(pes[r].high) {
            skip[r] = true;
        }
    }
    if skip.iter().all(|&s| s) {
        return 0;
    }
    let min_seed = opt.finalize.min_seed_len;
    let mut n = 0;
    for r in 0..4 {
        if skip[r] {
            continue;
        }
        let is_rev = (r >> 1) != (r & 1); // Partner auf dem anderen Strang
        let is_larger = (r >> 1) == 0; // Partner liegt weiter hinten
        let (low, high) = (i64::from(pes[r].low), i64::from(pes[r].high));
        let (mut rb, mut re) = if is_larger {
            (a_rb + low, a_rb + high)
        } else {
            (a_rb - high, a_rb - low)
        };
        // Platz für die Länge des Partners.
        if is_rev {
            rb -= l_ms;
        } else {
            re += l_ms;
        }
        rb = rb.max(0);
        re = re.min(l << 1);
        let mut found_contig = None;
        if rb < re {
            // Auf den Contig begrenzen, in dem die Mitte des Fensters liegt.
            if let Some((contig, beg, end)) = index.contig_text_range(((rb + re) >> 1) as u64) {
                rb = rb.max(beg as i64);
                re = re.min(end as i64);
                found_contig = Some(contig);
            }
        }
        if found_contig == Some(a.contig) && re - rb >= i64::from(min_seed) {
            let seq: Vec<u8> = if is_rev {
                mate.iter()
                    .rev()
                    .map(|&b| if b < 4 { 3 - b } else { 4 })
                    .collect()
            } else {
                mate.to_vec()
            };
            let target: Vec<u8> = (rb..re).map(|p| index.text_base(p as u64)).collect();
            let byte = l_ms * i64::from(sc.matching) < 250;
            let aln = local_align(&seq, &target, sc, byte, min_seed * sc.matching, true);
            if aln.score >= min_seed && aln.qb >= 0 {
                let (qb, qe, tb, te) = (
                    i64::from(aln.qb),
                    i64::from(aln.qe),
                    i64::from(aln.tb),
                    i64::from(aln.te),
                );
                let mut b = Region {
                    rb: 0,
                    re: 0,
                    qb: if is_rev { l_ms - (qe + 1) } else { qb } as i32,
                    qe: if is_rev { l_ms - qb } else { qe + 1 } as i32,
                    contig: a.contig,
                    score: aln.score,
                    truesc: 0,
                    w: 0,
                    seedcov: 0,
                    seedlen0: 0,
                    sub: 0,
                    sub_n: 0,
                    secondary: -1,
                    csub: aln.score2,
                    rescued: true,
                };
                if is_rev {
                    b.rb = ((l << 1) - (rb + te + 1)) as u64;
                    b.re = ((l << 1) - (rb + tb)) as u64;
                } else {
                    b.rb = (rb + tb) as u64;
                    b.re = (rb + te + 1) as u64;
                }
                b.seedcov = ((b.re - b.rb) as i32).min(b.qe - b.qb) >> 1;
                // Nach Score einordnen: vor den ersten Bereich mit kleinerem Score.
                let at = ma
                    .iter()
                    .position(|m| m.score < b.score)
                    .unwrap_or(ma.len());
                ma.insert(at, b);
            }
            n += 1;
        }
        if n > 0 {
            sort_dedup(ma, index, None, sc, &opt.finalize);
        }
    }
    n
}

/// Das beste Paar aus den Bereichen beider Reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BestPair {
    /// Score des Paares: Summe beider Scores, vermindert um einen Abzug, der mit der
    /// Abweichung des Abstands vom Mittel wächst.
    pub score: i32,
    /// Score des zweitbesten Paares (0 = keines) und Zahl der fast so guten.
    pub sub: i32,
    pub n_sub: i32,
    /// Index des gewählten Bereichs je Read.
    pub z: [usize; 2],
}

/// Sucht unter allen Kombinationen von Bereichen der beiden Reads das beste Paar in
/// erlaubter Richtung und erlaubtem Abstand (bwa `mem_pair`). `id` ist die laufende Nummer des
/// Paares; sie entscheidet zwischen gleich guten Paaren.
pub fn best_pair(
    index: &FmIndex,
    pes: &[InsertStats; 4],
    a: [&[Region]; 2],
    id: u64,
    opt: &AlignOptions,
) -> Option<BestPair> {
    let sc = &opt.region.scoring;
    let l = index.reference().len();
    // Alle Bereiche, nach Ort auf dem Vorwärtsstrang geordnet. y: Score, Index, Strang, Read.
    let mut v: Vec<(u64, u64)> = Vec::new();
    for (r, regs) in a.iter().enumerate() {
        for (i, e) in regs.iter().enumerate() {
            let fwd = if e.rb < l { e.rb } else { (l << 1) - 1 - e.rb };
            let offset = index.reference().contigs[e.contig as usize].offset;
            let x = u64::from(e.contig) << 32 | (fwd - offset);
            let y = (e.score as u64) << 32 | (i as u64) << 2 | u64::from(e.rb >= l) << 1 | r as u64;
            v.push((x, y));
        }
    }
    v.sort_unstable();
    // bwa führt die Paarnummer hier als 32-Bit-Zahl mit Vorzeichen; für gleiche Wahl bei sehr
    // vielen Paaren wird das nachgebildet.
    let id_bits = i64::from((id as i32).wrapping_shl(8)) as u64;
    let mut last = [-1i64; 4];
    let mut u: Vec<(u64, u64)> = Vec::new();
    for i in 0..v.len() {
        for r in 0..2usize {
            let dir = r << 1 | (v[i].1 >> 1 & 1) as usize;
            if pes[dir].failed {
                continue;
            }
            let which = (r << 1) as u64 | ((v[i].1 & 1) ^ 1);
            let mut k = last[which as usize];
            while k >= 0 {
                let vk = v[k as usize];
                k -= 1;
                if vk.1 & 3 != which {
                    continue;
                }
                let dist = v[i].0 as i64 - vk.0 as i64;
                if dist > i64::from(pes[dir].high) {
                    break;
                }
                if dist < i64::from(pes[dir].low) {
                    continue;
                }
                let ns = (dist as f64 - pes[dir].avg) / pes[dir].std;
                // 0,721 = 1 / ln 4: Die Wahrscheinlichkeit einer so großen Abweichung, in
                // Score-Punkte umgerechnet.
                let q = ((v[i].1 >> 32) as f64
                    + (vk.1 >> 32) as f64
                    + 0.721
                        * (2.0 * libm::erfc(ns.abs() * std::f64::consts::FRAC_1_SQRT_2)).ln()
                        * f64::from(sc.matching)
                    + 0.499) as i32;
                let q = q.max(0) as u64;
                let y = ((k + 1) as u64) << 32 | i as u64;
                u.push((q << 32 | (hash_64(y ^ id_bits) & 0xffff_ffff), y));
            }
        }
        last[(v[i].1 & 3) as usize] = i as i64;
    }
    if u.is_empty() {
        return None;
    }
    u.sort_unstable();
    let top = u[u.len() - 1];
    let (i, k) = ((top.1 >> 32) as usize, (top.1 & 0xffff_ffff) as usize);
    let mut z = [0usize; 2];
    z[(v[i].1 & 1) as usize] = ((v[i].1 & 0xffff_ffff) >> 2) as usize;
    z[(v[k].1 & 1) as usize] = ((v[k].1 & 0xffff_ffff) >> 2) as usize;
    let near = (sc.matching + sc.mismatch)
        .max(sc.o_del + sc.e_del)
        .max(sc.o_ins + sc.e_ins);
    let sub = if u.len() > 1 {
        (u[u.len() - 2].0 >> 32) as i32
    } else {
        0
    };
    let n_sub = u[..u.len() - 1]
        .iter()
        .filter(|p| sub - (p.0 >> 32) as i32 <= near)
        .count() as i32;
    Some(BestPair {
        score: (top.0 >> 32) as i32,
        sub,
        n_sub,
        z,
    })
}

fn raw_mapq(diff: i32, a: i32) -> i32 {
    (6.02 * f64::from(diff) / f64::from(a) + 0.499) as i32
}

/// Ein Read eines Paares mit seinen Bereichen aus [`crate::sam::find_regions`].
pub struct PairedRead<'a> {
    pub name: &'a str,
    pub read: &'a [u8],
    /// Leer: keine Qualitätswerte.
    pub qual: &'a [u8],
    pub regions: Vec<Region>,
    pub frac_rep: f32,
}

/// Schreibt die SAM-Zeilen eines Paares nach `out`; siehe [`pair_lines`].
pub fn write_sam_pair(
    out: &mut String,
    index: &FmIndex,
    pes: &[InsertStats; 4],
    id: u64,
    s: &mut [PairedRead<'_>; 2],
    opt: &AlignOptions,
) -> usize {
    let (lines, n) = pair_lines(index, pes, id, s, opt);
    for line in &lines {
        line.write_sam(out, index);
    }
    n
}

/// Die Ausgabezeilen eines Paares, erst Read 1, dann Read 2 (bwa `mem_sam_pe`), und die
/// Zahl der Fenster, in denen ein Partner gesucht wurde.
///
/// Ablauf: Partnersuche für die guten Bereiche beider Reads, Rangfolge je Read, dann das
/// beste Paar. Ist das Paar besser als die beiden besten Einzelorte abzüglich `pen_unpaired`,
/// werden die Orte des Paares gemeldet, sonst die besten Einzelorte. Hat ein Read mehrere
/// Teile an verschiedenen Orten, wird er wie ein einzelner Read gemeldet.
pub fn pair_lines(
    index: &FmIndex,
    pes: &[InsertStats; 4],
    id: u64,
    s: &mut [PairedRead<'_>; 2],
    opt: &AlignOptions,
) -> (Vec<Line>, usize) {
    let mut out: Vec<Line> = Vec::with_capacity(2);
    let sc = &opt.region.scoring;
    let l = index.reference().len() as i64;
    let mut n = 0;
    // Partnersuche, ausgehend von den Bereichen, die fast so gut sind wie der beste.
    let good: [Vec<Region>; 2] = [0, 1].map(|i| {
        let regs = &s[i].regions;
        regs.iter()
            .filter(|r| r.score >= regs[0].score - opt.pair.pen_unpaired)
            .copied()
            .collect()
    });
    for i in 0..2 {
        for b in good[i].iter().take(opt.pair.max_matesw) {
            let mate = s[1 - i].read;
            n += mate_rescue(index, pes, b, mate, &mut s[1 - i].regions, opt);
        }
    }
    for (i, p) in s.iter_mut().enumerate() {
        mark_primary(&mut p.regions, id << 1 | i as u64, sc, &opt.finalize);
    }
    let own_mapq =
        |p: &PairedRead<'_>, k: usize| mapq(&p.regions[k], p.frac_rep, sc, &opt.finalize);
    let mut extra_flag = 1u32;

    let pair = if s[0].regions.is_empty() || s[1].regions.is_empty() {
        None
    } else {
        best_pair(index, pes, [&s[0].regions, &s[1].regions], id, opt).filter(|p| p.score > 0)
    };
    // Mehrere Teile an verschiedenen Orten: nicht als Paar melden.
    let is_multi = |p: &PairedRead<'_>| {
        p.regions[1..]
            .iter()
            .any(|r| r.secondary < 0 && r.score >= opt.sam.min_score)
    };
    if let Some(pair) = pair.filter(|_| !is_multi(&s[0]) && !is_multi(&s[1])) {
        let o = pair.score;
        let score_un = s[0].regions[0].score + s[1].regions[0].score - opt.pair.pen_unpaired;
        let subo = pair.sub.max(score_un);
        let mut q_pe = raw_mapq(o - subo, sc.matching);
        if pair.n_sub > 0 {
            q_pe -= (4.343 * f64::from(pair.n_sub + 1).ln() + 0.499) as i32;
        }
        q_pe = q_pe.clamp(0, 60);
        let frac = |p: &PairedRead<'_>| {
            if p.regions[0].rescued {
                0.0
            } else {
                p.frac_rep
            }
        };
        q_pe =
            (f64::from(q_pe) * (1.0 - 0.5 * f64::from(frac(&s[0]) + frac(&s[1]))) + 0.499) as i32;
        let mut z = pair.z;
        let mut q_se = [0i32; 2];
        // Rangfolge vor dem Umstellen: Sie bestimmt, was in XA steht.
        let mut parents: [Vec<i32>; 2] =
            [0, 1].map(|i| s[i].regions.iter().map(|r| r.secondary).collect());
        if o > score_un {
            // Das Paar gewinnt.
            for i in 0..2 {
                let k = s[i].regions[z[i]].secondary;
                if k >= 0 {
                    s[i].regions[z[i]].sub = s[i].regions[k as usize].score;
                    s[i].regions[z[i]].secondary = -2;
                }
                q_se[i] = i32::from(own_mapq(&s[i], z[i]));
            }
            for i in 0..2 {
                // Der Ort wird durch den Partner sicherer, aber um höchstens 40.
                q_se[i] = if q_se[i] > q_pe {
                    q_se[i]
                } else {
                    q_pe.min(q_se[i] + 40)
                };
                // Obergrenze aus dem zweitbesten Score der Partnersuche (dicht benachbarte
                // Wiederholungen).
                let c = &s[i].regions[z[i]];
                q_se[i] = q_se[i].min(raw_mapq(c.score - c.csub, sc.matching));
            }
            extra_flag |= 2;
        } else {
            // Die besten Einzelorte gewinnen.
            z = [0, 0];
            for i in 0..2 {
                q_se[i] = i32::from(own_mapq(&s[i], 0));
            }
        }
        let mut h: Vec<Rec> = Vec::with_capacity(2);
        for i in 0..2 {
            // Der gewählte Bereich tritt an die Stelle seines bisher übergeordneten.
            let p = &mut parents[i];
            let k = p[z[i]];
            if k >= 0 {
                for (j, pj) in p.iter_mut().enumerate() {
                    if *pj == k || j == k as usize {
                        *pj = z[i] as i32;
                    }
                }
                p[z[i]] = -1;
            }
            let mut xa = gen_xa(index, s[i].read, &s[i].regions, p, opt);
            let mut rec = match crate::align::region_to_alignment(
                index,
                s[i].read,
                &s[i].regions[z[i]],
                sc,
            ) {
                Some(aln) => Rec::mapped(aln, q_se[i].max(0) as u8),
                None => Rec::unmapped(),
            };
            rec.flag |= 0x40 << i | extra_flag;
            rec.xa = std::mem::take(&mut xa[z[i]]);
            h.push(rec);
        }
        for i in 0..2 {
            out.push(make_line(
                index,
                s[i].name,
                s[i].read,
                s[i].qual,
                std::slice::from_ref(&h[i]),
                0,
                Some(&h[1 - i]),
            ));
        }
        return (out, n);
    }

    // Kein Paar: jeder Read wie ein einzelner, mit den Feldern des Partners.
    let h: Vec<Rec> = (0..2)
        .map(|i| {
            let p = &s[i];
            match p.regions.first() {
                Some(r) if r.score >= opt.sam.min_score => {
                    match crate::align::region_to_alignment(index, p.read, r, sc) {
                        Some(aln) => Rec::mapped(aln, own_mapq(p, 0)),
                        None => Rec::unmapped(),
                    }
                }
                _ => Rec::unmapped(),
            }
        })
        .collect();
    if h[0].contig.is_some() && h[0].contig == h[1].contig {
        // Liegen die besten Einzelorte wie ein richtiges Paar, wird das vermerkt.
        let (d, dist) = infer_dir(l, s[0].regions[0].rb as i64, s[1].regions[0].rb as i64);
        if !pes[d].failed && dist >= i64::from(pes[d].low) && dist <= i64::from(pes[d].high) {
            extra_flag |= 2;
        }
    }
    for i in 0..2 {
        let p = &s[i];
        let flag = (0x40 << i) | 1 | extra_flag;
        let recs = regions_to_recs(index, p.read, &p.regions, p.frac_rep, flag, opt);
        for which in 0..recs.len() {
            out.push(make_line(
                index,
                p.name,
                p.read,
                p.qual,
                &recs,
                which,
                Some(&h[1 - i]),
            ));
        }
    }
    (out, n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::find_regions;
    use gc_index::{Reference, Variant};

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn setup() -> (FmIndex, Vec<u8>) {
        let mut rng = 0x00C0_FFEE_1234_5678_u64;
        let s: Vec<u8> = (0..20_000)
            .map(|_| b"ACGT"[(xorshift(&mut rng) & 3) as usize])
            .collect();
        let mut r = Reference::new();
        r.push_contig("a", &s);
        let index = FmIndex::build(r, Variant::Kompakt, 0).unwrap().0;
        let codes = s
            .iter()
            .map(|&b| gc_index::reference::encode(b).unwrap())
            .collect();
        (index, codes)
    }

    fn revcomp(read: &[u8]) -> Vec<u8> {
        read.iter().rev().map(|&b| 3 - b).collect()
    }

    /// Paar aus einem Stück ab `pos` mit Länge `size`: Read 1 vorn, Read 2 hinten auf dem
    /// Gegenstrang, je 100 Basen.
    fn pair_at(codes: &[u8], pos: usize, size: usize) -> [Vec<u8>; 2] {
        [
            codes[pos..pos + 100].to_vec(),
            revcomp(&codes[pos + size - 100..pos + size]),
        ]
    }

    /// Richtet die Paare als einen Block aus und liefert je Paar die Felder der SAM-Zeilen.
    fn run(index: &FmIndex, pairs: &[[Vec<u8>; 2]]) -> (Vec<Vec<Vec<String>>>, [InsertStats; 4]) {
        let opt = AlignOptions::default();
        let mut regs = Vec::new();
        let mut fracs = Vec::new();
        for p in pairs {
            for read in p {
                let (r, f) = find_regions(index, read, &opt);
                regs.push(r);
                fracs.push(f);
            }
        }
        let pes = estimate_insert(index, &regs, &opt);
        let mut regs = regs.into_iter();
        let mut all = Vec::new();
        for (id, p) in pairs.iter().enumerate() {
            let mut s = [0, 1].map(|k| PairedRead {
                name: "p",
                read: &p[k],
                qual: &[],
                regions: regs.next().unwrap(),
                frac_rep: fracs[2 * id + k],
            });
            let mut out = String::new();
            write_sam_pair(&mut out, index, &pes, id as u64, &mut s, &opt);
            all.push(
                out.lines()
                    .map(|l| l.split('\t').map(str::to_string).collect())
                    .collect(),
            );
        }
        (all, pes)
    }

    /// 40 gewöhnliche Paare mit Stücklängen um 400.
    fn ordinary(codes: &[u8]) -> Vec<[Vec<u8>; 2]> {
        (0..40)
            .map(|i| pair_at(codes, 300 + i * 450, 380 + (i * 7) % 41))
            .collect()
    }

    fn tag(line: &[String], name: &str) -> Option<String> {
        line.iter()
            .find(|f| f.starts_with(name))
            .map(|f| f[5..].to_string())
    }

    #[test]
    fn direction_and_distance_of_two_places() {
        let l = 1000;
        // Read 1 vorwärts bei 100, Read 2 auf dem Gegenstrang, sein Anfang (rechtes Ende auf
        // dem Vorwärtsstrang) bei 499: Richtung FR, Abstand 399.
        assert_eq!(infer_dir(l, 100, 2 * l - 1 - 499), (1, 399));
        // Umgekehrt gesehen: Read 1 auf dem Gegenstrang, der Partner liegt davor.
        assert_eq!(infer_dir(l, 2 * l - 1 - 499, 100), (1, 399));
        // Beide vorwärts.
        assert_eq!(infer_dir(l, 100, 300), (0, 200));
        assert_eq!(infer_dir(l, 300, 100), (3, 200));
    }

    #[test]
    fn insert_size_is_estimated_from_the_block() {
        let (index, codes) = setup();
        let (_, pes) = run(&index, &ordinary(&codes));
        assert!(pes[0].failed && pes[2].failed && pes[3].failed);
        let fr = pes[1];
        assert!(!fr.failed);
        assert_eq!(fr.n, 40);
        // Abstand der Anfänge = Stücklänge - 1, Stücklängen 380 bis 420.
        assert!((fr.avg - 399.0).abs() < 3.0, "{fr:?}");
        assert!(fr.low < 379 && fr.high > 419 && fr.high < 600, "{fr:?}");
        // Unter zehn Paaren gibt es keine Schätzung.
        let (lines, pes) = run(&index, &ordinary(&codes)[..9]);
        assert!(pes.iter().all(|p| p.failed));
        assert_eq!(
            lines[0][0][1], "97",
            "ohne Schätzung kein richtiges Paar (Bit 2 fehlt)"
        );
    }

    #[test]
    fn an_ordinary_pair_carries_the_fields_of_its_mate() {
        let (index, codes) = setup();
        let (lines, _) = run(&index, &ordinary(&codes));
        // Paar 0: Stück ab 300 mit Länge 380.
        let (a, b) = (&lines[0][0], &lines[0][1]);
        assert_eq!(lines[0].len(), 2);
        assert_eq!(a[1..9], ["99", "a", "301", "60", "100M", "=", "581", "380"]);
        assert_eq!(
            b[1..9],
            ["147", "a", "581", "60", "100M", "=", "301", "-380"]
        );
        assert_eq!(tag(a, "MC:Z:").unwrap(), "100M");
        assert_eq!(tag(b, "MQ:i:").unwrap(), "60");
    }

    #[test]
    fn a_mate_without_seeds_is_found_near_its_partner() {
        let (index, codes) = setup();
        let mut pairs = ordinary(&codes);
        // Read 2 von Paar 5: jede 12. Base vertauscht. Kein Stück von 19 Basen stimmt mehr,
        // allein findet der Read keinen Ort.
        for i in (5..100).step_by(12) {
            pairs[5][1][i] = (pairs[5][1][i] + 1) & 3;
        }
        let opt = AlignOptions::default();
        assert!(find_regions(&index, &pairs[5][1], &opt).0.is_empty());
        let (lines, _) = run(&index, &pairs);
        let (a, b) = (&lines[5][0], &lines[5][1]);
        let size = 380 + 35;
        let pos1 = (300 + 5 * 450 + 1).to_string();
        let pos2 = (300 + 5 * 450 + size - 100 + 1).to_string();
        assert_eq!(a[1..4], ["99", "a", &pos1[..]]);
        assert_eq!(b[1..4], ["147", "a", &pos2[..]]);
        assert_eq!(b[5], "100M");
        assert_eq!(tag(b, "NM:i:").unwrap(), "8");
        assert_eq!(tag(b, "AS:i:").unwrap(), "60");
    }

    #[test]
    fn a_mate_without_a_place_is_listed_at_the_place_of_its_partner() {
        let (index, codes) = setup();
        let mut pairs = ordinary(&codes);
        let mut rng = 0x5151_AAAA_0F0F_1234_u64;
        pairs[7][1] = (0..100).map(|_| (xorshift(&mut rng) & 3) as u8).collect();
        let (lines, _) = run(&index, &pairs);
        let (a, b) = (&lines[7][0], &lines[7][1]);
        let pos = (300 + 7 * 450 + 1).to_string();
        // Read 1: Partner ohne Ort (8), erster im Paar (64).
        assert_eq!(
            a[1..9],
            ["73", "a", &pos[..], "60", "100M", "=", &pos[..], "0"]
        );
        // Read 2: ohne Ort (4), zweiter im Paar (128), steht beim Ort von Read 1.
        assert_eq!(
            b[1..9],
            ["133", "a", &pos[..], "0", "*", "=", &pos[..], "0"]
        );
        assert_eq!(tag(b, "MC:Z:").unwrap(), "100M");
        assert_eq!(tag(a, "MC:Z:"), None);
    }

    #[test]
    fn the_mate_decides_between_two_equally_good_places() {
        // Read 1 passt an zwei Stellen gleich gut, Read 2 nur neben eine davon.
        let mut rng = 0x1357_9BDF_2468_ACE0_u64;
        let mut s: Vec<u8> = (0..20_000)
            .map(|_| b"ACGT"[(xorshift(&mut rng) & 3) as usize])
            .collect();
        let copy = s[10_000..10_100].to_vec();
        s[15_000..15_100].copy_from_slice(&copy);
        let mut r = Reference::new();
        r.push_contig("a", &s);
        let index = FmIndex::build(r, Variant::Kompakt, 0).unwrap().0;
        let codes: Vec<u8> = s
            .iter()
            .map(|&b| gc_index::reference::encode(b).unwrap())
            .collect();
        let mut pairs: Vec<[Vec<u8>; 2]> = (0..40)
            .map(|i| pair_at(&codes, 300 + i * 200, 380 + (i * 7) % 41))
            .collect();
        pairs.push(pair_at(&codes, 15_000, 400));
        let (lines, _) = run(&index, &pairs);
        let a = &lines[40][0];
        assert_eq!(a[1..4], ["99", "a", "15001"]);
        // Allein hätte der Read MAPQ 0; der Partner macht den Ort sicher. Der andere Ort
        // steht in XA.
        assert!(a[4].parse::<u8>().unwrap() > 20, "{a:?}");
        assert_eq!(tag(a, "XA:Z:").unwrap(), "a,+10001,100M,0;");
    }
}
