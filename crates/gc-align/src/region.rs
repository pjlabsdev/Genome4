//! Aus einer Kette Alignment-Bereiche machen (BWA-MEM: `mem_chain2aln`).
//!
//! Für jeden Seed der Kette, den längsten zuerst, wird nach links und rechts verlängert. Das
//! Ergebnis ist ein Bereich im Read und in der Referenz mit Score, noch ohne CIGAR. Seeds, die
//! schon in einem früheren Bereich derselben Kette liegen und zu ihm passen, werden übersprungen.

use crate::chain::Chain;
use crate::extend::{Scoring, extend};
use crate::seed::Seed;
use gc_index::FmIndex;

/// Einstellungen. Die Vorgaben sind die von BWA-MEM.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionOptions {
    pub scoring: Scoring,
    /// Abzug für abgeschnittene Read-Enden links und rechts (`-L`). Der Read wird bis zum Ende
    /// ausgerichtet, wenn das höchstens so viele Punkte kostet.
    pub pen_clip5: i32,
    pub pen_clip3: i32,
}

impl Default for RegionOptions {
    fn default() -> Self {
        Self {
            scoring: Scoring::default(),
            pen_clip5: 5,
            pen_clip3: 5,
        }
    }
}

/// Ein Alignment-Bereich: `read[qb..qe]` gegen den Indextext `[rb, re)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Region {
    pub rb: u64,
    pub re: u64,
    pub qb: i32,
    pub qe: i32,
    pub contig: u32,
    /// Bester lokaler Score.
    pub score: i32,
    /// Score des gemeldeten Bereichs; kleiner als `score`, wenn bis zum Read-Ende verlängert wurde.
    pub truesc: i32,
    /// Verwendete Bandbreite.
    pub w: i32,
    /// Basen der Ketten-Seeds, die im Bereich liegen.
    pub seedcov: i32,
    /// Länge des Seeds, von dem die Verlängerung ausging.
    pub seedlen0: i32,
    /// Score des besten anderen Bereichs, der diesen im Read stark überlappt (0 = keiner).
    pub sub: i32,
    /// Zahl der überlappenden Bereiche, die fast so gut sind wie dieser.
    pub sub_n: i32,
    /// Index des Bereichs, dem dieser untergeordnet ist; -1 bei einem primären Bereich.
    pub secondary: i32,
}

/// Längste Lücke, die sich bei `qlen` Read-Basen noch lohnen kann, begrenzt auf das doppelte Band.
fn cal_max_gap(sc: &Scoring, qlen: i32) -> i32 {
    let l_del = (f64::from(qlen * sc.matching - sc.o_del) / f64::from(sc.e_del) + 1.0) as i32;
    let l_ins = (f64::from(qlen * sc.matching - sc.o_ins) / f64::from(sc.e_ins) + 1.0) as i32;
    l_del.max(l_ins).max(1).min(sc.band_width << 1)
}

/// Verlängert mit der Vorgabe-Bandbreite und, wenn das Maximum nahe am Rand des Bandes lag, noch
/// einmal mit der doppelten. `score` ist der bisherige Score des Bereichs und wird ersetzt.
fn extend_retry(
    query: &[u8],
    target: &[u8],
    sc: &Scoring,
    end_bonus: i32,
    h0: i32,
    score: &mut i32,
) -> (crate::extend::Extension, i32) {
    let mut i = 0;
    loop {
        let prev = *score;
        let w = sc.band_width << i;
        let band = Scoring {
            band_width: w,
            ..sc.clone()
        };
        let x = extend(query, target, &band, end_bonus, h0);
        *score = x.score;
        if i == 1 || *score == prev || x.max_off < (w >> 1) + (w >> 2) {
            return (x, w);
        }
        i += 1;
    }
}

/// Verlängert die Seeds einer Kette zu Bereichen und hängt sie an `out` an. `read` sind Codes
/// 0..3, alles andere N.
pub fn chain_to_regions(
    index: &FmIndex,
    read: &[u8],
    chain: &Chain,
    opt: &RegionOptions,
    out: &mut Vec<Region>,
) {
    let sc = &opt.scoring;
    let l_query = read.len() as i32;
    let l = index.reference().len();
    if chain.seeds.is_empty() {
        return;
    }
    // Größtes Stück Referenz, das eine Verlängerung erreichen kann.
    let (mut r0, mut r1) = (2 * l as i64, 0i64);
    for t in &chain.seeds {
        let (qbeg, len) = (t.qbeg as i32, t.len as i32);
        let b = t.rbeg as i64 - i64::from(qbeg + cal_max_gap(sc, qbeg));
        let rest = l_query - qbeg - len;
        let e = t.rbeg as i64 + i64::from(len) + i64::from(rest + cal_max_gap(sc, rest));
        r0 = r0.min(b);
        r1 = r1.max(e);
    }
    let first = chain.seeds[0].rbeg;
    let Some((_, c_beg, c_end)) = index.contig_text_range(first) else {
        return;
    };
    // Nie über Contig- oder Stranggrenze hinaus.
    let r0 = r0.max(c_beg as i64) as u64;
    let r1 = r1.min(c_end as i64) as u64;
    let rseq: Vec<u8> = (r0..r1).map(|p| index.text_base(p)).collect();

    // Längste Seeds zuerst; bei gleicher Länge der spätere in der Kette zuerst (wie bwa).
    let mut order: Vec<usize> = (0..chain.seeds.len()).collect();
    order.sort_by_key(|&i| (chain.seeds[i].len, i));
    let mut skipped = vec![false; order.len()];

    for k in (0..order.len()).rev() {
        let s: Seed = chain.seeds[order[k]];
        let (s_qbeg, s_len) = (s.qbeg as i32, s.len as i32);
        // Liegt der Seed schon in einem früheren Bereich und nahe an dessen Diagonale? Geprüft
        // werden wie in bwa alle bisherigen Bereiche des Reads, auch die aus anderen Ketten.
        let near = out.iter().any(|p| {
            if s.rbeg < p.rb
                || s.rbeg + u64::from(s.len) > p.re
                || s_qbeg < p.qb
                || s_qbeg + s_len > p.qe
            {
                return false;
            }
            if f64::from(s_len - p.seedlen0) > 0.1 * f64::from(l_query) {
                return false;
            }
            let qd = s_qbeg - p.qb;
            let rd = (s.rbeg - p.rb) as i32;
            let w = cal_max_gap(sc, qd.min(rd)).min(p.w);
            if qd - rd < w && rd - qd < w {
                return true;
            }
            let qd = p.qe - (s_qbeg + s_len);
            let rd = (p.re - (s.rbeg + u64::from(s.len))) as i32;
            let w = cal_max_gap(sc, qd.min(rd)).min(p.w);
            qd - rd < w && rd - qd < w
        });
        if near {
            // Trotzdem verlängern, wenn ein bereits verlängerter, etwa gleich langer Seed der
            // Kette mit diesem überlappt, aber auf einer anderen Diagonale liegt.
            let conflict = (k + 1..order.len()).any(|i| {
                if skipped[i] {
                    return false;
                }
                let t = chain.seeds[order[i]];
                if f64::from(t.len) < f64::from(s.len) * 0.95 {
                    return false;
                }
                let (tq, tl) = (t.qbeg as i32, t.len as i32);
                let dq = tq - s_qbeg;
                let dr = t.rbeg as i64 - s.rbeg as i64;
                (s_qbeg <= tq && s_qbeg + s_len - tq >= s_len >> 2 && i64::from(dq) != dr)
                    || (tq <= s_qbeg && tq + tl - s_qbeg >= s_len >> 2 && i64::from(dq) != dr)
            });
            if !conflict {
                skipped[k] = true;
                continue;
            }
        }

        let mut a = Region {
            rb: s.rbeg,
            re: s.rbeg + u64::from(s.len),
            qb: 0,
            qe: l_query,
            contig: chain.contig,
            score: -1,
            truesc: -1,
            w: sc.band_width,
            seedcov: 0,
            seedlen0: s_len,
            sub: 0,
            sub_n: 0,
            secondary: -1,
        };
        let mut aw = [sc.band_width; 2];
        // Nach links: Read und Referenz vor dem Seed, beide umgekehrt.
        if s_qbeg > 0 {
            let qs: Vec<u8> = read[..s.qbeg as usize].iter().rev().copied().collect();
            let rs: Vec<u8> = rseq[..(s.rbeg - r0) as usize]
                .iter()
                .rev()
                .copied()
                .collect();
            let (x, w) = extend_retry(
                &qs,
                &rs,
                sc,
                opt.pen_clip5,
                s_len * sc.matching,
                &mut a.score,
            );
            aw[0] = w;
            if x.global_score <= 0 || x.global_score <= a.score - opt.pen_clip5 {
                a.qb = s_qbeg - x.query_len;
                a.rb = s.rbeg - x.target_len as u64;
                a.truesc = a.score;
            } else {
                a.qb = 0;
                a.rb = s.rbeg - x.global_target_len as u64;
                a.truesc = x.global_score;
            }
        } else {
            a.score = s_len * sc.matching;
            a.truesc = a.score;
            a.qb = 0;
            a.rb = s.rbeg;
        }
        // Nach rechts.
        if s_qbeg + s_len != l_query {
            let qe = (s.qbeg + s.len) as usize;
            let re = (s.rbeg + u64::from(s.len) - r0) as usize;
            let sc0 = a.score;
            let (x, w) = extend_retry(
                &read[qe..],
                &rseq[re..],
                sc,
                opt.pen_clip3,
                sc0,
                &mut a.score,
            );
            aw[1] = w;
            if x.global_score <= 0 || x.global_score <= a.score - opt.pen_clip3 {
                a.qe = qe as i32 + x.query_len;
                a.re = r0 + re as u64 + x.target_len as u64;
                a.truesc += a.score - sc0;
            } else {
                a.qe = l_query;
                a.re = r0 + re as u64 + x.global_target_len as u64;
                a.truesc += x.global_score - sc0;
            }
        } else {
            a.qe = l_query;
            a.re = s.rbeg + u64::from(s.len);
        }
        a.seedcov = chain
            .seeds
            .iter()
            .filter(|t| {
                t.qbeg as i32 >= a.qb
                    && (t.qbeg + t.len) as i32 <= a.qe
                    && t.rbeg >= a.rb
                    && t.rbeg + u64::from(t.len) <= a.re
            })
            .map(|t| t.len as i32)
            .sum();
        a.w = aw[0].max(aw[1]);
        out.push(a);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chain::{ChainOptions, FilterOptions, chain_seeds, filter_chains};
    use crate::seed::{SeedOptions, collect_seeds};
    use gc_index::{Reference, Variant};

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn setup() -> (FmIndex, Vec<u8>) {
        let mut rng = 0x00C0_FFEE_1234_5678_u64;
        let s: Vec<u8> = (0..4000)
            .map(|_| b"ACGT"[(xorshift(&mut rng) & 3) as usize])
            .collect();
        let mut r = Reference::new();
        r.push_contig("a", &s[..2500]);
        r.push_contig("b", &s[2500..]);
        let index = FmIndex::build(r, Variant::Kompakt, 0).unwrap().0;
        let codes = s
            .iter()
            .map(|&b| gc_index::reference::encode(b).unwrap())
            .collect();
        (index, codes)
    }

    fn best(index: &FmIndex, read: &[u8]) -> Region {
        let seeds = collect_seeds(index, read, &SeedOptions::default());
        let chains = chain_seeds(&seeds, index.reference().len(), &ChainOptions::default());
        let chains = filter_chains(chains, &FilterOptions::default());
        let mut regions = Vec::new();
        for c in &chains {
            chain_to_regions(index, read, c, &RegionOptions::default(), &mut regions);
        }
        *regions.iter().max_by_key(|r| r.score).expect("ein Bereich")
    }

    #[test]
    fn exact_mismatch_and_gaps_give_the_expected_scores() {
        let (index, g) = setup();
        let l = index.reference().len();
        // exakt
        let read = g[1000..1100].to_vec();
        let r = best(&index, &read);
        assert_eq!((r.rb, r.re, r.qb, r.qe, r.score), (1000, 1100, 0, 100, 100));
        // ein Fehltreffer in der Mitte: 99 Treffer, 4 Abzug
        let mut read = g[1000..1100].to_vec();
        read[50] = (read[50] + 1) & 3;
        let r = best(&index, &read);
        assert_eq!((r.rb, r.re, r.qb, r.qe, r.score), (1000, 1100, 0, 100, 95));
        // Deletion von 3 Basen im Read: 97 Treffer, Lücke 6 + 3
        let mut read = g[1000..1050].to_vec();
        read.extend_from_slice(&g[1053..1100]);
        let r = best(&index, &read);
        assert_eq!((r.rb, r.re, r.qb, r.qe, r.score), (1000, 1100, 0, 97, 88));
        // Insertion von 2 Basen im Read: 100 Treffer, Lücke 6 + 2
        let mut read = g[1000..1050].to_vec();
        read.extend_from_slice(&[(g[1050] + 1) & 3, (g[1050] + 2) & 3]);
        read.extend_from_slice(&g[1050..1100]);
        let r = best(&index, &read);
        assert_eq!((r.rb, r.re, r.qb, r.qe, r.score), (1000, 1100, 0, 102, 92));
        // Gegenstrang: Lage im reversen Komplement des Indextexts
        let read: Vec<u8> = g[3000..3100].iter().rev().map(|c| 3 - c).collect();
        let r = best(&index, &read);
        assert_eq!((r.rb, r.re, r.score), (2 * l - 3100, 2 * l - 3000, 100));
        assert_eq!(r.contig, 1);
    }

    #[test]
    fn read_ends_are_clipped_or_kept_by_the_clip_penalty() {
        let (index, g) = setup();
        // Fehltreffer an der vorletzten Base: bis zum Ende ausrichten kostet 4 - 1 = 3 Punkte,
        // weniger als der Abzug 5 für das Abschneiden; der Read bleibt ganz.
        let mut read = g[1000..1100].to_vec();
        read[98] = (read[98] + 1) & 3;
        let r = best(&index, &read);
        assert_eq!((r.qb, r.qe, r.re), (0, 100, 1100));
        assert_eq!((r.score, r.truesc), (98, 95));
        // Zehn fremde Basen am Ende: werden abgeschnitten.
        let mut read = g[1000..1090].to_vec();
        read.extend(g[1090..1100].iter().map(|c| 3 - c));
        let r = best(&index, &read);
        assert_eq!((r.qb, r.qe, r.rb, r.re, r.score), (0, 90, 1000, 1090, 90));
        // Dasselbe am Anfang.
        let mut read: Vec<u8> = g[1000..1010].iter().map(|c| 3 - c).collect();
        read.extend_from_slice(&g[1010..1100]);
        let r = best(&index, &read);
        assert_eq!((r.qb, r.qe, r.rb, r.re, r.score), (10, 100, 1010, 1100, 90));
    }

    #[test]
    fn regions_never_cross_a_contig_boundary() {
        let (index, g) = setup();
        // Read über die Grenze bei 2500: jeder Bereich bleibt in seinem Contig.
        let read = g[2440..2560].to_vec();
        let seeds = collect_seeds(&index, &read, &SeedOptions::default());
        let chains = chain_seeds(&seeds, index.reference().len(), &ChainOptions::default());
        let mut regions = Vec::new();
        for c in &chains {
            chain_to_regions(&index, &read, c, &RegionOptions::default(), &mut regions);
        }
        assert!(!regions.is_empty());
        for r in &regions {
            assert!(r.re <= 2500 || r.rb >= 2500, "{r:?}");
        }
    }
}
