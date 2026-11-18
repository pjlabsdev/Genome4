//! Aus einem Bereich das fertige Alignment machen (BWA-MEM: `mem_reg2aln`): Ort auf dem
//! Vorwärtsstrang, CIGAR mit abgeschnittenen Read-Enden, `NM` und `MD`.

use crate::extend::Scoring;
use crate::global::{Cigar, Op, gen_cigar};
use crate::region::Region;
use gc_index::FmIndex;

/// Ein Alignment, wie es in eine SAM-Zeile geht.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub contig: u32,
    /// Linkes Ende im Contig, ab 0 gezählt.
    pub pos: u64,
    pub is_reverse: bool,
    /// In Richtung des Vorwärtsstrangs, mit `S` für abgeschnittene Enden.
    pub cigar: Cigar,
    pub nm: u32,
    pub md: String,
    /// Bester lokaler Score des Bereichs (`AS`).
    pub score: i32,
    /// Score des zweitbesten Bereichs (`XS`).
    pub sub: i32,
    pub secondary: bool,
}

/// Nötige Bandbreite, geschätzt aus Längen und Score (bwa `infer_bw`).
fn infer_bw(l1: i32, l2: i32, score: i32, a: i32, q: i32, r: i32) -> i32 {
    if l1 == l2 && l1 * a - score < (q + r - a) << 1 {
        return 0;
    }
    let w = (f64::from(l1.min(l2) * a - score - q) / f64::from(r) + 2.0) as i32;
    w.max((l1 - l2).abs())
}

/// Rechnet das Alignment eines Bereichs aus. Das Band beginnt so schmal, wie der Score des
/// Bereichs es nahelegt, und wird bis zu zweimal verdoppelt, solange das globale Alignment
/// schlechter ausfällt als das lokale. `None`, wenn der Bereich kein gültiges Stück Referenz
/// beschreibt.
pub fn region_to_alignment(
    index: &FmIndex,
    read: &[u8],
    r: &Region,
    sc: &Scoring,
) -> Option<Alignment> {
    let l = index.reference().len();
    let l_query = read.len() as i32;
    let (ql, rl) = (r.qe - r.qb, (r.re - r.rb) as i32);
    let mut w2 = infer_bw(ql, rl, r.truesc, sc.matching, sc.o_del, sc.e_del).max(infer_bw(
        ql,
        rl,
        r.truesc,
        sc.matching,
        sc.o_ins,
        sc.e_ins,
    ));
    if w2 > sc.band_width {
        w2 = w2.min(r.w);
    }
    let query = &read[r.qb as usize..r.qe as usize];
    let max_w = sc.band_width << 2;
    let mut last_score = -(1 << 30);
    let mut tries = 0;
    let g = loop {
        w2 = w2.min(max_w);
        let g = gen_cigar(index, query, r.rb, r.re, sc, w2)?;
        if g.score == last_score || w2 == max_w {
            break g;
        }
        last_score = g.score;
        w2 <<= 1;
        tries += 1;
        if !(tries < 3 && g.score < r.truesc - sc.matching) {
            break g;
        }
    };
    let is_reverse = r.rb >= l;
    let mut pos = if is_reverse { 2 * l - r.re } else { r.rb };
    let mut cigar = g.cigar;
    // Eine Deletion am Rand gehört nicht zum Alignment.
    if let Some(&(Op::Del, n)) = cigar.first() {
        pos += u64::from(n);
        cigar.remove(0);
    } else if let Some(&(Op::Del, _)) = cigar.last() {
        cigar.pop();
    }
    if r.qb != 0 || r.qe != l_query {
        let (clip5, clip3) = if is_reverse {
            (l_query - r.qe, r.qb)
        } else {
            (r.qb, l_query - r.qe)
        };
        if clip5 > 0 {
            cigar.insert(0, (Op::SoftClip, clip5 as u32));
        }
        if clip3 > 0 {
            cigar.push((Op::SoftClip, clip3 as u32));
        }
    }
    let offset = index.reference().contigs[r.contig as usize].offset;
    Some(Alignment {
        contig: r.contig,
        pos: pos - offset,
        is_reverse,
        cigar,
        nm: g.nm,
        md: g.md,
        score: r.score,
        sub: r.sub,
        secondary: r.secondary >= 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chain::{ChainOptions, FilterOptions, chain_seeds, filter_chains};
    use crate::finalize::{FinalizeOptions, mark_primary, sort_dedup};
    use crate::global::cigar_string;
    use crate::region::{RegionOptions, chain_to_regions};
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
        r.push_contig("a", &s[..1000]);
        r.push_contig("b", &s[1000..]);
        let index = FmIndex::build(r, Variant::Kompakt, 0).unwrap().0;
        let codes = s
            .iter()
            .map(|&b| gc_index::reference::encode(b).unwrap())
            .collect();
        (index, codes)
    }

    fn align(index: &FmIndex, read: &[u8]) -> Alignment {
        let seeds = collect_seeds(index, read, &SeedOptions::default());
        let chains = chain_seeds(&seeds, index.reference().len(), &ChainOptions::default());
        let chains = filter_chains(chains, &FilterOptions::default());
        let ro = RegionOptions::default();
        let mut regions = Vec::new();
        for c in &chains {
            chain_to_regions(index, read, c, &ro, &mut regions);
        }
        let fo = FinalizeOptions::default();
        sort_dedup(&mut regions, index, read, &ro.scoring, &fo);
        mark_primary(&mut regions, 0, &ro.scoring, &fo);
        region_to_alignment(index, read, &regions[0], &ro.scoring).expect("Alignment")
    }

    fn revcomp(read: &[u8]) -> Vec<u8> {
        read.iter().rev().map(|&b| 3 - b).collect()
    }

    fn summary(a: &Alignment) -> (u32, u64, bool, String, u32, String) {
        (
            a.contig,
            a.pos,
            a.is_reverse,
            cigar_string(&a.cigar),
            a.nm,
            a.md.clone(),
        )
    }

    #[test]
    fn exact_mismatch_deletion_and_insertion_on_both_strands() {
        let (index, codes) = setup();
        let letter = |c: u8| b"ACGT"[c as usize] as char;
        // Contig b beginnt bei 1000; Read ab 1500, also Position 500 in b.
        let exact = codes[1500..1600].to_vec();
        let want = (1, 500, false, "100M".to_string(), 0, "100".to_string());
        assert_eq!(summary(&align(&index, &exact)), want);
        let mut rc = summary(&align(&index, &revcomp(&exact)));
        rc.2 = false;
        assert_eq!(rc, want, "Gegenstrang: gleicher Ort, gleiches CIGAR");

        let mut mm = exact.clone();
        mm[40] = (mm[40] + 1) & 3;
        let md = format!("40{}59", letter(codes[1540]));
        for (read, rev) in [(mm.clone(), false), (revcomp(&mm), true)] {
            let a = align(&index, &read);
            assert_eq!(
                summary(&a),
                (1, 500, rev, "100M".to_string(), 1, md.clone())
            );
        }

        // 3 Basen der Referenz fehlen im Read (1550..1553).
        let del: Vec<u8> = [&codes[1500..1550], &codes[1553..1603]].concat();
        let deleted: String = codes[1550..1553].iter().map(|&c| letter(c)).collect();
        for (read, rev) in [(del.clone(), false), (revcomp(&del), true)] {
            let a = align(&index, &read);
            assert_eq!((a.contig, a.pos, a.is_reverse, a.nm), (1, 500, rev, 3));
            let ops: Vec<Op> = a.cigar.iter().map(|x| x.0).collect();
            assert_eq!(ops, [Op::Match, Op::Del, Op::Match]);
            assert_eq!(a.cigar[1].1, 3);
            // Die Lücke darf in einer Wiederholung nach links rutschen; ohne Wiederholung
            // steht sie genau an der Stelle.
            if a.cigar[0].1 == 50 {
                assert_eq!(a.md, format!("50^{deleted}50"));
            }
        }

        // 2 zusätzliche Basen im Read.
        let mut ins = exact.clone();
        ins.splice(60..60, [(codes[1560] + 1) & 3, (codes[1560] + 2) & 3]);
        for (read, rev) in [(ins.clone(), false), (revcomp(&ins), true)] {
            let a = align(&index, &read);
            assert_eq!((a.contig, a.pos, a.is_reverse, a.nm), (1, 500, rev, 2));
            assert_eq!(a.md, "100");
            let ops: Vec<Op> = a.cigar.iter().map(|x| x.0).collect();
            assert_eq!(ops, [Op::Match, Op::Ins, Op::Match]);
        }
    }

    #[test]
    fn foreign_read_ends_become_soft_clips_in_forward_direction() {
        let (index, codes) = setup();
        // 12 fremde Basen am Anfang des Reads.
        let mut read: Vec<u8> = (0..12).map(|i| (codes[1500 + i] + 2) & 3).collect();
        read.extend_from_slice(&codes[1512..1600]);
        let a = align(&index, &read);
        assert_eq!((a.pos, cigar_string(&a.cigar).as_str()), (512, "12S88M"));
        // Derselbe Read vom Gegenstrang: das abgeschnittene Stück steht in Vorwärtsrichtung
        // weiter links, also wieder vorn.
        let a = align(&index, &revcomp(&read));
        assert_eq!(
            (a.pos, a.is_reverse, cigar_string(&a.cigar).as_str()),
            (512, true, "12S88M")
        );
    }
}
