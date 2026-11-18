//! Einen einzeln ausgerichteten Read als SAM-Zeilen schreiben (BWA-MEM: `mem_reg2sam`,
//! `mem_aln2sam`, `mem_gen_alt`), dazu der ganze Weg vom Read zu seinen Bereichen.
//!
//! Noch ohne Paare und ohne ALT-Contigs.

use crate::align::{Alignment, region_to_alignment};
use crate::chain::{ChainOptions, FilterOptions, chain_seeds, filter_chains};
use crate::finalize::{FinalizeOptions, mapq, mark_primary, sort_dedup};
use crate::global::{Cigar, Op};
use crate::pair::PairOptions;
use crate::region::{Region, RegionOptions, chain_to_regions};
use crate::seed::{SeedOptions, collect_intervals, locate_seeds, repeat_fraction};
use gc_index::FmIndex;
use std::fmt::Write;

/// Alle Einstellungen des Ausrichtens. Die Vorgaben sind die von BWA-MEM.
#[derive(Debug, Clone, Default)]
pub struct AlignOptions {
    pub seed: SeedOptions,
    pub chain: ChainOptions,
    pub filter: FilterOptions,
    pub region: RegionOptions,
    pub finalize: FinalizeOptions,
    pub sam: SamOptions,
    pub pair: PairOptions,
}

/// Einstellungen der Ausgabe.
#[derive(Debug, Clone, PartialEq)]
pub struct SamOptions {
    /// Bereiche unter diesem Score werden nicht gemeldet (`-T`).
    pub min_score: i32,
    /// Ein untergeordneter Bereich unter diesem Anteil des Scores seines übergeordneten wird
    /// auch mit "alle melden" nicht gemeldet (`-D`).
    pub drop_ratio: f32,
    /// Untergeordnete Bereiche ab diesem Anteil des Scores kommen in `XA`. Bewusst mit
    /// einfacher Genauigkeit wie in bwa: 0,8 ist dort etwas größer als 0,8, ein Bereich mit
    /// genau 80 % des Scores zählt deshalb nicht.
    pub xa_drop_ratio: f32,
    /// Gibt es mehr solche Bereiche als hier, entfällt `XA` ganz (`-h`).
    pub max_xa_hits: usize,
}

impl Default for SamOptions {
    fn default() -> Self {
        Self {
            min_score: 30,
            drop_ratio: 0.5,
            xa_drop_ratio: 0.8,
            max_xa_hits: 5,
        }
    }
}

/// Vom Read (Codes 0..3, sonst N) zu seinen bereinigten, nach Rang geordneten Bereichen.
/// `read_id` ist die laufende Nummer des Reads; sie entscheidet zwischen gleich guten Orten.
/// Zweiter Wert: Anteil des Reads in sehr häufigen Seeds.
pub fn align_read(
    index: &FmIndex,
    read: &[u8],
    read_id: u64,
    opt: &AlignOptions,
) -> (Vec<Region>, f32) {
    let (mut regions, frac_rep) = find_regions(index, read, opt);
    mark_primary(&mut regions, read_id, &opt.region.scoring, &opt.finalize);
    (regions, frac_rep)
}

/// Wie [`align_read`], aber ohne die Rangfolge festzulegen: Die Bereiche sind bereinigt und
/// nach Score, dann Lage geordnet (bwa `mem_align1_core`). So gehen sie bei Paaren in die
/// Schätzung der Insertgröße und in die Partnersuche.
pub fn find_regions(index: &FmIndex, read: &[u8], opt: &AlignOptions) -> (Vec<Region>, f32) {
    let l = index.reference().len();
    let intervals = collect_intervals(index, read, &opt.seed);
    let frac_rep = repeat_fraction(&intervals, opt.seed.max_occ, read.len());
    let seeds = locate_seeds(index, &intervals, &opt.seed);
    let chains = filter_chains(chain_seeds(&seeds, l, &opt.chain), &opt.filter);
    let mut regions = Vec::new();
    for c in &chains {
        chain_to_regions(index, read, c, &opt.region, &mut regions);
    }
    let sc = &opt.region.scoring;
    sort_dedup(&mut regions, index, Some(read), sc, &opt.finalize);
    (regions, frac_rep)
}

/// Kopfzeilen: je Contig eine `@SQ`-Zeile, dazu `@PG`.
pub fn sam_header(index: &FmIndex, program_version: &str, command_line: &str) -> String {
    let mut s = String::new();
    for c in &index.reference().contigs {
        let _ = writeln!(s, "@SQ\tSN:{}\tLN:{}", c.name, c.len);
    }
    let _ = writeln!(
        s,
        "@PG\tID:genome\tPN:genome\tVN:{program_version}\tCL:{command_line}"
    );
    s
}

/// Ein Ausgabesatz, wie er in eine SAM-Zeile geht (bwa `mem_aln_t`).
#[derive(Debug, Clone)]
pub(crate) struct Rec {
    /// `None`: nicht ausgerichtet.
    pub contig: Option<u32>,
    pub pos: u64,
    pub is_rev: bool,
    /// Leer: kein Alignment (auch bei einem Read, der nur den Ort seines Partners trägt).
    pub cigar: Cigar,
    pub nm: u32,
    pub md: String,
    pub score: i32,
    /// -1: kein `XS`.
    pub sub: i32,
    pub mapq: u8,
    pub flag: u32,
    pub xa: String,
}

impl Rec {
    pub(crate) fn unmapped() -> Rec {
        Rec {
            contig: None,
            pos: 0,
            is_rev: false,
            cigar: Cigar::new(),
            nm: 0,
            md: String::new(),
            score: 0,
            sub: 0,
            mapq: 0,
            flag: 0x4,
            xa: String::new(),
        }
    }

    pub(crate) fn mapped(a: Alignment, mapq: u8) -> Rec {
        Rec {
            contig: Some(a.contig),
            pos: a.pos,
            is_rev: a.is_reverse,
            cigar: a.cigar,
            nm: a.nm,
            md: a.md,
            score: a.score,
            sub: a.sub,
            mapq,
            flag: if a.secondary { 0x100 } else { 0 },
            xa: String::new(),
        }
    }

    /// Länge auf der Referenz.
    fn ref_len(&self) -> i64 {
        self.cigar
            .iter()
            .filter(|(op, _)| matches!(op, Op::Match | Op::Del))
            .map(|&(_, n)| i64::from(n))
            .sum()
    }
}

/// CIGAR schreiben; `*`, wenn es keinen gibt. In ergänzenden Zeilen (`hard`) stehen harte
/// statt weicher Schnitte.
fn push_cigar(out: &mut String, cigar: &[(Op, u32)], hard: bool) {
    if cigar.is_empty() {
        out.push('*');
    }
    for &(op, n) in cigar {
        let c = if op == Op::SoftClip && hard {
            'H'
        } else {
            op.letter()
        };
        let _ = write!(out, "{n}{c}");
    }
}

/// `XA` je Bereich: die untergeordneten Bereiche, die fast so gut sind wie ihr übergeordneter
/// (bwa `mem_gen_alt`). `parents[i]` ist der Bereich, dem `i` untergeordnet ist, oder negativ.
pub(crate) fn gen_xa(
    index: &FmIndex,
    read: &[u8],
    regions: &[Region],
    parents: &[i32],
    opt: &AlignOptions,
) -> Vec<String> {
    let sc = &opt.region.scoring;
    let parent = |i: usize| -> Option<usize> {
        let k = parents[i];
        (k >= 0
            && f64::from(regions[i].score)
                >= f64::from(regions[k as usize].score) * f64::from(opt.sam.xa_drop_ratio))
        .then_some(k as usize)
    };
    let mut count = vec![0usize; regions.len()];
    for i in 0..regions.len() {
        if let Some(k) = parent(i) {
            count[k] += 1;
        }
    }
    let mut xa = vec![String::new(); regions.len()];
    for (i, region) in regions.iter().enumerate() {
        let Some(k) = parent(i) else { continue };
        if count[k] > opt.sam.max_xa_hits {
            continue;
        }
        let Some(a) = region_to_alignment(index, read, region, sc) else {
            continue;
        };
        let contig = &index.reference().contigs[a.contig as usize].name;
        let strand = if a.is_reverse { '-' } else { '+' };
        let _ = write!(xa[k], "{contig},{strand}{},", a.pos + 1);
        push_cigar(&mut xa[k], &a.cigar, false);
        let _ = write!(xa[k], ",{};", a.nm);
    }
    xa
}

/// Die Sätze eines Reads (bwa `mem_reg2sam`): der beste Bereich und jeder weitere, der einen
/// anderen Teil des Reads abdeckt (ergänzendes Alignment). Ohne brauchbaren Bereich ein Satz
/// "nicht ausgerichtet". `extra_flag` kommt in jeden Satz.
pub(crate) fn regions_to_recs(
    index: &FmIndex,
    read: &[u8],
    regions: &[Region],
    frac_rep: f32,
    extra_flag: u32,
    opt: &AlignOptions,
) -> Vec<Rec> {
    let sc = &opt.region.scoring;
    let parents: Vec<i32> = regions.iter().map(|r| r.secondary).collect();
    let mut xa = gen_xa(index, read, regions, &parents, opt);
    let mut recs: Vec<Rec> = Vec::new();
    for (k, r) in regions.iter().enumerate() {
        if r.score < opt.sam.min_score || r.secondary >= 0 {
            continue;
        }
        let Some(aln) = region_to_alignment(index, read, r, sc) else {
            continue;
        };
        let mut rec = Rec::mapped(aln, mapq(r, frac_rep, sc, &opt.finalize));
        rec.flag |= extra_flag;
        if let Some(first) = recs.first() {
            rec.flag |= 0x800;
            rec.mapq = rec.mapq.min(first.mapq);
        }
        rec.xa = std::mem::take(&mut xa[k]);
        recs.push(rec);
    }
    if recs.is_empty() {
        let mut rec = Rec::unmapped();
        rec.flag |= extra_flag;
        recs.push(rec);
    }
    recs
}

/// Schreibt Satz `which` aus `list` als SAM-Zeile (bwa `mem_aln2sam`). `mate` ist der Satz des
/// Partners bei Paaren.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_rec(
    out: &mut String,
    index: &FmIndex,
    name: &str,
    read: &[u8],
    qual: &[u8],
    list: &[Rec],
    which: usize,
    mate: Option<&Rec>,
) {
    let mut p = list[which].clone();
    let mut m = mate.cloned();
    if let Some(m) = m.as_mut() {
        p.flag |= 0x1;
        if m.contig.is_none() {
            p.flag |= 0x8;
        }
        // Ein Read ohne Ort steht beim Ort seines Partners, und umgekehrt.
        if p.contig.is_none() && m.contig.is_some() {
            p.contig = m.contig;
            p.pos = m.pos;
            p.is_rev = m.is_rev;
            p.cigar.clear();
        }
        if m.contig.is_none() && p.contig.is_some() {
            m.contig = p.contig;
            m.pos = p.pos;
            m.is_rev = p.is_rev;
            m.cigar.clear();
        }
        if m.is_rev {
            p.flag |= 0x20;
        }
    }
    if p.is_rev {
        p.flag |= 0x10;
    }
    let contigs = &index.reference().contigs;
    let _ = write!(out, "{name}\t{}\t", p.flag);
    // Ergänzende Alignments bekommen harte Schnitte: die abgeschnittenen Basen stehen schon in
    // der ersten Zeile.
    let hard = which > 0;
    match p.contig {
        Some(c) => {
            let _ = write!(
                out,
                "{}\t{}\t{}\t",
                contigs[c as usize].name,
                p.pos + 1,
                p.mapq
            );
            push_cigar(out, &p.cigar, hard);
        }
        None => out.push_str("*\t0\t0\t*"),
    }
    out.push('\t');
    match m.as_ref().and_then(|m| m.contig.map(|c| (m, c))) {
        Some((m, c)) => {
            if p.contig == Some(c) {
                out.push('=');
            } else {
                out.push_str(&contigs[c as usize].name);
            }
            let _ = write!(out, "\t{}\t", m.pos + 1);
            if p.contig == Some(c) && !m.cigar.is_empty() && !p.cigar.is_empty() {
                // Abstand der äußeren Enden: bei einem Read auf dem Gegenstrang zählt sein
                // rechtes Ende.
                let p0 = p.pos as i64 + if p.is_rev { p.ref_len() - 1 } else { 0 };
                let p1 = m.pos as i64 + if m.is_rev { m.ref_len() - 1 } else { 0 };
                let _ = write!(out, "{}", -(p0 - p1 + (p0 - p1).signum()));
            } else {
                out.push('0');
            }
        }
        None => out.push_str("*\t0\t0"),
    }
    out.push('\t');
    // Anteil des Reads, der in dieser Zeile steht, in Leserichtung des Reads.
    let clip = |end: Option<&(Op, u32)>| -> usize {
        match end {
            Some(&(Op::SoftClip, n)) if hard => n as usize,
            _ => 0,
        }
    };
    let (first, last) = (clip(p.cigar.first()), clip(p.cigar.last()));
    let (qb, qe) = if p.is_rev {
        (last, read.len() - first)
    } else {
        (first, read.len() - last)
    };
    let letters = |codes: &mut dyn Iterator<Item = u8>, out: &mut String, table: &[u8; 5]| {
        out.extend(codes.map(|c| table[usize::from(c.min(4))] as char));
    };
    if p.is_rev {
        letters(&mut read[qb..qe].iter().rev().copied(), out, b"TGCAN");
    } else {
        letters(&mut read[qb..qe].iter().copied(), out, b"ACGTN");
    }
    out.push('\t');
    if qual.is_empty() {
        out.push('*');
    } else if p.is_rev {
        out.extend(qual[qb..qe].iter().rev().map(|&b| b as char));
    } else {
        out.extend(qual[qb..qe].iter().map(|&b| b as char));
    }
    if !p.cigar.is_empty() {
        let _ = write!(out, "\tNM:i:{}\tMD:Z:{}", p.nm, p.md);
    }
    if let Some(m) = &m {
        if !m.cigar.is_empty() {
            out.push_str("\tMC:Z:");
            // Wie in bwa richtet sich der Schnitt-Buchstabe nach der eigenen Zeile.
            push_cigar(out, &m.cigar, hard);
        }
        let _ = write!(out, "\tMQ:i:{}", m.mapq);
    }
    let _ = write!(out, "\tAS:i:{}", p.score);
    if p.sub >= 0 {
        let _ = write!(out, "\tXS:i:{}", p.sub);
    }
    if list.len() > 1 {
        out.push_str("\tSA:Z:");
        for (i, o) in list.iter().enumerate() {
            if i == which {
                continue;
            }
            let Some(c) = o.contig else { continue };
            let strand = if o.is_rev { '-' } else { '+' };
            let _ = write!(out, "{},{},{strand},", contigs[c as usize].name, o.pos + 1);
            push_cigar(out, &o.cigar, false);
            let _ = write!(out, ",{},{};", o.mapq, o.nm);
        }
    }
    if !p.xa.is_empty() {
        let _ = write!(out, "\tXA:Z:{}", p.xa);
    }
    out.push('\n');
}

/// Schreibt die SAM-Zeilen eines einzelnen Reads nach `out`. `regions` kommt aus
/// [`align_read`]. `qual` leer heißt: keine Qualitätswerte (`*`).
#[allow(clippy::too_many_arguments)]
pub fn write_sam(
    out: &mut String,
    index: &FmIndex,
    name: &str,
    read: &[u8],
    qual: &[u8],
    regions: &[Region],
    frac_rep: f32,
    opt: &AlignOptions,
) {
    let recs = regions_to_recs(index, read, regions, frac_rep, 0, opt);
    for which in 0..recs.len() {
        write_rec(out, index, name, read, qual, &recs, which, None);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gc_index::{Reference, Variant};

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn setup() -> (FmIndex, Vec<u8>, Vec<u8>) {
        let mut rng = 0x00C0_FFEE_1234_5678_u64;
        let s: Vec<u8> = (0..6000)
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
        (index, codes, s)
    }

    fn sam(index: &FmIndex, read: &[u8], qual: &[u8]) -> Vec<Vec<String>> {
        let opt = AlignOptions::default();
        let (regions, frac_rep) = align_read(index, read, 0, &opt);
        let mut out = String::new();
        write_sam(&mut out, index, "r1", read, qual, &regions, frac_rep, &opt);
        out.lines()
            .map(|l| l.split('\t').map(str::to_string).collect())
            .collect()
    }

    fn revcomp(read: &[u8]) -> Vec<u8> {
        read.iter().rev().map(|&b| 3 - b).collect()
    }

    #[test]
    fn a_plain_read_gives_one_line_on_either_strand() {
        let (index, codes, letters) = setup();
        let read = codes[1500..1600].to_vec();
        let text = String::from_utf8(letters[1500..1600].to_vec()).unwrap();
        let qual: Vec<u8> = (0..100).map(|i| b'!' + (i % 40) as u8).collect();
        let lines = sam(&index, &read, &qual);
        assert_eq!(lines.len(), 1);
        let want = [
            "r1",
            "0",
            "b",
            "501",
            "60",
            "100M",
            "*",
            "0",
            "0",
            &text,
            std::str::from_utf8(&qual).unwrap(),
            "NM:i:0",
            "MD:Z:100",
            "AS:i:100",
            "XS:i:0",
        ];
        assert_eq!(lines[0], want);
        // Gegenstrang: Flag 16, Basen und Qualitäten wieder in Vorwärtsrichtung.
        let rq: Vec<u8> = qual.iter().rev().copied().collect();
        let lines = sam(&index, &revcomp(&read), &rq);
        let mut want_rev = want.map(str::to_string);
        want_rev[1] = "16".into();
        assert_eq!(lines[0], want_rev);
    }

    #[test]
    fn a_read_without_a_place_is_written_as_unmapped() {
        let (index, _, _) = setup();
        let read = vec![0u8; 60]; // 60 mal A kommt in der Zufallsfolge nicht vor
        let lines = sam(&index, &read, &[]);
        assert_eq!(lines.len(), 1);
        assert_eq!(&lines[0][1..9], ["4", "*", "0", "0", "*", "*", "0", "0"]);
        assert_eq!(lines[0][10], "*");
        assert_eq!(&lines[0][11..], ["AS:i:0", "XS:i:0"]);
    }

    #[test]
    fn a_split_read_gives_a_primary_and_a_supplementary_line_that_point_to_each_other() {
        let (index, codes, _) = setup();
        // 90 Basen von einer Stelle, 60 von einer anderen, gut 2.000 Basen weiter. Die Naht
        // so gewählt, dass keines der beiden Stücke zufällig über sie hinaus passt.
        let j = (3600..3700)
            .find(|&j| codes[j] != codes[1590] && codes[j - 1] != codes[1589])
            .unwrap();
        let read: Vec<u8> = [&codes[1500..1590], &codes[j..j + 60]].concat();
        let lines = sam(&index, &read, &[]);
        assert_eq!(lines.len(), 2);
        let (p, s) = (&lines[0], &lines[1]);
        let second = (j - 1000 + 1).to_string();
        assert_eq!((&p[1][..], &p[3][..], &p[5][..]), ("0", "501", "90M60S"));
        assert_eq!(p[9].len(), 150);
        // Ergänzende Zeile: Flag 2048, harter Schnitt, nur die 60 eigenen Basen.
        assert_eq!(
            (&s[1][..], &s[3][..], &s[5][..]),
            ("2048", &second[..], "90H60M")
        );
        assert_eq!(s[9].len(), 60);
        let tag = |l: &Vec<String>, t: &str| l.iter().find(|f| f.starts_with(t)).cloned().unwrap();
        let mq = |l: &Vec<String>| l[4].clone();
        assert_eq!(
            tag(p, "SA:Z:"),
            format!("SA:Z:b,{second},+,90S60M,{},0;", mq(s))
        );
        assert_eq!(tag(s, "SA:Z:"), format!("SA:Z:b,501,+,90M60S,{},0;", mq(p)));
        assert!(mq(s).parse::<u8>().unwrap() <= mq(p).parse::<u8>().unwrap());
    }

    #[test]
    fn two_equally_good_places_give_mapq_0_and_the_other_place_in_xa() {
        let mut rng = 0x1357_9BDF_2468_ACE0_u64;
        let mut s: Vec<u8> = (0..5000)
            .map(|_| b"ACGT"[(xorshift(&mut rng) & 3) as usize])
            .collect();
        let copy = s[500..600].to_vec();
        s[3000..3100].copy_from_slice(&copy);
        let mut r = Reference::new();
        r.push_contig("a", &s);
        let index = FmIndex::build(r, Variant::Kompakt, 0).unwrap().0;
        let read: Vec<u8> = copy
            .iter()
            .map(|&b| gc_index::reference::encode(b).unwrap())
            .collect();
        let lines = sam(&index, &read, &[]);
        assert_eq!(lines.len(), 1);
        let l = &lines[0];
        assert_eq!((&l[4][..], &l[5][..]), ("0", "100M"));
        let (here, there) = if l[3] == "501" {
            ("501", "3001")
        } else {
            ("3001", "501")
        };
        assert_eq!(l[3], here);
        assert!(l.contains(&"XS:i:100".to_string()));
        assert!(l.contains(&format!("XA:Z:a,+{there},100M,0;")));
    }
}
