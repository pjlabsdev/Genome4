//! Einen einzeln ausgerichteten Read als SAM-Zeilen schreiben (BWA-MEM: `mem_reg2sam`,
//! `mem_aln2sam`, `mem_gen_alt`), dazu der ganze Weg vom Read zu seinen Bereichen.
//!
//! Noch ohne Paare und ohne ALT-Contigs.

use crate::align::{Alignment, region_to_alignment};
use crate::chain::{ChainOptions, FilterOptions, chain_seeds, filter_chains};
use crate::fast::{FastMode, check_against_full, fast_region};
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
    /// Schnellweg für leichte Reads ([`crate::fast`]); aus heißt Ausgabe wie bwa.
    pub fast: FastMode,
    /// Rechenkern für die gesammelten Ausrichtungen der Partnersuche (zum Beispiel auf der
    /// Grafikeinheit). Ohne rechnet der Prozessor; die Ausgabe ist in beiden Fällen dieselbe.
    pub rescue_kernel: Option<std::sync::Arc<dyn crate::local::RescueKernel>>,
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
    /// Kennung der Lesegruppe; jede Zeile bekommt sie als `RG` (`-R`).
    pub read_group: Option<String>,
    /// Paare: Angaben zum Partner wie nach `samtools fixmate -m` (`ms` an beiden Hauptzeilen,
    /// `MC` auch an einem Read ohne Ort). Voraussetzung für das Markieren von Duplikaten.
    pub mate_tags: bool,
}

impl Default for SamOptions {
    fn default() -> Self {
        Self {
            min_score: 30,
            drop_ratio: 0.5,
            xa_drop_ratio: 0.8,
            max_xa_hits: 5,
            read_group: None,
            mate_tags: false,
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
    match opt.fast {
        FastMode::Off => find_regions_full(index, read, opt),
        FastMode::On | FastMode::Simple => {
            match fast_region(index, read, opt, opt.fast == FastMode::On) {
                Some(r) => (vec![r], 0.0),
                None => find_regions_full(index, read, opt),
            }
        }
        FastMode::Check => {
            let full = find_regions_full(index, read, opt);
            if let Some(r) = fast_region(index, read, opt, true) {
                check_against_full(&r, &full.0, opt);
            }
            full
        }
    }
}

/// Der volle Weg nach BWA-MEM: Seeds, Ketten, Verlängern, Bereiche bereinigen.
pub fn find_regions_full(index: &FmIndex, read: &[u8], opt: &AlignOptions) -> (Vec<Region>, f32) {
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
pub fn sam_header(
    index: &FmIndex,
    read_group: Option<&ReadGroup>,
    program_version: &str,
    command_line: &str,
) -> String {
    let mut s = String::new();
    for c in &index.reference().contigs {
        let _ = writeln!(s, "@SQ\tSN:{}\tLN:{}", c.name, c.len);
    }
    if let Some(rg) = read_group {
        let _ = writeln!(s, "{}", rg.line);
    }
    let _ = writeln!(
        s,
        "@PG\tID:genome\tPN:genome\tVN:{program_version}\tCL:{command_line}"
    );
    s
}

/// Eine Lesegruppe: die Kopfzeile `@RG` und ihre Kennung.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadGroup {
    /// Die Kopfzeile mit echten Tabulatoren, ohne Zeilenende.
    pub line: String,
    /// Der Wert von `ID`.
    pub id: String,
}

impl ReadGroup {
    /// Liest die Angabe, wie bwa sie mit `-R` nimmt: `@RG\tID:x\tSM:y`, die Tabulatoren als
    /// die zwei Zeichen `\t` geschrieben. Geprüft wird wie in bwa (`bwa_set_rg`), dazu, dass
    /// jedes Feld die Form `XX:Wert` hat und kein Zeilenende vorkommt.
    pub fn parse(text: &str) -> Result<Self, String> {
        if !text.starts_with("@RG") {
            return Err("die Lesegruppe muss mit @RG beginnen".to_string());
        }
        if text.contains('\t') {
            return Err(
                "die Lesegruppe enthält echte Tabulatoren; bitte als \\t schreiben".to_string(),
            );
        }
        let line = text.replace("\\t", "\t");
        if line.contains(['\n', '\r']) || line.contains("\\n") || line.contains("\\r") {
            return Err("die Lesegruppe darf kein Zeilenende enthalten".to_string());
        }
        let mut id = None;
        for field in line.split('\t').skip(1) {
            let bytes = field.as_bytes();
            if bytes.len() < 4 || bytes[2] != b':' {
                return Err(format!(
                    "Feld \"{field}\" der Lesegruppe hat nicht die Form XX:Wert"
                ));
            }
            if &field[..3] == "ID:" && id.is_none() {
                id = Some(field[3..].to_string());
            }
        }
        if line.split('\t').next() != Some("@RG") {
            return Err("nach @RG muss ein \\t folgen".to_string());
        }
        let Some(id) = id else {
            return Err("die Lesegruppe hat kein Feld ID".to_string());
        };
        if id.len() > 255 {
            return Err("die Kennung der Lesegruppe ist länger als 255 Zeichen".to_string());
        }
        Ok(ReadGroup { line, id })
    }
}

/// Summe der Basenqualitäten ab 15 (ASCII wie in FASTQ), wie `samtools fixmate -m` sie als `ms`
/// beim Partner einträgt und `samtools markdup` sie zum Auswählen des besten Paares benutzt.
pub fn quality_score(qual: &[u8]) -> i64 {
    qual.iter()
        .map(|&q| i64::from(q.saturating_sub(33)))
        .filter(|&q| q >= 15)
        .sum()
}

/// Ergänzt bei den beiden Hauptzeilen eines Paares die Angaben zum Partner wie
/// `samtools fixmate -m`: `ms` (Qualitätssumme des Partners) und, wo sie fehlt, `MC` (CIGAR des
/// Partners; hat der Partner keinen Ort, schreibt bwa sie nicht und fixmate schreibt `*`). Zeilen ohne Qualitäten bekommen
/// kein `ms`. Ergänzende Zeilen (Flag 2048) übernehmen Ort, Strang, `MQ` und `MC` von der
/// Hauptzeile des Partners.
pub fn add_mate_tags(lines: &mut [Line]) {
    let main = |lines: &[Line], bit: u32| {
        lines
            .iter()
            .position(|l| l.flag & bit != 0 && l.flag & (0x100 | 0x800) == 0)
    };
    let (Some(a), Some(b)) = (main(lines, 0x40), main(lines, 0x80)) else {
        return;
    };
    for (me, mate) in [(a, b), (b, a)] {
        let mut mate_cigar: String = lines[mate]
            .cigar
            .iter()
            .map(|&(op, len)| format!("{len}{}", op as char))
            .collect();
        if mate_cigar.is_empty() {
            // Partner ohne Ort: fixmate schreibt dann `MC:Z:*`.
            mate_cigar.push('*');
        }
        let score = (!lines[mate].qual.is_empty()).then(|| quality_score(&lines[mate].qual));
        let line = &mut lines[me];
        if !line.tags.iter().any(|(t, _)| t == b"MC") {
            line.tags.push((*b"MC", TagValue::Text(mate_cigar)));
        }
        if let Some(score) = score {
            line.tags.push((*b"ms", TagValue::Int(score)));
        }
    }
    // Ergänzende Zeilen bekommen die Angaben der Hauptzeile des Partners (Ort, Strang, `MQ`,
    // `MC`), wie fixmate es macht. bwa schreibt dort den Ort der eigenen Zeile, wenn der Partner
    // keinen Ort hat, und `MC` mit hartem statt weichem Kürzen.
    for (bit, mate) in [(0x40, b), (0x80, a)] {
        let mate_line = &lines[mate];
        let unmapped = mate_line.flag & 0x4 != 0;
        let (contig, pos, reverse, mapq) = (
            mate_line.contig,
            mate_line.pos,
            mate_line.flag & 0x10 != 0,
            mate_line.mapq,
        );
        let mut cigar: String = mate_line
            .cigar
            .iter()
            .map(|&(op, len)| format!("{len}{}", op as char))
            .collect();
        if cigar.is_empty() {
            cigar.push('*');
        }
        for line in lines
            .iter_mut()
            .filter(|l| l.flag & bit != 0 && l.flag & 0x800 != 0)
        {
            line.mate_contig = contig;
            line.mate_pos = pos;
            line.flag = if reverse {
                line.flag | 0x20
            } else {
                line.flag & !0x20
            };
            if unmapped {
                line.flag |= 0x8;
            } else {
                set_tag(line, *b"MQ", TagValue::Int(i64::from(mapq)));
            }
            set_tag(line, *b"MC", TagValue::Text(cigar.clone()));
        }
    }
}

fn set_tag(line: &mut Line, tag: [u8; 2], value: TagValue) {
    match line.tags.iter_mut().find(|(t, _)| *t == tag) {
        Some(entry) => entry.1 = value,
        None => line.tags.push((tag, value)),
    }
}

impl Line {
    /// Trägt die Lesegruppe ein: `RG` steht wie bei bwa hinter `XS` und vor `SA` und `XA`.
    pub fn set_read_group(&mut self, id: &str) {
        let at = self
            .tags
            .iter()
            .position(|(tag, _)| tag == b"SA" || tag == b"XA")
            .unwrap_or(self.tags.len());
        self.tags
            .insert(at, (*b"RG", TagValue::Text(id.to_string())));
    }
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

/// Wert einer Zusatzangabe.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TagValue {
    Int(i64),
    Text(String),
}

/// Eine fertige Ausgabezeile, unabhängig davon, ob sie als SAM-Text oder als BAM-Satz
/// geschrieben wird.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Line {
    pub name: String,
    pub flag: u32,
    /// `None`: ohne Contig und Position.
    pub contig: Option<u32>,
    /// Ab 0 gezählt; gilt nur mit `contig`.
    pub pos: u64,
    pub mapq: u8,
    /// Buchstabe (`M`, `I`, `D`, `S`, `H`) und Länge; leer heißt `*`.
    pub cigar: Vec<(u8, u32)>,
    pub mate_contig: Option<u32>,
    pub mate_pos: u64,
    pub tlen: i64,
    /// Basen als Codes 0..3, 4 = N, schon in Richtung des Vorwärtsstrangs.
    pub seq: Vec<u8>,
    /// Qualitäten als ASCII wie in FASTQ; leer heißt `*`.
    pub qual: Vec<u8>,
    pub tags: Vec<([u8; 2], TagValue)>,
}

fn cigar_letters(cigar: &[(Op, u32)], hard: bool) -> Vec<(u8, u32)> {
    cigar
        .iter()
        .map(|&(op, n)| {
            let c = if op == Op::SoftClip && hard {
                'H'
            } else {
                op.letter()
            };
            (c as u8, n)
        })
        .collect()
}

fn cigar_text(cigar: &[(u8, u32)]) -> String {
    if cigar.is_empty() {
        return "*".to_string();
    }
    let mut s = String::new();
    for &(c, n) in cigar {
        let _ = write!(s, "{n}{}", c as char);
    }
    s
}

/// Baut Zeile `which` aus `list` (bwa `mem_aln2sam`). `mate` ist der Satz des Partners bei
/// Paaren.
pub(crate) fn make_line(
    index: &FmIndex,
    name: &str,
    read: &[u8],
    qual: &[u8],
    list: &[Rec],
    which: usize,
    mate: Option<&Rec>,
) -> Line {
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
    // Ergänzende Alignments bekommen harte Schnitte: die abgeschnittenen Basen stehen schon in
    // der ersten Zeile.
    let hard = which > 0;
    let mut line = Line {
        name: name.to_string(),
        flag: p.flag,
        contig: p.contig,
        pos: p.pos,
        mapq: if p.contig.is_some() { p.mapq } else { 0 },
        cigar: if p.contig.is_some() {
            cigar_letters(&p.cigar, hard)
        } else {
            Vec::new()
        },
        mate_contig: None,
        mate_pos: 0,
        tlen: 0,
        seq: Vec::new(),
        qual: Vec::new(),
        tags: Vec::new(),
    };
    if let Some((m, c)) = m.as_ref().and_then(|m| m.contig.map(|c| (m, c))) {
        line.mate_contig = Some(c);
        line.mate_pos = m.pos;
        if p.contig == Some(c) && !m.cigar.is_empty() && !p.cigar.is_empty() {
            // Abstand der äußeren Enden: bei einem Read auf dem Gegenstrang zählt sein
            // rechtes Ende.
            let p0 = p.pos as i64 + if p.is_rev { p.ref_len() - 1 } else { 0 };
            let p1 = m.pos as i64 + if m.is_rev { m.ref_len() - 1 } else { 0 };
            line.tlen = -(p0 - p1 + (p0 - p1).signum());
        }
    }
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
    if p.is_rev {
        line.seq = read[qb..qe]
            .iter()
            .rev()
            .map(|&c| if c < 4 { 3 - c } else { 4 })
            .collect();
        if !qual.is_empty() {
            line.qual = qual[qb..qe].iter().rev().copied().collect();
        }
    } else {
        line.seq = read[qb..qe].iter().map(|&c| c.min(4)).collect();
        if !qual.is_empty() {
            line.qual = qual[qb..qe].to_vec();
        }
    }
    let tags = &mut line.tags;
    if !p.cigar.is_empty() {
        tags.push((*b"NM", TagValue::Int(i64::from(p.nm))));
        tags.push((*b"MD", TagValue::Text(p.md.clone())));
    }
    if let Some(m) = &m {
        if !m.cigar.is_empty() {
            // Wie in bwa richtet sich der Schnitt-Buchstabe nach der eigenen Zeile.
            tags.push((
                *b"MC",
                TagValue::Text(cigar_text(&cigar_letters(&m.cigar, hard))),
            ));
        }
        tags.push((*b"MQ", TagValue::Int(i64::from(m.mapq))));
    }
    tags.push((*b"AS", TagValue::Int(i64::from(p.score))));
    if p.sub >= 0 {
        tags.push((*b"XS", TagValue::Int(i64::from(p.sub))));
    }
    if list.len() > 1 {
        let mut sa = String::new();
        for (i, o) in list.iter().enumerate() {
            if i == which {
                continue;
            }
            let Some(c) = o.contig else { continue };
            let strand = if o.is_rev { '-' } else { '+' };
            let _ = write!(sa, "{},{},{strand},", contigs[c as usize].name, o.pos + 1);
            push_cigar(&mut sa, &o.cigar, false);
            let _ = write!(sa, ",{},{};", o.mapq, o.nm);
        }
        tags.push((*b"SA", TagValue::Text(sa)));
    }
    if !p.xa.is_empty() {
        tags.push((*b"XA", TagValue::Text(p.xa.clone())));
    }
    line
}

/// BAM-Behälter für ein Stück `[beg, end)` der Referenz (SAM-Spezifikation, Abschnitt 5.3).
fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    for (shift, offset) in [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> shift == end >> shift {
            return (offset + (beg >> shift)) as u16;
        }
    }
    0
}

impl Line {
    /// Als BAM-Satz (mit Längenfeld vorn), angehängt an `out`. Ganze Zahlen in den
    /// Zusatzangaben bekommen den kleinsten passenden Typ, wie bei samtools.
    pub fn write_bam(&self, out: &mut Vec<u8>) {
        let start = out.len();
        out.extend_from_slice(&[0; 4]);
        let ref_id = self.contig.map_or(-1, |c| c as i32);
        let pos = if self.contig.is_some() {
            self.pos as i32
        } else {
            -1
        };
        let ref_len: i64 = self
            .cigar
            .iter()
            .filter(|(c, _)| matches!(c, b'M' | b'D'))
            .map(|&(_, n)| i64::from(n))
            .sum();
        let bin = reg2bin(i64::from(pos), i64::from(pos) + ref_len.max(1));
        out.extend_from_slice(&ref_id.to_le_bytes());
        out.extend_from_slice(&pos.to_le_bytes());
        out.push((self.name.len() + 1) as u8);
        out.push(self.mapq);
        out.extend_from_slice(&bin.to_le_bytes());
        out.extend_from_slice(&(self.cigar.len() as u16).to_le_bytes());
        out.extend_from_slice(&(self.flag as u16).to_le_bytes());
        out.extend_from_slice(&(self.seq.len() as u32).to_le_bytes());
        let mate_id = self.mate_contig.map_or(-1, |c| c as i32);
        let mate_pos = if self.mate_contig.is_some() {
            self.mate_pos as i32
        } else {
            -1
        };
        out.extend_from_slice(&mate_id.to_le_bytes());
        out.extend_from_slice(&mate_pos.to_le_bytes());
        out.extend_from_slice(&(self.tlen as i32).to_le_bytes());
        out.extend_from_slice(self.name.as_bytes());
        out.push(0);
        for &(c, n) in &self.cigar {
            let op: u32 = match c {
                b'M' => 0,
                b'I' => 1,
                b'D' => 2,
                b'S' => 4,
                _ => 5, // H
            };
            out.extend_from_slice(&(n << 4 | op).to_le_bytes());
        }
        // Basen: zwei je Byte, Codes A=1, C=2, G=4, T=8, N=15.
        const CODE: [u8; 5] = [1, 2, 4, 8, 15];
        for pair in self.seq.chunks(2) {
            let hi = CODE[usize::from(pair[0])];
            let lo = pair.get(1).map_or(0, |&b| CODE[usize::from(b)]);
            out.push(hi << 4 | lo);
        }
        if self.qual.is_empty() {
            out.extend(std::iter::repeat_n(0xFF, self.seq.len()));
        } else {
            out.extend(self.qual.iter().map(|&q| q.saturating_sub(33)));
        }
        for (tag, value) in &self.tags {
            out.extend_from_slice(tag);
            match value {
                TagValue::Text(v) => {
                    out.push(b'Z');
                    out.extend_from_slice(v.as_bytes());
                    out.push(0);
                }
                &TagValue::Int(v) => {
                    if v < 0 {
                        if v >= i64::from(i8::MIN) {
                            out.push(b'c');
                            out.push(v as i8 as u8);
                        } else if v >= i64::from(i16::MIN) {
                            out.push(b's');
                            out.extend_from_slice(&(v as i16).to_le_bytes());
                        } else {
                            out.push(b'i');
                            out.extend_from_slice(&(v as i32).to_le_bytes());
                        }
                    } else if v <= i64::from(u8::MAX) {
                        out.push(b'C');
                        out.push(v as u8);
                    } else if v <= i64::from(u16::MAX) {
                        out.push(b'S');
                        out.extend_from_slice(&(v as u16).to_le_bytes());
                    } else {
                        out.push(b'I');
                        out.extend_from_slice(&(v as u32).to_le_bytes());
                    }
                }
            }
        }
        let size = (out.len() - start - 4) as u32;
        out[start..start + 4].copy_from_slice(&size.to_le_bytes());
    }

    /// Als SAM-Textzeile.
    pub fn write_sam(&self, out: &mut String, index: &FmIndex) {
        let contigs = &index.reference().contigs;
        let _ = write!(out, "{}\t{}\t", self.name, self.flag);
        match self.contig {
            Some(c) => {
                let _ = write!(
                    out,
                    "{}\t{}\t{}\t{}",
                    contigs[c as usize].name,
                    self.pos + 1,
                    self.mapq,
                    cigar_text(&self.cigar)
                );
            }
            None => out.push_str("*\t0\t0\t*"),
        }
        out.push('\t');
        match self.mate_contig {
            Some(c) => {
                if self.contig == Some(c) {
                    out.push('=');
                } else {
                    out.push_str(&contigs[c as usize].name);
                }
                let _ = write!(out, "\t{}\t{}", self.mate_pos + 1, self.tlen);
            }
            None => out.push_str("*\t0\t0"),
        }
        out.push('\t');
        out.extend(self.seq.iter().map(|&c| b"ACGTN"[usize::from(c)] as char));
        out.push('\t');
        if self.qual.is_empty() {
            out.push('*');
        } else {
            out.extend(self.qual.iter().map(|&b| b as char));
        }
        for (tag, value) in &self.tags {
            let (a, b) = (tag[0] as char, tag[1] as char);
            match value {
                TagValue::Int(v) => {
                    let _ = write!(out, "\t{a}{b}:i:{v}");
                }
                TagValue::Text(v) => {
                    let _ = write!(out, "\t{a}{b}:Z:{v}");
                }
            }
        }
        out.push('\n');
    }
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
    for line in read_lines(index, name, read, qual, regions, frac_rep, opt) {
        line.write_sam(out, index);
    }
}

/// Die Ausgabezeilen eines einzelnen Reads.
pub fn read_lines(
    index: &FmIndex,
    name: &str,
    read: &[u8],
    qual: &[u8],
    regions: &[Region],
    frac_rep: f32,
    opt: &AlignOptions,
) -> Vec<Line> {
    let recs = regions_to_recs(index, read, regions, frac_rep, 0, opt);
    let mut lines: Vec<Line> = (0..recs.len())
        .map(|which| make_line(index, name, read, qual, &recs, which, None))
        .collect();
    if let Some(id) = &opt.sam.read_group {
        for line in &mut lines {
            line.set_read_group(id);
        }
    }
    lines
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_group_is_parsed_and_checked_like_bwa() {
        let rg = ReadGroup::parse("@RG\\tID:lauf1\\tSM:probe 1\\tPL:ILLUMINA").unwrap();
        assert_eq!(rg.id, "lauf1");
        assert_eq!(rg.line, "@RG\tID:lauf1\tSM:probe 1\tPL:ILLUMINA");
        // ID muss nicht vorn stehen
        assert_eq!(ReadGroup::parse("@RG\\tSM:x\\tID:7").unwrap().id, "7");
        for bad in [
            "RG\\tID:a",           // beginnt nicht mit @RG
            "@RG\tID:a",           // echter Tabulator
            "@RG\\tSM:x",          // keine ID
            "@RG ID:a",            // kein \t nach @RG
            "@RG\\tID:a\\tfalsch", // Feld ohne XX:
            "@RG\\tID:a\\nSM:x",   // Zeilenende
            "@RG\\tID:",           // leere Kennung
        ] {
            assert!(ReadGroup::parse(bad).is_err(), "{bad:?}");
        }
        let long = format!("@RG\\tID:{}", "x".repeat(256));
        assert!(ReadGroup::parse(&long).is_err());
    }

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

    /// Ergänzende Zeilen tragen nach `add_mate_tags` die Angaben der Hauptzeile des Partners.
    #[test]
    fn supplementary_lines_take_the_mate_primary() {
        let line =
            |flag: u32, pos: u64, cigar: Vec<(u8, u32)>, tags: Vec<([u8; 2], TagValue)>| Line {
                name: "r".into(),
                flag,
                contig: Some(0),
                pos,
                mapq: 60,
                cigar,
                mate_contig: Some(0),
                mate_pos: 0,
                tlen: 0,
                seq: vec![0; 4],
                qual: b"IIII".to_vec(),
                tags,
            };
        let text = |s: &str| TagValue::Text(s.into());
        // Read 1: Hauptzeile vorwärts und eine ergänzende Zeile. Read 2: rückwärts, weich gekürzt.
        let mut lines = vec![
            line(
                0x1 | 0x20 | 0x40,
                100,
                vec![(b'M', 4)],
                vec![(*b"MC", text("1S3M"))],
            ),
            line(
                0x1 | 0x40 | 0x800,
                900,
                vec![(b'H', 2), (b'M', 2)],
                vec![(*b"MC", text("1H3M")), (*b"MQ", TagValue::Int(7))],
            ),
            line(
                0x1 | 0x10 | 0x80,
                300,
                vec![(b'S', 1), (b'M', 3)],
                vec![(*b"MC", text("4M"))],
            ),
        ];
        lines[1].mate_pos = 900;
        add_mate_tags(&mut lines);
        let sup = &lines[1];
        assert_eq!((sup.mate_contig, sup.mate_pos), (Some(0), 300));
        assert_eq!(sup.flag, 0x1 | 0x20 | 0x40 | 0x800);
        assert!(sup.tags.contains(&(*b"MC", text("1S3M"))));
        assert!(sup.tags.contains(&(*b"MQ", TagValue::Int(60))));
        assert!(!sup.tags.iter().any(|(t, _)| t == b"ms"));

        // Partner ohne Ort: Bit 8 gesetzt, Bit 32 nach dem Strang seiner Hauptzeile, `MC:Z:*`.
        let mut lines = vec![
            line(0x1 | 0x8 | 0x40, 100, vec![(b'M', 4)], vec![]),
            line(
                0x1 | 0x8 | 0x20 | 0x40 | 0x800,
                900,
                vec![(b'H', 2), (b'M', 2)],
                vec![],
            ),
            line(0x1 | 0x4 | 0x80, 100, vec![], vec![]),
        ];
        add_mate_tags(&mut lines);
        let sup = &lines[1];
        assert_eq!(sup.mate_pos, 100);
        assert_eq!(sup.flag, 0x1 | 0x8 | 0x40 | 0x800);
        assert!(sup.tags.contains(&(*b"MC", text("*"))));
        assert!(!sup.tags.iter().any(|(t, _)| t == b"MQ"));
    }
}
