//! Qualitätszahlen einer nach Ort sortierten BAM in einem Durchlauf über die rohen Datensätze:
//! Zählung nach Flags, Tags, MAPQ, Insertgröße, Read-Länge, Basenqualität, Fehlerrate und
//! Abdeckung je Base (Verteilung, je Contig und in Fächern entlang des Contigs).
//!
//! Die Abdeckung entsteht ohne Feld über das ganze Genom: Jeder Datensatz trägt an Anfang und Ende
//! seiner ausgerichteten Stücke +1 und −1 in einen Ring ein. Weil die Datensätze nach Ort
//! kommen, sind alle Stellen vor dem Anfang des aktuellen Datensatzes fertig und werden
//! ausgezählt. Ein Ring von 4 Mio. Stellen reicht für jedes Alignment aus kurzen Reads.

use crate::flagstat::FlagStats;

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const SECONDARY: u16 = 0x100;
const QC_FAIL: u16 = 0x200;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// Größte Insertgröße mit eigenem Fach; größere landen im letzten Fach.
pub const MAX_INSERT: usize = 2000;
/// Größte Read-Länge mit eigenem Fach.
pub const MAX_READ_LENGTH: usize = 1000;
/// Größte Abdeckung mit eigenem Fach in der Verteilung.
pub const MAX_DEPTH: usize = 1000;
/// Breite der Fächer entlang eines Contigs.
pub const BIN: u64 = 100_000;
/// Schwellen, für die je Contig der Anteil der Stellen mit mindestens so viel Abdeckung zählt.
pub const DEPTH_STEPS: [u64; 4] = [1, 10, 20, 30];

const RING_BITS: u32 = 22;
const RING: u64 = 1 << RING_BITS;
const MASK: u64 = RING - 1;

/// Zahlen je Contig.
#[derive(Debug, Clone, Default)]
pub struct ContigQc {
    pub name: String,
    pub length: u64,
    /// Hauptzeilen mit Ort auf diesem Contig, Duplikate eingeschlossen.
    pub reads: u64,
    /// Davon als Duplikat markiert.
    pub duplicates: u64,
    /// Ausgerichtete Basen (M, =, X) der Sätze, die zur Abdeckung zählen.
    pub bases: u64,
    /// Stellen mit Abdeckung mindestens `DEPTH_STEPS[i]`.
    pub covered: [u64; 4],
    /// Summe der Abdeckung je Fach von `BIN` Stellen.
    pub bins: Vec<u64>,
    /// Je Fach (nach dem Anfang des Reads): Hauptzeilen mit Ort, davon Duplikate, und Hauptzeilen
    /// ohne Duplikate mit MAPQ 0.
    pub bin_reads: Vec<u32>,
    pub bin_dups: Vec<u32>,
    pub bin_mapq0: Vec<u32>,
    /// Verteilung der Abdeckung je Stelle; nur für Chromosomen 1–22, X, Y und M, sonst leer.
    pub depth_hist: Vec<u64>,
    /// Summe der MAPQ der Hauptzeilen mit Ort ohne Duplikate (für den Mittelwert).
    pub mapq_sum: u64,
    pub mapq_reads: u64,
    /// Chromosom 1 bis 22 (mit oder ohne `chr`).
    pub autosome: bool,
}

/// Alle Zahlen einer BAM.
#[derive(Debug, Clone)]
pub struct BamQc {
    pub contigs: Vec<ContigQc>,
    pub flags: FlagStats,
    pub records: u64,
    /// Hauptzeilen ohne Ort.
    pub unmapped: u64,
    /// MAPQ der Hauptzeilen mit Ort, ohne Duplikate.
    pub mapq: Vec<u64>,
    /// Insertgröße (TLEN > 0, beide mit Ort auf demselben Contig, Hauptzeile, ohne Duplikate).
    pub insert: Vec<u64>,
    /// Read-Länge der Hauptzeilen.
    pub read_length: Vec<u64>,
    /// Basenqualität der Hauptzeilen: Summe, Zahl der Basen, Basen mit Qualität ≥ 30.
    pub quality_sum: u64,
    pub quality_bases: u64,
    pub q30_bases: u64,
    /// Sätze mit diesem Tag, Index = die zwei Buchstaben als u16.
    pub tags: Vec<u64>,
    /// Ab hier: Sätze, die zur Abdeckung zählen (mit Ort, keine Nebenzeile, kein Duplikat,
    /// Prüfung bestanden; ergänzende Alignments zählen mit).
    pub aligned_bases: u64,
    pub mismatches: u64,
    pub inserted_bases: u64,
    pub deleted_bases: u64,
    pub soft_clipped_bases: u64,
    /// Sätze mit NM-Tag (nur diese gehen in die Fehlerrate ein) und deren ausgerichtete Basen.
    pub nm_records: u64,
    pub nm_aligned_bases: u64,
    /// Wie `samtools stats`: Summe der NM-Werte und Basen aus M, I, = und X, über alle
    /// Hauptzeilen mit Ort (Duplikate eingeschlossen, ohne ergänzende Alignments).
    pub nm_sum: u64,
    pub nm_query_bases: u64,
    /// Verteilung der Abdeckung je Stelle auf den Autosomen.
    pub depth_autosomes: Vec<u64>,
    /// Ausgerichtete Stücke, die länger als der Ring waren und gekürzt wurden.
    pub clipped_blocks: u64,
    // Zustand der Abdeckung.
    ring: Vec<i32>,
    tid: i64,
    next: u64,
    depth: i64,
    pending_end: u64,
}

/// Fehler: Die Datensätze kommen nicht nach Ort sortiert.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NotSorted {
    pub tid: i32,
    pub pos: i64,
}

impl std::fmt::Display for NotSorted {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "BAM ist nicht nach Ort sortiert (Satz auf Contig {} an {})",
            self.tid, self.pos
        )
    }
}

impl std::error::Error for NotSorted {}

/// Chromosom 1–22, X, Y oder M (mit oder ohne `chr`).
pub fn is_main(name: &str) -> bool {
    let n = name.strip_prefix("chr").unwrap_or(name);
    is_autosome(name) || matches!(n, "X" | "Y" | "M" | "MT")
}

fn is_autosome(name: &str) -> bool {
    let n = name.strip_prefix("chr").unwrap_or(name);
    matches!(n.parse::<u32>(), Ok(1..=22)) && !n.starts_with('0')
}

impl BamQc {
    /// `contigs`: Namen und Längen aus dem Kopf, in der Reihenfolge der Contig-Nummern.
    pub fn new(contigs: Vec<(String, u64)>) -> Self {
        BamQc {
            contigs: contigs
                .into_iter()
                .map(|(name, length)| {
                    let n_bins = length.div_ceil(BIN) as usize;
                    ContigQc {
                        autosome: is_autosome(&name),
                        depth_hist: if is_main(&name) {
                            vec![0; MAX_DEPTH + 1]
                        } else {
                            Vec::new()
                        },
                        bins: vec![0; n_bins],
                        bin_reads: vec![0; n_bins],
                        bin_dups: vec![0; n_bins],
                        bin_mapq0: vec![0; n_bins],
                        name,
                        length,
                        ..ContigQc::default()
                    }
                })
                .collect(),
            flags: FlagStats::default(),
            records: 0,
            unmapped: 0,
            mapq: vec![0; 256],
            insert: vec![0; MAX_INSERT + 1],
            read_length: vec![0; MAX_READ_LENGTH + 1],
            quality_sum: 0,
            quality_bases: 0,
            q30_bases: 0,
            tags: vec![0; 1 << 16],
            aligned_bases: 0,
            mismatches: 0,
            inserted_bases: 0,
            deleted_bases: 0,
            soft_clipped_bases: 0,
            nm_records: 0,
            nm_aligned_bases: 0,
            nm_sum: 0,
            nm_query_bases: 0,
            depth_autosomes: vec![0; MAX_DEPTH + 1],
            clipped_blocks: 0,
            ring: vec![0; RING as usize],
            tid: -1,
            next: 0,
            depth: 0,
            pending_end: 0,
        }
    }

    /// Zählt einen rohen Datensatz (ohne Längenfeld), in der Reihenfolge der Datei.
    pub fn push(&mut self, rec: &[u8]) -> Result<(), NotSorted> {
        if rec.len() < 32 {
            return Ok(());
        }
        self.records += 1;
        self.flags.add(rec);
        let i32_at =
            |at: usize| i32::from_le_bytes([rec[at], rec[at + 1], rec[at + 2], rec[at + 3]]);
        let tid = i32_at(0);
        let pos = i32_at(4);
        let l_name = rec[8] as usize;
        let mapq = rec[9];
        let n_cigar = u16::from_le_bytes([rec[12], rec[13]]) as usize;
        let flag = u16::from_le_bytes([rec[14], rec[15]]);
        let l_seq = u32::from_le_bytes([rec[16], rec[17], rec[18], rec[19]]) as usize;
        let mate_tid = i32_at(20);
        let tlen = i32_at(28);
        let cigar_at = 32 + l_name;
        let seq_at = cigar_at + 4 * n_cigar;
        let qual_at = seq_at + l_seq.div_ceil(2);
        let tags_at = qual_at + l_seq;
        if tags_at > rec.len() {
            return Ok(());
        }
        let mapped = flag & UNMAPPED == 0 && tid >= 0;
        let primary = flag & (SECONDARY | SUPPLEMENTARY) == 0;
        let duplicate = flag & DUPLICATE != 0;
        let failed = flag & QC_FAIL != 0;

        // Tags zählen, NM merken.
        let mut nm: Option<i64> = None;
        let tags = &rec[tags_at..];
        let mut at = 0usize;
        while at + 3 <= tags.len() {
            let key = u16::from_le_bytes([tags[at], tags[at + 1]]);
            let kind = tags[at + 2];
            self.tags[key as usize] += 1;
            let value = at + 3;
            let size = match kind {
                b'A' | b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                b'Z' | b'H' => match tags[value..].iter().position(|&b| b == 0) {
                    Some(end) => end + 1,
                    None => break,
                },
                b'B' => {
                    if value + 5 > tags.len() {
                        break;
                    }
                    let each = match tags[value] {
                        b'c' | b'C' => 1,
                        b's' | b'S' => 2,
                        _ => 4,
                    };
                    let n = u32::from_le_bytes([
                        tags[value + 1],
                        tags[value + 2],
                        tags[value + 3],
                        tags[value + 4],
                    ]) as usize;
                    5 + each * n
                }
                _ => break,
            };
            if value + size > tags.len() {
                break;
            }
            if &tags[at..at + 2] == b"NM" {
                let v = &tags[value..value + size];
                nm = Some(match kind {
                    b'c' => v[0] as i8 as i64,
                    b'C' => v[0] as i64,
                    b's' => i16::from_le_bytes([v[0], v[1]]) as i64,
                    b'S' => u16::from_le_bytes([v[0], v[1]]) as i64,
                    b'i' => i32::from_le_bytes([v[0], v[1], v[2], v[3]]) as i64,
                    b'I' => u32::from_le_bytes([v[0], v[1], v[2], v[3]]) as i64,
                    _ => 0,
                });
            }
            at = value + size;
        }

        if primary {
            self.read_length[l_seq.min(MAX_READ_LENGTH)] += 1;
            let qual = &rec[qual_at..tags_at];
            if qual.first().is_some_and(|&q| q != 0xff) {
                let (mut sum, mut q30) = (0u64, 0u64);
                for &q in qual {
                    sum += q as u64;
                    q30 += (q >= 30) as u64;
                }
                self.quality_sum += sum;
                self.q30_bases += q30;
                self.quality_bases += qual.len() as u64;
            }
            if !mapped {
                self.unmapped += 1;
            }
        }

        if !mapped {
            return Ok(());
        }
        let pos = pos.max(0) as u64;
        // Fertige Stellen bis zum Anfang dieses Satzes auszählen.
        if (tid as i64) != self.tid {
            if (tid as i64) < self.tid {
                return Err(NotSorted {
                    tid,
                    pos: pos as i64,
                });
            }
            self.finish_contig();
            self.skip_contigs(tid as i64);
            self.tid = tid as i64;
            self.next = 0;
            self.depth = 0;
            self.pending_end = 0;
        }
        if pos < self.next {
            return Err(NotSorted {
                tid,
                pos: pos as i64,
            });
        }
        let contig = tid as usize;
        if contig >= self.contigs.len() {
            return Ok(());
        }
        self.advance(pos);

        let c = &mut self.contigs[contig];
        if primary {
            c.reads += 1;
            c.duplicates += duplicate as u64;
            let bin = ((pos / BIN) as usize).min(c.bin_reads.len().saturating_sub(1));
            if let Some(n) = c.bin_reads.get_mut(bin) {
                *n += 1;
                c.bin_dups[bin] += duplicate as u32;
                c.bin_mapq0[bin] += (!duplicate && mapq == 0) as u32;
            }
            if !duplicate {
                self.mapq[mapq as usize] += 1;
                c.mapq_sum += mapq as u64;
                c.mapq_reads += 1;
                if flag & PAIRED != 0
                    && flag & MATE_UNMAPPED == 0
                    && mate_tid == tid
                    && tlen > 0
                    && !failed
                {
                    self.insert[(tlen as usize).min(MAX_INSERT)] += 1;
                }
            }
        }
        // Fehlerrate wie `samtools stats`: Hauptzeilen mit Ort, Duplikate eingeschlossen, NM je
        // Basen aus M, I, = und X.
        if primary && let Some(nm) = nm {
            let mut query = 0u64;
            for k in 0..n_cigar {
                let o = cigar_at + 4 * k;
                let op = u32::from_le_bytes([rec[o], rec[o + 1], rec[o + 2], rec[o + 3]]);
                if matches!(op & 0xf, 0 | 1 | 7 | 8) {
                    query += (op >> 4) as u64;
                }
            }
            self.nm_sum += nm.max(0) as u64;
            self.nm_query_bases += query;
        }
        if duplicate || failed || flag & SECONDARY != 0 {
            return Ok(());
        }
        // Ausgerichtete Stücke eintragen.
        let length = c.length;
        let mut ref_at = pos;
        let (mut aligned, mut ins, mut del) = (0u64, 0u64, 0u64);
        for k in 0..n_cigar {
            let o = cigar_at + 4 * k;
            let op = u32::from_le_bytes([rec[o], rec[o + 1], rec[o + 2], rec[o + 3]]);
            let len = (op >> 4) as u64;
            match op & 0xf {
                0 | 7 | 8 => {
                    let start = ref_at.min(length);
                    let mut end = (ref_at + len).min(length);
                    if end > start {
                        if end - self.next >= RING {
                            end = self.next + RING - 1;
                            self.clipped_blocks += 1;
                        }
                        if end > start {
                            self.ring[(start & MASK) as usize] += 1;
                            // Am Ende des Contigs braucht es kein −1: Dort wird nicht mehr gezählt.
                            if end < length {
                                self.ring[(end & MASK) as usize] -= 1;
                                self.pending_end = self.pending_end.max(end + 1);
                            } else {
                                self.pending_end = self.pending_end.max(end);
                            }
                        }
                    }
                    aligned += len;
                    ref_at += len;
                }
                1 => ins += len,
                2 => {
                    del += len;
                    ref_at += len;
                }
                3 => ref_at += len,
                4 => self.soft_clipped_bases += len,
                _ => {}
            }
        }
        self.contigs[contig].bases += aligned;
        self.aligned_bases += aligned;
        self.inserted_bases += ins;
        self.deleted_bases += del;
        if let Some(nm) = nm {
            self.mismatches += (nm - ins as i64 - del as i64).max(0) as u64;
            self.nm_records += 1;
            self.nm_aligned_bases += aligned;
        }
        Ok(())
    }

    /// Zählt die Stellen `[self.next, to)` des aktuellen Contigs aus.
    fn advance(&mut self, to: u64) {
        if self.tid < 0 {
            return;
        }
        let c = &mut self.contigs[self.tid as usize];
        let to = to.min(c.length);
        let mut p = self.next;
        while p < to {
            if self.depth == 0 && p >= self.pending_end {
                // Nichts mehr offen: der Rest bis `to` hat Abdeckung 0.
                if c.autosome {
                    self.depth_autosomes[0] += to - p;
                }
                if let Some(h) = c.depth_hist.first_mut() {
                    *h += to - p;
                }
                p = to;
                break;
            }
            let slot = (p & MASK) as usize;
            self.depth += self.ring[slot] as i64;
            self.ring[slot] = 0;
            let d = self.depth.max(0) as u64;
            if d > 0 {
                for (i, step) in DEPTH_STEPS.iter().enumerate() {
                    c.covered[i] += (d >= *step) as u64;
                }
                c.bins[(p / BIN) as usize] += d;
            }
            if c.autosome {
                self.depth_autosomes[(d as usize).min(MAX_DEPTH)] += 1;
            }
            if !c.depth_hist.is_empty() {
                c.depth_hist[(d as usize).min(MAX_DEPTH)] += 1;
            }
            p += 1;
        }
        self.next = p.max(self.next);
    }

    fn finish_contig(&mut self) {
        if self.tid >= 0 && (self.tid as usize) < self.contigs.len() {
            let length = self.contigs[self.tid as usize].length;
            // Danach ist der Ring leer: Alle Einträge liegen vor dem Ende des Contigs.
            self.advance(length);
        }
    }

    /// Contigs zwischen dem aktuellen und `to` haben keine Sätze: Abdeckung überall 0.
    fn skip_contigs(&mut self, to: i64) {
        let from = (self.tid + 1).max(0);
        for t in from..to.min(self.contigs.len() as i64) {
            let c = &mut self.contigs[t as usize];
            if c.autosome {
                self.depth_autosomes[0] += c.length;
            }
            if let Some(h) = c.depth_hist.first_mut() {
                *h += c.length;
            }
        }
    }

    /// Am Ende der Datei: offene Stellen auszählen und Contigs ohne Sätze nachtragen.
    pub fn finish(&mut self) {
        self.finish_contig();
        self.skip_contigs(self.contigs.len() as i64);
        self.tid = self.contigs.len() as i64;
        self.ring = Vec::new();
    }

    /// Name eines Tags aus seinem Index in `tags`.
    pub fn tag_name(index: usize) -> [u8; 2] {
        (index as u16).to_le_bytes()
    }

    /// Kennzahlen einer Verteilung (Fach = Wert): Mittel, Median, Streuung.
    pub fn summary(hist: &[u64]) -> HistSummary {
        let n: u64 = hist.iter().sum();
        if n == 0 {
            return HistSummary::default();
        }
        let sum: f64 = hist
            .iter()
            .enumerate()
            .map(|(v, c)| v as f64 * *c as f64)
            .sum();
        let mean = sum / n as f64;
        let var: f64 = hist
            .iter()
            .enumerate()
            .map(|(v, c)| (v as f64 - mean).powi(2) * *c as f64)
            .sum::<f64>()
            / n as f64;
        let quantile = |q: f64| {
            let goal = (q * n as f64).ceil().max(1.0) as u64;
            let mut seen = 0;
            for (v, c) in hist.iter().enumerate() {
                seen += c;
                if seen >= goal {
                    return v as u64;
                }
            }
            (hist.len() - 1) as u64
        };
        HistSummary {
            n,
            mean,
            sd: var.sqrt(),
            median: quantile(0.5),
            p05: quantile(0.05),
            p95: quantile(0.95),
        }
    }
}

/// Kennzahlen einer Verteilung.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct HistSummary {
    pub n: u64,
    pub mean: f64,
    pub sd: f64,
    pub median: u64,
    pub p05: u64,
    pub p95: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Baut einen rohen Datensatz mit Name "r", Folge aus A und Qualität 30.
    #[allow(clippy::too_many_arguments)]
    fn record(
        tid: i32,
        pos: i32,
        mapq: u8,
        flag: u16,
        cigar: &[(u32, u32)],
        mate_tid: i32,
        tlen: i32,
        tags: &[u8],
    ) -> Vec<u8> {
        let l_seq: u32 = cigar
            .iter()
            .filter(|(op, _)| matches!(op, 0 | 1 | 4 | 7 | 8))
            .map(|(_, n)| n)
            .sum();
        let mut r = Vec::new();
        r.extend_from_slice(&tid.to_le_bytes());
        r.extend_from_slice(&pos.to_le_bytes());
        r.push(2); // l_read_name
        r.push(mapq);
        r.extend_from_slice(&0u16.to_le_bytes());
        r.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
        r.extend_from_slice(&flag.to_le_bytes());
        r.extend_from_slice(&l_seq.to_le_bytes());
        r.extend_from_slice(&mate_tid.to_le_bytes());
        r.extend_from_slice(&0i32.to_le_bytes());
        r.extend_from_slice(&tlen.to_le_bytes());
        r.extend_from_slice(b"r\0");
        for (op, n) in cigar {
            r.extend_from_slice(&((n << 4) | op).to_le_bytes());
        }
        r.extend(std::iter::repeat_n(0x11u8, (l_seq as usize).div_ceil(2)));
        r.extend(std::iter::repeat_n(30u8, l_seq as usize));
        r.extend_from_slice(tags);
        r
    }

    fn nm(v: u8) -> Vec<u8> {
        vec![b'N', b'M', b'C', v]
    }

    #[test]
    fn abdeckung_je_stelle_stimmt() {
        // Contig 1 (Länge 100), Contig 2 (Länge 50, ohne Sätze).
        let mut qc = BamQc::new(vec![("chr1".into(), 100), ("chr2".into(), 50)]);
        // Zwei Reads: [10,20) und [15,25) mit einer Deletion [18,20) im zweiten nicht.
        let mut tags = nm(1);
        tags.extend_from_slice(b"RGZx\0");
        qc.push(&record(
            0,
            10,
            60,
            0x1 | 0x2 | 0x40,
            &[(0, 10)],
            0,
            15,
            &tags,
        ))
        .unwrap();
        qc.push(&record(
            0,
            15,
            60,
            0x1 | 0x2 | 0x80,
            &[(0, 3), (2, 2), (0, 5)],
            0,
            -15,
            &nm(3),
        ))
        .unwrap();
        // Ein Duplikat zählt nicht zur Abdeckung.
        qc.push(&record(0, 15, 60, 0x1 | 0x400, &[(0, 10)], 0, 0, &[]))
            .unwrap();
        qc.finish();
        let c = &qc.contigs[0];
        // Stellen mit Abdeckung: [10,15) 1, [15,18) 2, [18,20) 1, [20,25) 1.
        assert_eq!(c.covered[0], 15);
        assert_eq!(c.bases, 18);
        assert_eq!(c.bins[0], 5 + 6 + 2 + 5);
        assert_eq!(qc.depth_autosomes[0], 100 - 15 + 50);
        assert_eq!(qc.depth_autosomes[1], 12);
        assert_eq!(qc.depth_autosomes[2], 3);
        assert_eq!(c.reads, 3);
        assert_eq!(c.duplicates, 1);
        // NM 1 ohne Indel = 1 Fehler, NM 3 mit 2 deletierten Basen = 1 Fehler.
        assert_eq!(qc.mismatches, 2);
        assert_eq!(qc.deleted_bases, 2);
        assert_eq!(qc.nm_sum, 4);
        assert_eq!(c.bin_reads[0], 3);
        assert_eq!(c.bin_dups[0], 1);
        assert_eq!(c.bin_mapq0[0], 0);
        assert_eq!(c.depth_hist[1], 12);
        assert_eq!(c.depth_hist[0], 85);
        assert_eq!(qc.contigs[1].depth_hist[0], 50);
        assert_eq!(qc.insert[15], 1);
        assert_eq!(qc.tags[u16::from_le_bytes(*b"NM") as usize], 2);
        assert_eq!(qc.tags[u16::from_le_bytes(*b"RG") as usize], 1);
        assert_eq!(qc.q30_bases, qc.quality_bases);
        assert_eq!(qc.flags.passed.duplicates, 1);
    }

    #[test]
    fn unsortiert_wird_erkannt() {
        let mut qc = BamQc::new(vec![("chr1".into(), 100)]);
        qc.push(&record(0, 50, 60, 0, &[(0, 10)], -1, 0, &[]))
            .unwrap();
        assert!(
            qc.push(&record(0, 10, 60, 0, &[(0, 10)], -1, 0, &[]))
                .is_err()
        );
    }

    #[test]
    fn kennzahlen_einer_verteilung() {
        let mut h = vec![0u64; 11];
        h[2] = 1;
        h[4] = 2;
        h[10] = 1;
        let s = BamQc::summary(&h);
        assert_eq!(s.n, 4);
        assert!((s.mean - 5.0).abs() < 1e-9);
        assert_eq!(s.median, 4);
    }

    #[test]
    fn autosomen_werden_erkannt() {
        assert!(is_autosome("chr1"));
        assert!(is_autosome("22"));
        assert!(!is_autosome("chrX"));
        assert!(!is_autosome("chr23"));
        assert!(!is_autosome("chr1_KI270706v1_random"));
        assert!(!is_autosome("chr01"));
    }
}
