//! Zählung der Datensätze nach ihren Flags, wie `samtools flagstat`.
//!
//! Jede Zahl gibt es zweimal: für Datensätze, die die Qualitätsprüfung des Geräts bestanden
//! haben, und für die übrigen (Flag 512).

const PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ1: u16 = 0x40;
const READ2: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const QC_FAIL: u16 = 0x200;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// Zähler einer Gruppe (Qualitätsprüfung bestanden oder nicht). Die Namen der Felder folgen den
/// Zeilen von `samtools flagstat`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Counts {
    pub total: u64,
    pub primary: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub duplicates: u64,
    pub primary_duplicates: u64,
    pub mapped: u64,
    pub primary_mapped: u64,
    /// Ab hier nur Hauptzeilen von Reads aus Paaren.
    pub paired: u64,
    pub read1: u64,
    pub read2: u64,
    pub properly_paired: u64,
    pub both_mapped: u64,
    /// Read mit Ort, Partner ohne.
    pub singletons: u64,
    pub mate_on_other_contig: u64,
    /// Wie zuvor, MAPQ mindestens 5.
    pub mate_on_other_contig_mapq5: u64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FlagStats {
    pub passed: Counts,
    pub failed: Counts,
}

impl FlagStats {
    /// Zählt einen rohen Datensatz (ohne Längenfeld). Zu kurze Sätze werden nicht gezählt und
    /// als `false` gemeldet.
    pub fn add(&mut self, record: &[u8]) -> bool {
        if record.len() < 32 {
            return false;
        }
        let int = |at: usize| {
            i32::from_le_bytes([record[at], record[at + 1], record[at + 2], record[at + 3]])
        };
        let flag = u16::from_le_bytes([record[14], record[15]]);
        self.add_fields(flag, int(0), int(20), record[9]);
        true
    }

    /// Zählt einen Datensatz aus seinen Feldern (samtools `flagstat_loop`).
    pub fn add_fields(&mut self, flag: u16, tid: i32, mate_tid: i32, mapq: u8) {
        let c = if flag & QC_FAIL != 0 {
            &mut self.failed
        } else {
            &mut self.passed
        };
        let mapped = flag & UNMAPPED == 0;
        c.total += 1;
        if flag & SECONDARY != 0 {
            c.secondary += 1;
        } else if flag & SUPPLEMENTARY != 0 {
            c.supplementary += 1;
        } else {
            c.primary += 1;
            if flag & PAIRED != 0 {
                c.paired += 1;
                c.properly_paired += u64::from(flag & PROPER_PAIR != 0 && mapped);
                c.read1 += u64::from(flag & READ1 != 0);
                c.read2 += u64::from(flag & READ2 != 0);
                let mate_mapped = flag & MATE_UNMAPPED == 0;
                c.singletons += u64::from(mapped && !mate_mapped);
                if mapped && mate_mapped {
                    c.both_mapped += 1;
                    if mate_tid != tid {
                        c.mate_on_other_contig += 1;
                        c.mate_on_other_contig_mapq5 += u64::from(mapq >= 5);
                    }
                }
            }
            c.primary_mapped += u64::from(mapped);
            c.primary_duplicates += u64::from(flag & DUPLICATE != 0);
        }
        c.mapped += u64::from(mapped);
        c.duplicates += u64::from(flag & DUPLICATE != 0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_line_counts_what_its_name_says() {
        let mut s = FlagStats::default();
        // richtiges Paar, beide mit Ort
        s.add_fields(0x1 | 0x2 | 0x40, 0, 0, 60);
        s.add_fields(0x1 | 0x2 | 0x80 | 0x400, 0, 0, 60);
        // Paar über zwei Contigs, einer mit MAPQ 3
        s.add_fields(0x1 | 0x40, 0, 1, 60);
        s.add_fields(0x1 | 0x80, 1, 0, 3);
        // Read mit Ort, Partner ohne; der Partner selbst
        s.add_fields(0x1 | 0x8 | 0x40, 2, 2, 20);
        s.add_fields(0x1 | 0x4 | 0x80, 2, 2, 0);
        // zusätzliche Zeilen zählen nicht als Paar
        s.add_fields(0x800 | 0x1 | 0x40, 0, 0, 60);
        s.add_fields(0x100 | 0x400, 0, -1, 0);
        // ein einzelner Read ohne Ort, einer mit Flag 512
        s.add_fields(0x4, -1, -1, 0);
        s.add_fields(0x200, 0, -1, 30);
        let p = s.passed;
        assert_eq!(
            (p.total, p.primary, p.secondary, p.supplementary),
            (9, 7, 1, 1)
        );
        assert_eq!((p.duplicates, p.primary_duplicates), (2, 1));
        assert_eq!((p.mapped, p.primary_mapped), (7, 5));
        assert_eq!((p.paired, p.read1, p.read2), (6, 3, 3));
        assert_eq!((p.properly_paired, p.both_mapped, p.singletons), (2, 4, 1));
        assert_eq!(
            (p.mate_on_other_contig, p.mate_on_other_contig_mapq5),
            (2, 1)
        );
        assert_eq!((s.failed.total, s.failed.mapped), (1, 1));
        assert!(!s.add(&[0; 10]));
    }
}
