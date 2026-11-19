//! BAM-Index (`.bai`) für eine nach Ort sortierte BAM.
//!
//! Der Index sagt, an welchen Stellen der Datei die Datensätze eines Bereichs stehen
//! (SAM-Spezifikation, Abschnitt 5). Er hat je Contig zwei Teile:
//!
//! - Fächer: Das Contig ist in sechs Ebenen in immer kleinere Abschnitte geteilt (512 Mbp, 64 Mbp,
//!   8 Mbp, 1 Mbp, 128 kbp, 16 kbp). Jeder Datensatz gehört in das kleinste Fach, das ihn ganz
//!   enthält. Je Fach sind die Dateistücke vermerkt, in denen seine Datensätze stehen.
//! - Linearer Index: für jedes Fenster von 16.384 Basen die Stelle des ersten Datensatzes, der
//!   das Fenster berührt.
//!
//! Der Bau folgt htslib (`hts_idx_push`, `hts_idx_finish`), damit das Ergebnis mit dem von
//! `samtools index` verglichen werden kann: Kleine Fächer, deren Stücke nahe beieinander liegen,
//! gehen im übergeordneten Fach auf, benachbarte Stücke werden verbunden, und ein Sonderfach
//! (37450) trägt Anfang und Ende der Daten des Contigs und die Zahl der Reads mit und ohne Ort.

use std::collections::BTreeMap;
use std::io::{self, Write};

/// Kleinstes Fenster: 2^14 Basen.
const MIN_SHIFT: u32 = 14;
const LEVELS: u32 = 5;
/// Zahl der Fächer aller Ebenen; das Sonderfach hat diese Nummer plus 1.
const N_BINS: u32 = ((1 << (3 * LEVELS + 3)) - 1) / 7;
const META_BIN: u32 = N_BINS + 1;
/// Größte Position, die dieses Format fassen kann (2^29).
const MAX_POS: i64 = 1 << (MIN_SHIFT + 3 * LEVELS);
/// Fächer, deren Stücke weniger als so viele komprimierte Bytes überspannen, gehen im
/// übergeordneten Fach auf (htslib `HTS_MIN_MARKER_DIST`).
const MIN_MARKER_DIST: u64 = 0x10000;

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// Fach eines Bereichs `[beg, end)` (htslib `hts_reg2bin`).
pub fn reg2bin(beg: i64, end: i64) -> u32 {
    let end = end - 1;
    let mut shift = MIN_SHIFT;
    let mut first = ((1i64 << (3 * LEVELS)) - 1) / 7;
    let mut level = LEVELS;
    while level > 0 {
        if beg >> shift == end >> shift {
            return (first + (beg >> shift)) as u32;
        }
        level -= 1;
        shift += 3;
        first -= 1 << (3 * level);
    }
    0
}

fn first_bin_of_level(level: u32) -> u32 {
    ((1 << (3 * level)) - 1) / 7
}

fn parent(bin: u32) -> u32 {
    (bin - 1) >> 3
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Reference {
    /// Fach zu Dateistücken (Anfang, Ende) als virtuelle Positionen. Enthält auch das Sonderfach.
    pub bins: BTreeMap<u32, Vec<(u64, u64)>>,
    pub linear: Vec<u64>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Bai {
    pub references: Vec<Reference>,
    /// Datensätze ohne Contig am Ende der Datei.
    pub unplaced: u64,
}

/// Baut den Index, während die Datensätze in der Reihenfolge der Datei vorbeikommen.
pub struct BaiBuilder {
    refs: Vec<Reference>,
    seen: Vec<bool>,
    unplaced: u64,
    last_tid: i32,
    last_bin: u32,
    last_off: u64,
    last_coor: i64,
    save_tid: i32,
    save_bin: u32,
    save_off: u64,
    off_beg: u64,
    n_mapped: u64,
    n_unmapped: u64,
}

const NONE: u32 = u32::MAX;

impl BaiBuilder {
    /// `n_ref`: Zahl der Contigs im Kopf. `first_offset`: virtuelle Position des ersten
    /// Datensatzes (hinter dem Kopf).
    pub fn new(n_ref: usize, first_offset: u64) -> Self {
        BaiBuilder {
            refs: vec![Reference::default(); n_ref],
            seen: vec![false; n_ref],
            unplaced: 0,
            last_tid: -1,
            last_bin: NONE,
            last_off: first_offset,
            last_coor: -1,
            save_tid: -1,
            save_bin: NONE,
            save_off: first_offset,
            off_beg: first_offset,
            n_mapped: 0,
            n_unmapped: 0,
        }
    }

    fn add_chunk(&mut self, tid: i32, bin: u32, beg: u64, end: u64) {
        if tid >= 0 {
            self.refs[tid as usize]
                .bins
                .entry(bin)
                .or_default()
                .push((beg, end));
        }
    }

    /// Nächster Datensatz: Contig, Bereich `[beg, end)` ab 0 gezählt, virtuelle Position
    /// direkt hinter dem Datensatz, und ob er einen Ort hat (Flag 4 nicht gesetzt).
    pub fn push(
        &mut self,
        tid: i32,
        beg: i64,
        end: i64,
        end_offset: u64,
        mapped: bool,
    ) -> io::Result<()> {
        let (mut beg, mut end) = if tid < 0 { (-1, 0) } else { (beg, end) };
        if tid >= 0 && (beg > MAX_POS || end > MAX_POS) {
            return Err(invalid(format!(
                "Position {end} passt nicht in einen BAI-Index (höchstens {MAX_POS})"
            )));
        }
        if tid >= self.refs.len() as i32 {
            return Err(invalid(
                "Datensatz nennt ein Contig, das nicht im Kopf steht",
            ));
        }
        if self.last_tid != tid || (self.last_tid >= 0 && tid < 0) {
            if tid >= 0 && self.unplaced > 0 {
                return Err(invalid(
                    "nicht nach Ort sortiert: Datensätze ohne Contig stehen nicht am Ende",
                ));
            }
            if tid >= 0 && self.seen[tid as usize] {
                return Err(invalid(
                    "nicht nach Ort sortiert: ein Contig kommt in zwei getrennten Stücken",
                ));
            }
            self.last_tid = tid;
            self.last_bin = NONE;
        } else if tid >= 0 && self.last_coor > beg {
            return Err(invalid(format!(
                "nicht nach Ort sortiert: Position {} nach {}",
                beg + 1,
                self.last_coor + 1
            )));
        }
        if end < beg {
            return Err(invalid("Datensatz endet vor seinem Anfang"));
        }
        if tid >= 0 {
            self.seen[tid as usize] = true;
            beg = beg.max(0);
            if end <= 0 {
                end = 1;
            }
            // Linearer Index: jedes berührte Fenster merkt sich den ersten Datensatz.
            let linear = &mut self.refs[tid as usize].linear;
            let (from, to) = (
                (beg >> MIN_SHIFT) as usize,
                ((end - 1) >> MIN_SHIFT) as usize,
            );
            if linear.len() < to + 1 {
                linear.resize(to + 1, u64::MAX);
            }
            for slot in &mut linear[from..=to] {
                if *slot == u64::MAX {
                    *slot = self.last_off;
                }
            }
        } else {
            self.unplaced += 1;
        }
        let bin = reg2bin(beg, end);
        if self.last_bin != bin {
            if self.save_bin != NONE {
                self.add_chunk(self.save_tid, self.save_bin, self.save_off, self.last_off);
            }
            if self.last_bin == NONE && self.save_bin != NONE {
                // Wechsel des Contigs: Kennzahlen des vorigen ablegen.
                self.add_chunk(self.save_tid, META_BIN, self.off_beg, self.last_off);
                self.add_chunk(self.save_tid, META_BIN, self.n_mapped, self.n_unmapped);
                self.n_mapped = 0;
                self.n_unmapped = 0;
                self.off_beg = self.last_off;
            }
            self.save_off = self.last_off;
            self.save_bin = bin;
            self.last_bin = bin;
            self.save_tid = tid;
        }
        if mapped {
            self.n_mapped += 1;
        } else {
            self.n_unmapped += 1;
        }
        self.last_off = end_offset;
        self.last_coor = beg;
        Ok(())
    }

    /// Nächster Datensatz roh (ohne Längenfeld), wie er in der BAM steht.
    pub fn push_record(&mut self, record: &[u8], end_offset: u64) -> io::Result<()> {
        let v = crate::markdup::View::new(record)?;
        self.push(
            v.tid(),
            v.pos(),
            v.end_pos(),
            end_offset,
            v.flag() & 0x4 == 0,
        )
    }

    /// Schließt den Index ab. `final_offset`: virtuelle Position hinter dem letzten Datensatz.
    pub fn finish(mut self, final_offset: u64) -> Bai {
        if self.save_tid >= 0 {
            self.add_chunk(self.save_tid, self.save_bin, self.save_off, final_offset);
            self.add_chunk(self.save_tid, META_BIN, self.off_beg, final_offset);
            self.add_chunk(self.save_tid, META_BIN, self.n_mapped, self.n_unmapped);
        }
        for reference in &mut self.refs {
            // Leere Fenster zeigen auf den nächsten Datensatz dahinter.
            for i in (0..reference.linear.len().saturating_sub(1)).rev() {
                if reference.linear[i] == u64::MAX {
                    reference.linear[i] = reference.linear[i + 1];
                }
            }
            compress(&mut reference.bins);
        }
        Bai {
            references: self.refs,
            unplaced: self.unplaced,
        }
    }
}

/// Kleine Fächer gehen im übergeordneten auf, benachbarte Stücke werden verbunden
/// (htslib `compress_binning`).
fn compress(bins: &mut BTreeMap<u32, Vec<(u64, u64)>>) {
    for level in (1..=LEVELS).rev() {
        let start = first_bin_of_level(level);
        // Von hinten: tiefere Ebenen zuerst.
        let keys: Vec<u32> = bins.range(start..N_BINS).map(|(&k, _)| k).rev().collect();
        for key in keys {
            let list = bins.get_mut(&key).expect("Fach vorhanden");
            if level < LEVELS && list.len() > 1 {
                list.sort_unstable_by_key(|c| c.0);
            }
            let span = (list[list.len() - 1].1 >> 16).wrapping_sub(list[0].0 >> 16);
            if span < MIN_MARKER_DIST && bins.contains_key(&parent(key)) {
                let list = bins.remove(&key).expect("Fach vorhanden");
                bins.get_mut(&parent(key))
                    .expect("übergeordnetes Fach vorhanden")
                    .extend(list);
            }
        }
    }
    if let Some(list) = bins.get_mut(&0) {
        list.sort_unstable_by_key(|c| c.0);
    }
    for (_, list) in bins.range_mut(..N_BINS) {
        let mut m = 0;
        for l in 1..list.len() {
            if list[m].1 >> 16 >= list[l].0 >> 16 {
                if list[m].1 < list[l].1 {
                    list[m].1 = list[l].1;
                }
            } else {
                m += 1;
                list[m] = list[l];
            }
        }
        list.truncate(m + 1);
    }
}

impl Bai {
    /// Schreibt das Format `.bai` (SAM-Spezifikation 5.2), Fächer in aufsteigender Nummer.
    pub fn write(&self, out: &mut impl Write) -> io::Result<()> {
        let too_many = |_| invalid("zu viele Einträge für einen BAI-Index");
        out.write_all(b"BAI\x01")?;
        out.write_all(
            &i32::try_from(self.references.len())
                .map_err(too_many)?
                .to_le_bytes(),
        )?;
        for reference in &self.references {
            out.write_all(
                &i32::try_from(reference.bins.len())
                    .map_err(too_many)?
                    .to_le_bytes(),
            )?;
            for (bin, chunks) in &reference.bins {
                out.write_all(&bin.to_le_bytes())?;
                out.write_all(&i32::try_from(chunks.len()).map_err(too_many)?.to_le_bytes())?;
                for (beg, end) in chunks {
                    out.write_all(&beg.to_le_bytes())?;
                    out.write_all(&end.to_le_bytes())?;
                }
            }
            out.write_all(
                &i32::try_from(reference.linear.len())
                    .map_err(too_many)?
                    .to_le_bytes(),
            )?;
            for offset in &reference.linear {
                out.write_all(&offset.to_le_bytes())?;
            }
        }
        out.write_all(&self.unplaced.to_le_bytes())
    }

    /// Liest eine `.bai`-Datei (auch eine von samtools), zum Vergleichen und Prüfen.
    pub fn read(data: &[u8]) -> io::Result<Bai> {
        let cut = || invalid("BAI-Datei ist abgeschnitten");
        let mut rest = data;
        let mut take = |n: usize| -> io::Result<&[u8]> {
            if rest.len() < n {
                return Err(cut());
            }
            let (part, tail) = rest.split_at(n);
            rest = tail;
            Ok(part)
        };
        let int = |b: &[u8]| -> io::Result<usize> {
            usize::try_from(i32::from_le_bytes([b[0], b[1], b[2], b[3]])).map_err(|_| cut())
        };
        let long = |b: &[u8]| u64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
        if take(4)? != b"BAI\x01" {
            return Err(invalid("keine BAI-Datei"));
        }
        let n_ref = int(take(4)?)?;
        let mut references = Vec::new();
        for _ in 0..n_ref {
            let mut reference = Reference::default();
            for _ in 0..int(take(4)?)? {
                let b = take(4)?;
                let bin = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
                let mut chunks = Vec::new();
                for _ in 0..int(take(4)?)? {
                    let b = take(16)?;
                    chunks.push((long(&b[..8]), long(&b[8..])));
                }
                reference.bins.insert(bin, chunks);
            }
            for _ in 0..int(take(4)?)? {
                reference.linear.push(long(take(8)?));
            }
            references.push(reference);
        }
        let unplaced = take(8).map(long).unwrap_or(0);
        Ok(Bai {
            references,
            unplaced,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_follow_the_specification() {
        // Beispiele aus der SAM-Spezifikation (Abschnitt 5.3): kleinstes Fach ab 4681.
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(16383, 16384), 4681);
        assert_eq!(reg2bin(16384, 16385), 4682);
        assert_eq!(reg2bin(16383, 16385), 585);
        assert_eq!(reg2bin(0, 1 << 29), 0);
        assert_eq!(reg2bin(-1, 0), 4680);
        assert_eq!(parent(4681), 585);
        assert_eq!(first_bin_of_level(5), 4681);
    }

    #[test]
    fn index_round_trips_and_keeps_counts() {
        // Virtuelle Positionen: Block (Dateistelle) << 16 | Stelle im Block.
        let at = |block: u64, inner: u64| block << 16 | inner;
        let mut b = BaiBuilder::new(3, at(0, 100));
        b.push(0, 10, 60, at(0, 200), true).unwrap();
        b.push(0, 20, 70, at(0, 300), true).unwrap();
        b.push(0, 40_000, 40_050, at(0, 400), false).unwrap();
        // Contig 1 bleibt leer.
        b.push(2, 5, 55, at(70_000, 0), true).unwrap();
        b.push(-1, -1, 0, at(70_000, 100), false).unwrap();
        let bai = b.finish(at(70_000, 100));
        assert_eq!(bai.unplaced, 1);
        assert!(bai.references[1].bins.is_empty() && bai.references[1].linear.is_empty());
        let r0 = &bai.references[0];
        assert_eq!(r0.bins[&META_BIN], [(at(0, 100), at(0, 400)), (2, 1)]);
        // Fenster 0 und 2 belegt, Fenster 1 zeigt auf den nächsten Datensatz.
        assert_eq!(r0.linear, [at(0, 100), at(0, 300), at(0, 300)]);
        assert_eq!(bai.references[2].bins[&META_BIN][1], (1, 0));
        let mut bytes = Vec::new();
        bai.write(&mut bytes).unwrap();
        assert_eq!(Bai::read(&bytes).unwrap(), bai);
    }

    #[test]
    fn wrong_order_is_rejected() {
        let mut b = BaiBuilder::new(2, 0);
        b.push(0, 100, 150, 1, true).unwrap();
        assert!(b.push(0, 99, 150, 2, true).is_err());
        let mut b = BaiBuilder::new(2, 0);
        b.push(1, 100, 150, 1, true).unwrap();
        b.push(0, 100, 150, 2, true).unwrap();
        assert!(b.push(1, 200, 250, 3, true).is_err());
        let mut b = BaiBuilder::new(2, 0);
        b.push(-1, -1, 0, 1, false).unwrap();
        assert!(b.push(0, 1, 2, 2, true).is_err());
    }
}
