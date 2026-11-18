//! Referenzgenom als 2-Bit-Folge mit Contig-Tabelle und Liste der N-Strecken.

/// Basen-Codes: A=0, C=1, G=2, T=3.
pub const BASES: &[u8; 4] = b"ACGT";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub len: u64,
    /// Anfang im verketteten Vorwärtsstrang.
    pub offset: u64,
}

/// Eine Strecke aus N (oder anderen Mehrdeutigkeits-Codes) im verketteten Vorwärtsstrang.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hole {
    pub offset: u64,
    pub len: u64,
}

/// Code einer Base, `None` für N und alle Mehrdeutigkeits-Codes.
#[inline]
pub fn encode(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// 2-Bit-gepackte Folge, 32 Basen je Wort, Base `i` in den Bits `2*(i%32)`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Packed {
    pub len: u64,
    pub words: Vec<u64>,
}

impl Packed {
    pub fn push(&mut self, code: u8) {
        let i = self.len;
        if i.is_multiple_of(32) {
            self.words.push(0);
        }
        *self.words.last_mut().expect("Wort vorhanden") |= u64::from(code & 3) << (2 * (i % 32));
        self.len += 1;
    }

    #[inline]
    pub fn get(&self, i: u64) -> u8 {
        (self.words[(i / 32) as usize] >> (2 * (i % 32)) & 3) as u8
    }
}

/// Das Referenzgenom, wie es in den Index eingeht.
#[derive(Debug, Clone, Default)]
pub struct Reference {
    pub contigs: Vec<Contig>,
    pub holes: Vec<Hole>,
    /// Vorwärtsstrang aller Contigs hintereinander. N-Strecken sind mit festen Pseudozufallsbasen
    /// gefüllt (wie bei BWA), damit sie keine künstlichen Treffer erzeugen; `holes` merkt sie sich.
    pub forward: Packed,
    filler: u64,
}

impl PartialEq for Reference {
    /// Gleich heißt: gleiche Contigs, N-Strecken und Basen. Der Zustand des Füllgenerators zählt nicht.
    fn eq(&self, other: &Self) -> bool {
        self.contigs == other.contigs && self.holes == other.holes && self.forward == other.forward
    }
}

impl Eq for Reference {}

impl Reference {
    pub fn new() -> Self {
        Reference {
            filler: 0x2545_F491_4F6C_DD1D,
            ..Default::default()
        }
    }

    pub fn len(&self) -> u64 {
        self.forward.len
    }

    pub fn is_empty(&self) -> bool {
        self.forward.len == 0
    }

    fn fill_base(&mut self) -> u8 {
        // xorshift64: gleiche Füllung auf jedem System und bei jedem Lauf.
        self.filler ^= self.filler << 13;
        self.filler ^= self.filler >> 7;
        self.filler ^= self.filler << 17;
        (self.filler >> 33 & 3) as u8
    }

    /// Hängt ein Contig an.
    pub fn push_contig(&mut self, name: &str, seq: &[u8]) {
        let offset = self.forward.len;
        for (i, &b) in seq.iter().enumerate() {
            match encode(b) {
                Some(code) => self.forward.push(code),
                None => {
                    let pos = offset + i as u64;
                    match self.holes.last_mut() {
                        Some(h) if h.offset + h.len == pos && h.offset >= offset => h.len += 1,
                        _ => self.holes.push(Hole {
                            offset: pos,
                            len: 1,
                        }),
                    }
                    let code = self.fill_base();
                    self.forward.push(code);
                }
            }
        }
        self.contigs.push(Contig {
            name: name.to_string(),
            len: seq.len() as u64,
            offset,
        });
    }

    /// Contig, in dem die Vorwärtsposition liegt.
    pub fn contig_at(&self, pos: u64) -> Option<usize> {
        if pos >= self.len() {
            return None;
        }
        let i = self.contigs.partition_point(|c| c.offset <= pos);
        // Contigs der Länge 0 teilen sich den Anfang mit dem Nachfolger; das letzte passende gilt.
        Some(i - 1)
    }

    /// Liegt eine N-Strecke im Bereich `[pos, pos + len)`?
    pub fn overlaps_hole(&self, pos: u64, len: u64) -> bool {
        let i = self.holes.partition_point(|h| h.offset + h.len <= pos);
        self.holes.get(i).is_some_and(|h| h.offset < pos + len)
    }

    /// Text für das Suffix-Array: Vorwärtsstrang, dann sein reverses Komplement, dann Endzeichen.
    /// Zeichen sind Code + 1, das Endzeichen ist 0.
    pub fn index_text(&self) -> Vec<u8> {
        let l = self.len();
        let mut text = Vec::with_capacity(2 * l as usize + 1);
        for i in 0..l {
            text.push(self.forward.get(i) + 1);
        }
        for i in (0..l).rev() {
            text.push(3 - self.forward.get(i) + 1);
        }
        text.push(0);
        text
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contigs_holes_and_text() {
        let mut r = Reference::new();
        r.push_contig("a", b"ACGTNNAC");
        r.push_contig("b", b"NGg");
        assert_eq!(r.len(), 11);
        assert_eq!(
            r.holes,
            [Hole { offset: 4, len: 2 }, Hole { offset: 8, len: 1 }]
        );
        assert_eq!(r.contig_at(7), Some(0));
        assert_eq!(r.contig_at(8), Some(1));
        assert_eq!(r.contig_at(11), None);
        assert!(r.overlaps_hole(2, 3));
        assert!(!r.overlaps_hole(0, 4));
        assert!(!r.overlaps_hole(6, 2));
        assert!(r.overlaps_hole(6, 3));
        let t = r.index_text();
        assert_eq!(t.len(), 23);
        assert_eq!(&t[..4], [1, 2, 3, 4]);
        assert_eq!(t[22], 0);
        // reverses Komplement: letztes Vorwärtszeichen g wird zu c am Anfang der zweiten Hälfte
        assert_eq!(t[11], 2);
        assert_eq!(t[12], 2);
        // gleiche Füllung bei Wiederholung
        let mut r2 = Reference::new();
        r2.push_contig("a", b"ACGTNNAC");
        r2.push_contig("b", b"NGg");
        assert_eq!(r, r2);
    }
}
