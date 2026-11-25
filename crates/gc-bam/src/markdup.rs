//! Duplikate markieren in einer nach Ort sortierten BAM.
//!
//! Das Verfahren ist das von `samtools markdup` (Vorgabe-Einstellungen, Modus "template"),
//! Schritt für Schritt nachgebaut, damit das Ergebnis Satz für Satz verglichen werden kann:
//!
//! - Ein Paar gilt als Duplikat eines anderen, wenn beide Reads an denselben Enden liegen. Als
//!   Ende zählt der Anfang des Moleküls einschließlich abgeschnittener Basen (bei einem Read auf
//!   dem Gegenstrang also sein rechtes Ende), dazu Richtung und die Rolle im Paar. Der Ort des
//!   Partners kommt aus `MC` (CIGAR des Partners), seine Qualität aus `ms`.
//! - Von mehreren gleichen Paaren bleibt das mit der höchsten Summe der Basenqualitäten ab 15
//!   (beide Reads); bei Gleichstand das mit dem im Alphabet früheren Namen.
//! - Ein Read ohne Partner (oder mit Partner ohne Ort) ist Duplikat, wenn an seinem Ende schon
//!   ein Read eines Paares liegt, sonst entscheidet die Qualität.
//! - Zusätzliche Zeilen (Flags 256 und 2048), Reads ohne Ort und Reads mit Flag 512 werden nicht
//!   bewertet und nicht markiert.
//!
//! Die Sätze laufen durch ein Fenster: Ein Satz wird geschrieben, sobald die Eingabe weiter als
//! `max_length` hinter seinem Ende ist. Der Speicherbedarf hängt deshalb nur von der Tiefe der
//! Abdeckung ab, nicht von der Größe der Datei.

use std::collections::{HashMap, VecDeque};
use std::io;

const FLAG_PAIRED: u16 = 0x1;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_MATE_REVERSE: u16 = 0x20;
const FLAG_READ1: u16 = 0x40;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QC_FAIL: u16 = 0x200;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Basenqualitäten unter diesem Wert zählen nicht (samtools `MD_MIN_QUALITY`).
const MIN_QUALITY: u8 = 15;

// Richtungen und Rollen mit den Zahlen von samtools.
const O_FF: u8 = 2;
const O_RR: u8 = 3;
const O_FR: u8 = 5;
const O_RF: u8 = 7;
const R_LE: u8 = 11;
const R_RI: u8 = 13;

#[derive(Debug, Clone, Copy)]
pub struct DupOptions {
    /// Längster erwarteter Read einschließlich abgeschnittener Basen (samtools `-l`, 300).
    /// Bestimmt, wie lange ein Satz im Fenster bleibt.
    pub max_length: i64,
    /// Optische Duplikate zählen; `None` = nicht zählen (wie samtools ohne `-d`). Das Flag der
    /// Sätze hängt nicht davon ab, nur die Zahlen und die Schätzung der Bibliotheksgröße.
    pub optical: Option<OpticalRule>,
}

impl Default for DupOptions {
    fn default() -> Self {
        DupOptions {
            max_length: 300,
            optical: Some(OpticalRule::default()),
        }
    }
}

/// Wann ein Duplikat als optisch gilt (Nachbar auf derselben Fläche statt Kopie aus der PCR).
/// Immer nötig: dieselbe Flowcell, Lane und Kachel bzw. dasselbe Bildfeld.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpticalRule {
    /// Illumina: höchster Abstand in x und in y, in Pixeln (Picard-Vorgabe 100).
    pub pixel: u32,
    /// MGI/DNBSEQ: höchster Abstand der laufenden Nummer im Bildfeld. Die Namen tragen keine
    /// Pixel; `C` und `R` sind die Nummern des Bildfelds.
    pub mgi_serial: u64,
}

impl Default for OpticalRule {
    fn default() -> Self {
        OpticalRule {
            pixel: 100,
            mgi_serial: 1000,
        }
    }
}

impl OpticalRule {
    /// Ob zwei Reads nach ihrem Namen Nachbarn auf derselben Fläche sind. Unbekannte
    /// Namensschemata sind nie optisch.
    pub fn neighbours(&self, a: &[u8], b: &[u8]) -> bool {
        let (Some(a), Some(b)) = (gc_io::readname::parse(a), gc_io::readname::parse(b)) else {
            return false;
        };
        if a.platform != b.platform
            || a.flowcell != b.flowcell
            || a.lane != b.lane
            || a.tile != b.tile
        {
            return false;
        }
        match (a.xy, b.xy, a.serial, b.serial) {
            (Some((ax, ay)), Some((bx, by)), _, _) => {
                ax.abs_diff(bx) <= self.pixel && ay.abs_diff(by) <= self.pixel
            }
            (_, _, Some(sa), Some(sb)) => sa.abs_diff(sb) <= self.mgi_serial,
            _ => false,
        }
    }
}

/// Zahlen wie in der Statistik von `samtools markdup -f`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DupStats {
    pub read: u64,
    pub written: u64,
    /// Nicht bewertet: ohne Ort, zusätzliche Zeilen, Flag 512.
    pub excluded: u64,
    pub examined: u64,
    /// Bewertete Reads mit Partner (Reads, nicht Paare).
    pub paired: u64,
    pub single: u64,
    /// Als Duplikat markierte Reads mit Partner.
    pub duplicate_pair: u64,
    pub duplicate_single: u64,
    /// Davon optisch: Reads mit Partner, deren Paar neben dem behaltenen Paar auf derselben
    /// Fläche lag (nur mit [`OpticalRule`]).
    pub optical_pair: u64,
    /// Größte Zahl von Sätzen, die gleichzeitig im Fenster lagen.
    pub window_peak: u64,
}

impl DupStats {
    /// Geschätzte Zahl verschiedener Moleküle in der Bibliothek (Formel aus Picard
    /// `DuplicationMetrics`, wie samtools sie übernimmt); `None`, wenn sie sich nicht schätzen lässt.
    /// Optische Duplikate zählen wie bei Picard nicht als Paare der Bibliothek; ohne
    /// [`OpticalRule`] ist ihre Zahl 0 und das Ergebnis das von samtools.
    pub fn estimated_library_size(&self) -> Option<u64> {
        let pairs = (self.paired - self.optical_pair.min(self.paired)) / 2;
        let unique = (self.paired - self.duplicate_pair) / 2;
        let duplicates = self.duplicate_pair / 2;
        if pairs == 0 || duplicates == 0 || unique == 0 || pairs <= unique {
            return None;
        }
        let (c, n) = (unique as f64, pairs as f64);
        let f = |x: f64| c / x - 1.0 + (-n / x).exp();
        let (mut lo, mut hi) = (1.0f64, 100.0f64);
        if f(lo * c) < 0.0 {
            return None;
        }
        while f(hi * c) > 0.0 {
            hi *= 10.0;
        }
        for _ in 0..40 {
            let r = (lo + hi) / 2.0;
            let u = f(r * c);
            if u > 0.0 {
                lo = r;
            } else if u < 0.0 {
                hi = r;
            } else {
                break;
            }
        }
        Some((c * (lo + hi) / 2.0) as u64)
    }
}

/// Schlüssel eines Reads: eigenes Ende, bei Paaren auch das des Partners.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Key {
    this_ref: i32,
    this_coord: i64,
    other_ref: i32,
    other_coord: i64,
    leftmost: u8,
    orientation: u8,
}

struct Entry {
    /// Der Satz mit Längenfeld.
    record: Vec<u8>,
    /// Ende des Reads (bewertete Reads) oder Position, für das Fenster.
    pos: i64,
    tid: i32,
    has_mate: bool,
    /// Schlüssel, die dieser Satz als Erster eingetragen hat; sie verfallen mit ihm.
    pair_key: Option<Key>,
    single_key: Option<Key>,
}

pub struct DupMarker {
    opt: DupOptions,
    queue: VecDeque<Entry>,
    /// Laufende Nummer des vordersten Satzes im Fenster.
    front: u64,
    /// Schlüssel zur laufenden Nummer des bisher besten Satzes.
    pairs: HashMap<Key, u64>,
    singles: HashMap<Key, u64>,
    /// Mitglieder von Gruppen gleicher Paare; optisch eingestuft wird erst, wenn die Gruppe das
    /// Fenster verlässt und das behaltene Paar feststeht.
    groups: HashMap<Key, Vec<u64>>,
    prev_tid: i32,
    prev_pos: i64,
    stats: DupStats,
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// Felder eines rohen Satzes (ohne Längenfeld).
pub(crate) struct View<'a> {
    b: &'a [u8],
}

impl<'a> View<'a> {
    pub(crate) fn new(b: &'a [u8]) -> io::Result<Self> {
        let v = View { b };
        if b.len() < 32 || v.aux_start() > b.len() {
            return Err(invalid("BAM-Satz ist kürzer als seine Felder"));
        }
        Ok(v)
    }
    fn i32_at(&self, at: usize) -> i32 {
        i32::from_le_bytes([self.b[at], self.b[at + 1], self.b[at + 2], self.b[at + 3]])
    }
    pub(crate) fn tid(&self) -> i32 {
        self.i32_at(0)
    }
    pub(crate) fn pos(&self) -> i64 {
        i64::from(self.i32_at(4))
    }
    fn name_len(&self) -> usize {
        self.b[8] as usize
    }
    fn n_cigar(&self) -> usize {
        u16::from_le_bytes([self.b[12], self.b[13]]) as usize
    }
    pub(crate) fn flag(&self) -> u16 {
        u16::from_le_bytes([self.b[14], self.b[15]])
    }
    fn l_seq(&self) -> usize {
        self.i32_at(16).max(0) as usize
    }
    fn mtid(&self) -> i32 {
        self.i32_at(20)
    }
    fn mpos(&self) -> i64 {
        i64::from(self.i32_at(24))
    }
    /// Name ohne das abschließende Nullbyte.
    fn name(&self) -> &'a [u8] {
        let end = 32 + self.name_len();
        let name = &self.b[32..end.min(self.b.len())];
        name.strip_suffix(&[0]).unwrap_or(name)
    }
    fn cigar_start(&self) -> usize {
        32 + self.name_len()
    }
    fn qual_start(&self) -> usize {
        self.cigar_start() + 4 * self.n_cigar() + self.l_seq().div_ceil(2)
    }
    fn aux_start(&self) -> usize {
        self.qual_start() + self.l_seq()
    }
    /// CIGAR als (Buchstabe, Länge).
    fn cigar(&self) -> impl DoubleEndedIterator<Item = (u8, i64)> + 'a {
        let start = self.cigar_start();
        self.b[start..start + 4 * self.n_cigar()]
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| {
                let v = u32::from_le_bytes(*c);
                (
                    *b"MIDNSHP=X".get((v & 0xF) as usize).unwrap_or(&b'?'),
                    i64::from(v >> 4),
                )
            })
    }
    /// Anfang einschließlich abgeschnittener Basen, ab 1 gezählt (samtools `unclipped_start`).
    fn unclipped_start(&self) -> i64 {
        let clipped: i64 = self
            .cigar()
            .take_while(|&(op, _)| op == b'S' || op == b'H')
            .map(|(_, len)| len)
            .sum();
        self.pos() - clipped + 1
    }
    /// Position hinter der letzten Base auf der Referenz (htslib `bam_endpos`): Ein Satz ohne
    /// Ort oder ohne CIGAR belegt eine Base.
    pub(crate) fn end_pos(&self) -> i64 {
        let reference: i64 = if self.flag() & FLAG_UNMAPPED != 0 {
            0
        } else {
            self.cigar()
                .filter(|&(op, _)| matches!(op, b'M' | b'D' | b'N' | b'=' | b'X'))
                .map(|(_, len)| len)
                .sum()
        };
        self.pos() + reference.max(1)
    }
    /// Ende einschließlich abgeschnittener Basen (samtools `unclipped_end`): die Position hinter
    /// der letzten Base, ab 0 gezählt. Ein Satz ohne CIGAR belegt wie in htslib eine Base.
    fn unclipped_end(&self) -> i64 {
        let reference: i64 = self
            .cigar()
            .filter(|&(op, _)| matches!(op, b'M' | b'D' | b'N' | b'=' | b'X'))
            .map(|(_, len)| len)
            .sum();
        let clipped: i64 = self
            .cigar()
            .rev()
            .take_while(|&(op, _)| op == b'S' || op == b'H')
            .map(|(_, len)| len)
            .sum();
        self.pos() + reference.max(1) + clipped
    }
    /// Summe der Basenqualitäten ab 15.
    fn score(&self) -> i64 {
        let start = self.qual_start();
        self.b[start..start + self.l_seq()]
            .iter()
            .filter(|&&q| q >= MIN_QUALITY)
            .map(|&q| i64::from(q))
            .sum()
    }
    /// Sucht eine Zusatzangabe und gibt Typ und Wert roh zurück.
    fn aux(&self, tag: [u8; 2]) -> io::Result<Option<(u8, &'a [u8])>> {
        let cut = || invalid("Zusatzangaben eines BAM-Satzes sind abgeschnitten");
        let b = self.b;
        let mut at = self.aux_start();
        while at < b.len() {
            if at + 3 > b.len() {
                return Err(cut());
            }
            let (name, kind) = ([b[at], b[at + 1]], b[at + 2]);
            at += 3;
            let fixed = |kind: u8| match kind {
                b'A' | b'c' | b'C' => Some(1),
                b's' | b'S' => Some(2),
                b'i' | b'I' | b'f' => Some(4),
                _ => None,
            };
            let len = match kind {
                b'Z' | b'H' => b[at..].iter().position(|&c| c == 0).ok_or_else(cut)? + 1,
                b'B' => {
                    if at + 5 > b.len() {
                        return Err(cut());
                    }
                    let size = fixed(b[at]).ok_or_else(cut)?;
                    let n = u32::from_le_bytes([b[at + 1], b[at + 2], b[at + 3], b[at + 4]]);
                    5 + size * n as usize
                }
                other => fixed(other).ok_or_else(cut)?,
            };
            if at + len > b.len() {
                return Err(cut());
            }
            if name == tag {
                return Ok(Some((kind, &b[at..at + len])));
            }
            at += len;
        }
        Ok(None)
    }
    /// Ganzzahl einer Zusatzangabe.
    fn aux_int(&self, tag: [u8; 2]) -> io::Result<Option<i64>> {
        Ok(match self.aux(tag)? {
            Some((b'c', v)) => Some(i64::from(v[0] as i8)),
            Some((b'C', v)) => Some(i64::from(v[0])),
            Some((b's', v)) => Some(i64::from(i16::from_le_bytes([v[0], v[1]]))),
            Some((b'S', v)) => Some(i64::from(u16::from_le_bytes([v[0], v[1]]))),
            Some((b'i', v)) => Some(i64::from(i32::from_le_bytes([v[0], v[1], v[2], v[3]]))),
            Some((b'I', v)) => Some(i64::from(u32::from_le_bytes([v[0], v[1], v[2], v[3]]))),
            _ => None,
        })
    }
    fn has_mate(&self) -> bool {
        let f = self.flag();
        f & FLAG_PAIRED != 0
            && f & FLAG_MATE_UNMAPPED == 0
            && !(self.mtid() == -1 && self.mpos() == -1)
    }
}

/// Zerlegt einen CIGAR-Text in (Buchstabe, Länge); eine fehlende Zahl zählt als 1, `*` beendet
/// (wie samtools `unclipped_other_start`).
fn text_cigar(text: &[u8]) -> Vec<(u8, i64)> {
    let mut out = Vec::new();
    let mut at = 0;
    while at < text.len() && text[at] != b'*' {
        let mut num = 0i64;
        let digits = text[at..].iter().take_while(|c| c.is_ascii_digit()).count();
        if digits == 0 {
            num = 1;
        } else {
            for &c in &text[at..at + digits] {
                num = num.saturating_mul(10).saturating_add(i64::from(c - b'0'));
            }
            at += digits;
        }
        if at >= text.len() {
            break;
        }
        out.push((text[at], num));
        at += 1;
    }
    out
}

/// Anfang und Ende des Partners einschließlich abgeschnittener Basen, aus seiner Position und
/// seinem CIGAR (samtools `unclipped_other_start`, `unclipped_other_end`).
fn mate_ends(mpos: i64, cigar: &[(u8, i64)]) -> (i64, i64) {
    let leading: i64 = cigar
        .iter()
        .take_while(|&&(op, _)| op == b'S' || op == b'H')
        .map(|&(_, len)| len)
        .sum();
    let mut reference = 0;
    let mut skip = true;
    for &(op, len) in cigar {
        match op {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                reference += len;
                skip = false;
            }
            b'S' | b'H' if !skip => reference += len,
            _ => {}
        }
    }
    (mpos - leading + 1, mpos + reference)
}

/// Schlüssel eines Reads ohne Partner: Ende und Richtung.
fn single_key(v: &View<'_>) -> Key {
    let reverse = v.flag() & FLAG_REVERSE != 0;
    Key {
        this_ref: v.tid() + 1,
        this_coord: if reverse {
            v.unclipped_end()
        } else {
            v.unclipped_start()
        },
        other_ref: 0,
        other_coord: 0,
        leftmost: 0,
        orientation: if reverse { O_RR } else { O_FF },
    }
}

/// Schlüssel eines Reads mit Partner (samtools `make_pair_key`, Modus "template").
fn pair_key(v: &View<'_>) -> io::Result<Key> {
    let missing = || {
        invalid(format!(
            "Read {} hat keine Angabe MC: Die Eingabe braucht die Partner-Angaben (genome index sam --partner-angaben oder samtools fixmate -m)",
            String::from_utf8_lossy(v.name())
        ))
    };
    let mc = match v.aux(*b"MC")? {
        Some((b'Z', text)) => text_cigar(&text[..text.len() - 1]),
        _ => return Err(missing()),
    };
    let flag = v.flag();
    let reverse = flag & FLAG_REVERSE != 0;
    let mate_reverse = flag & FLAG_MATE_REVERSE != 0;
    let first = flag & FLAG_READ1 != 0;
    let (this_ref, other_ref) = (v.tid() + 1, v.mtid() + 1);
    let (mut this_coord, this_end) = (v.unclipped_start(), v.unclipped_end());
    let (mut other_coord, other_end) = mate_ends(v.mpos(), &mc);

    let leftmost = if this_ref != other_ref {
        this_ref < other_ref
    } else if reverse == mate_reverse {
        if reverse {
            this_end <= other_end
        } else {
            this_coord <= other_coord
        }
    } else if reverse {
        this_end <= other_coord
    } else {
        this_coord <= other_end
    };

    let orientation;
    if leftmost {
        if reverse == mate_reverse {
            other_coord = other_end;
            orientation = if reverse != first { O_FF } else { O_RR };
        } else if reverse {
            orientation = O_RF;
            this_coord = this_end;
        } else {
            orientation = O_FR;
            other_coord = other_end;
        }
    } else if reverse == mate_reverse {
        this_coord = this_end;
        orientation = if reverse != first { O_RR } else { O_FF };
    } else if reverse {
        orientation = O_FR;
        this_coord = this_end;
    } else {
        orientation = O_RF;
        other_coord = other_end;
    }
    Ok(Key {
        this_ref,
        this_coord,
        other_ref,
        other_coord,
        leftmost: if leftmost { R_LE } else { R_RI },
        orientation,
    })
}

impl DupMarker {
    pub fn new(opt: DupOptions) -> Self {
        DupMarker {
            opt,
            queue: VecDeque::new(),
            front: 0,
            pairs: HashMap::new(),
            singles: HashMap::new(),
            groups: HashMap::new(),
            prev_tid: 0,
            prev_pos: 0,
            stats: DupStats::default(),
        }
    }

    pub fn stats(&self) -> DupStats {
        self.stats
    }

    fn entry(&mut self, nr: u64) -> &mut Entry {
        &mut self.queue[(nr - self.front) as usize]
    }

    fn mark(&mut self, nr: u64) {
        let record = &mut self.entry(nr).record;
        let flag = u16::from_le_bytes([record[18], record[19]]) | FLAG_DUPLICATE;
        record[18..20].copy_from_slice(&flag.to_le_bytes());
    }

    /// Summe der Qualitäten eines Satzes im Fenster, mit der seines Partners aus `ms`.
    fn pair_score(&self, nr: u64) -> io::Result<i64> {
        let v = View {
            b: &self.queue[(nr - self.front) as usize].record[4..],
        };
        let mate = v.aux_int(*b"ms")?.ok_or_else(|| {
            invalid(format!(
                "Read {} hat keine Angabe ms: Die Eingabe braucht die Partner-Angaben (genome index sam --partner-angaben oder samtools fixmate -m)",
                String::from_utf8_lossy(v.name())
            ))
        })?;
        Ok(v.score() + mate)
    }

    /// Nimmt den nächsten Satz (mit Längenfeld) und gibt alle Sätze an `sink`, die das Fenster
    /// verlassen. Die Reihenfolge der Sätze bleibt die der Eingabe.
    pub fn push(
        &mut self,
        record: &[u8],
        sink: &mut impl FnMut(&[u8]) -> io::Result<()>,
    ) -> io::Result<()> {
        if record.len() < 4
            || u32::from_le_bytes([record[0], record[1], record[2], record[3]]) as usize
                != record.len() - 4
        {
            return Err(invalid("BAM-Satz passt nicht zu seinem Längenfeld"));
        }
        let v = View::new(&record[4..])?;
        let (tid, pos, flag) = (v.tid(), v.pos(), v.flag());
        if tid >= 0 && (tid < self.prev_tid || (tid == self.prev_tid && pos < self.prev_pos)) {
            return Err(invalid(format!(
                "Eingabe ist nicht nach Ort sortiert (Read {})",
                String::from_utf8_lossy(v.name())
            )));
        }
        self.prev_tid = tid;
        self.prev_pos = pos;
        self.stats.read += 1;

        let nr = self.front + self.queue.len() as u64;
        let mut entry = Entry {
            record: Vec::new(),
            pos,
            tid,
            has_mate: v.has_mate(),
            pair_key: None,
            single_key: None,
        };
        let exclude = FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_UNMAPPED | FLAG_QC_FAIL;
        let keys = if flag & exclude != 0 {
            self.stats.excluded += 1;
            None
        } else {
            self.stats.examined += 1;
            let single = single_key(&v);
            entry.pos = single.this_coord;
            let pair = if entry.has_mate {
                Some(pair_key(&v)?)
            } else {
                None
            };
            Some((single, pair))
        };
        let name_is_smaller = |other: &[u8]| v.name() < View { b: &other[4..] }.name();
        entry.record = record.to_vec();
        // Erst ins Fenster, dann vergleichen: Markiert wird immer ein Satz im Fenster.
        let mut smaller_name_than: Option<bool> = None;
        if let Some((_, Some(pair))) = &keys
            && let Some(&old) = self.pairs.get(pair)
        {
            smaller_name_than = Some(name_is_smaller(
                &self.queue[(old - self.front) as usize].record,
            ));
        }
        self.queue.push_back(entry);
        self.stats.window_peak = self.stats.window_peak.max(self.queue.len() as u64);

        match keys {
            Some((single, Some(pair))) => {
                self.stats.paired += 1;
                // Gegen Reads ohne Partner: Der ohne Partner ist immer das Duplikat.
                match self.singles.get(&single).copied() {
                    None => {
                        self.singles.insert(single, nr);
                        self.entry(nr).single_key = Some(single);
                    }
                    Some(old) => {
                        if !self.entry(old).has_mate {
                            self.mark(old);
                            self.singles.insert(single, nr);
                            self.stats.duplicate_single += 1;
                        }
                    }
                }
                match self.pairs.get(&pair).copied() {
                    None => {
                        self.pairs.insert(pair, nr);
                        self.entry(nr).pair_key = Some(pair);
                    }
                    Some(old) => {
                        let (old_score, new_score) = (self.pair_score(old)?, self.pair_score(nr)?);
                        let tie = if new_score != old_score {
                            0
                        } else if smaller_name_than == Some(true) {
                            1
                        } else {
                            -1
                        };
                        if self.opt.optical.is_some() {
                            let members = self.groups.entry(pair).or_default();
                            if members.is_empty() {
                                members.push(old);
                            }
                            members.push(nr);
                        }
                        if new_score + tie > old_score {
                            self.mark(old);
                            self.pairs.insert(pair, nr);
                        } else {
                            self.mark(nr);
                        }
                        self.stats.duplicate_pair += 1;
                    }
                }
            }
            Some((single, None)) => {
                self.stats.single += 1;
                match self.singles.get(&single).copied() {
                    None => {
                        self.singles.insert(single, nr);
                        self.entry(nr).single_key = Some(single);
                    }
                    Some(old) => {
                        if self.entry(old).has_mate {
                            self.mark(nr);
                        } else {
                            let score = |m: &Self, n: u64| {
                                View {
                                    b: &m.queue[(n - m.front) as usize].record[4..],
                                }
                                .score()
                            };
                            if score(self, nr) > score(self, old) {
                                self.mark(old);
                                self.singles.insert(single, nr);
                            } else {
                                self.mark(nr);
                            }
                        }
                        self.stats.duplicate_single += 1;
                    }
                }
            }
            None => {}
        }

        // Alles schreiben, was nicht mehr getroffen werden kann.
        while let Some(first) = self.queue.front() {
            if first.pos + self.opt.max_length > pos && first.tid == tid && (tid != -1 || pos != -1)
            {
                break;
            }
            self.pop(sink)?;
        }
        Ok(())
    }

    /// Zählt die optischen Duplikate einer Gruppe gegen ihr endgültig behaltenes Paar. Alle
    /// Mitglieder liegen noch im Fenster: Der Satz, der die Gruppe begann, ist ihr frühestes.
    fn classify_optical(&mut self, key: &Key) {
        let (Some(rule), Some(members)) = (self.opt.optical, self.groups.remove(key)) else {
            return;
        };
        let Some(&kept) = self.pairs.get(key) else {
            return;
        };
        let name = |n: u64| View {
            b: &self.queue[(n - self.front) as usize].record[4..],
        }
        .name();
        let optical = members
            .iter()
            .filter(|&&m| m != kept && rule.neighbours(name(m), name(kept)))
            .count();
        self.stats.optical_pair += optical as u64;
    }

    fn pop(&mut self, sink: &mut impl FnMut(&[u8]) -> io::Result<()>) -> io::Result<()> {
        if let Some(key) = self.queue.front().and_then(|e| e.pair_key) {
            self.classify_optical(&key);
        }
        if let Some(entry) = self.queue.pop_front() {
            self.front += 1;
            if let Some(key) = entry.pair_key {
                self.pairs.remove(&key);
            }
            if let Some(key) = entry.single_key {
                self.singles.remove(&key);
            }
            sink(&entry.record)?;
            self.stats.written += 1;
        }
        Ok(())
    }

    /// Schreibt den Rest des Fensters.
    pub fn finish(
        mut self,
        sink: &mut impl FnMut(&[u8]) -> io::Result<()>,
    ) -> io::Result<DupStats> {
        while !self.queue.is_empty() {
            self.pop(sink)?;
        }
        Ok(self.stats)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Baut einen Satz mit Längenfeld. `cigar` als Text, Qualität für alle Basen gleich.
    #[allow(clippy::too_many_arguments)]
    fn rec(
        name: &str,
        flag: u16,
        tid: i32,
        pos: i32,
        cigar: &str,
        mpos: i32,
        qual: u8,
        tags: &[(&str, &str)],
    ) -> Vec<u8> {
        let ops = text_cigar(cigar.as_bytes());
        let l_seq: i64 = ops
            .iter()
            .filter(|(op, _)| matches!(op, b'M' | b'I' | b'S' | b'=' | b'X'))
            .map(|&(_, n)| n)
            .sum();
        let mut b = Vec::new();
        b.extend_from_slice(&tid.to_le_bytes());
        b.extend_from_slice(&pos.to_le_bytes());
        b.push(name.len() as u8 + 1);
        b.push(60);
        b.extend_from_slice(&0u16.to_le_bytes());
        b.extend_from_slice(&(ops.len() as u16).to_le_bytes());
        b.extend_from_slice(&flag.to_le_bytes());
        b.extend_from_slice(&(l_seq as u32).to_le_bytes());
        b.extend_from_slice(&(if flag & FLAG_PAIRED != 0 { tid } else { -1 }).to_le_bytes());
        b.extend_from_slice(&mpos.to_le_bytes());
        b.extend_from_slice(&0i32.to_le_bytes());
        b.extend_from_slice(name.as_bytes());
        b.push(0);
        for (op, n) in ops {
            let code = b"MIDNSHP=X".iter().position(|&c| c == op).unwrap() as u32;
            b.extend_from_slice(&((n as u32) << 4 | code).to_le_bytes());
        }
        b.extend(std::iter::repeat_n(0x11, (l_seq as usize).div_ceil(2)));
        b.extend(std::iter::repeat_n(qual, l_seq as usize));
        for (tag, value) in tags {
            b.extend_from_slice(tag.as_bytes());
            match value.parse::<i32>() {
                Ok(n) => {
                    b.push(b'i');
                    b.extend_from_slice(&n.to_le_bytes());
                }
                Err(_) => {
                    b.push(b'Z');
                    b.extend_from_slice(value.as_bytes());
                    b.push(0);
                }
            }
        }
        let mut out = (b.len() as u32).to_le_bytes().to_vec();
        out.extend_from_slice(&b);
        out
    }

    /// Ein Paar vorwärts/rückwärts: Read 1 bei `pos`, Read 2 bei `mpos`, je 50 Basen.
    fn pair(name: &str, pos: i32, cigar1: &str, mpos: i32, qual: u8) -> [Vec<u8>; 2] {
        let ms = (50 * i32::from(qual)).to_string();
        [
            rec(
                name,
                0x1 | 0x2 | 0x20 | 0x40,
                0,
                pos,
                cigar1,
                mpos,
                qual,
                &[("MC", "50M"), ("ms", &ms)],
            ),
            rec(
                name,
                0x1 | 0x2 | 0x10 | 0x80,
                0,
                mpos,
                "50M",
                pos,
                qual,
                &[("MC", cigar1), ("ms", &ms)],
            ),
        ]
    }

    /// Lässt die Sätze nach Ort sortiert durchlaufen und gibt je Satz (Name, Flag) zurück.
    fn run(mut records: Vec<Vec<u8>>) -> io::Result<(Vec<(String, u16)>, DupStats)> {
        records.sort_by_key(|r| crate::sort_key(&r[4..]));
        let mut out = Vec::new();
        let mut sink = |r: &[u8]| {
            let v = View { b: &r[4..] };
            out.push((String::from_utf8_lossy(v.name()).into_owned(), v.flag()));
            Ok(())
        };
        let mut marker = DupMarker::new(DupOptions::default());
        for r in &records {
            marker.push(r, &mut sink)?;
        }
        let stats = marker.finish(&mut sink)?;
        Ok((out, stats))
    }

    fn marked(out: &[(String, u16)]) -> Vec<&str> {
        let mut names: Vec<&str> = out
            .iter()
            .filter(|(_, f)| f & FLAG_DUPLICATE != 0)
            .map(|(n, _)| n.as_str())
            .collect();
        names.sort_unstable();
        names
    }

    #[test]
    fn the_pair_with_the_better_quality_stays() {
        let mut records = Vec::new();
        records.extend(pair("a", 1000, "50M", 1300, 30));
        records.extend(pair("b", 1000, "50M", 1300, 40));
        records.extend(pair("c", 1000, "50M", 1300, 20));
        // anderes Ende des Partners: kein Duplikat
        records.extend(pair("d", 1000, "50M", 1301, 10));
        let (out, stats) = run(records).unwrap();
        assert_eq!(marked(&out), ["a", "a", "c", "c"]);
        assert_eq!(out.len(), 8);
        assert_eq!(
            (stats.paired, stats.duplicate_pair, stats.single),
            (8, 4, 0)
        );
        assert_eq!(stats.written, 8);
    }

    #[test]
    fn a_tie_is_decided_by_the_name() {
        let mut records = Vec::new();
        records.extend(pair("y", 500, "50M", 800, 30));
        records.extend(pair("x", 500, "50M", 800, 30));
        records.extend(pair("z", 500, "50M", 800, 30));
        let (out, _) = run(records).unwrap();
        assert_eq!(marked(&out), ["y", "y", "z", "z"]);
    }

    #[test]
    fn clipped_bases_count_for_the_end() {
        // "b" beginnt wegen 5 abgeschnittener Basen 5 weiter rechts, das Molekül beginnt am
        // selben Ort. "c" hat denselben Ort wie "b", aber nichts abgeschnitten.
        let mut records = Vec::new();
        records.extend(pair("a", 1000, "50M", 1300, 40));
        records.extend(pair("b", 1005, "5S45M", 1300, 30));
        records.extend(pair("c", 1005, "50M", 1300, 30));
        let (out, _) = run(records).unwrap();
        assert_eq!(marked(&out), ["b", "b"]);
    }

    #[test]
    fn a_read_without_mate_loses_against_a_pair() {
        let mut records = Vec::new();
        records.extend(pair("paar", 2000, "50M", 2300, 20));
        // ohne Partner, bessere Qualität, gleicher Anfang: trotzdem das Duplikat
        records.push(rec("allein", 0, 0, 2000, "50M", -1, 40, &[]));
        // zwei ohne Partner am selben Ort: der bessere bleibt
        records.push(rec("e1", 0, 0, 5000, "50M", -1, 20, &[]));
        records.push(rec("e2", 0, 0, 5000, "50M", -1, 40, &[]));
        // Gegenstrang zählt mit seinem rechten Ende und ist kein Duplikat von e1/e2
        records.push(rec("e3", 0x10, 0, 5000, "50M", -1, 10, &[]));
        let (out, stats) = run(records).unwrap();
        assert_eq!(marked(&out), ["allein", "e1"]);
        assert_eq!((stats.single, stats.duplicate_single), (4, 2));
    }

    #[test]
    fn extra_lines_and_unplaced_reads_are_not_touched_and_order_is_kept() {
        let mut records = Vec::new();
        records.extend(pair("a", 1000, "50M", 1300, 40));
        records.extend(pair("b", 1000, "50M", 1300, 30));
        records.push(rec(
            "b",
            0x800 | 0x1 | 0x40,
            0,
            1000,
            "20M30H",
            1300,
            30,
            &[("MC", "50M"), ("ms", "1500")],
        ));
        records.push(rec("u", 0x4, -1, -1, "", -1, 30, &[]));
        let mut sorted = records.clone();
        sorted.sort_by_key(|r| crate::sort_key(&r[4..]));
        let (out, stats) = run(records).unwrap();
        assert_eq!(stats.excluded, 2);
        let names: Vec<String> = sorted
            .iter()
            .map(|r| String::from_utf8_lossy(View { b: &r[4..] }.name()).into_owned())
            .collect();
        assert_eq!(
            out.iter().map(|(n, _)| n.clone()).collect::<Vec<_>>(),
            names
        );
        let extra = out.iter().find(|(_, f)| f & 0x800 != 0).unwrap();
        assert_eq!(extra.1 & FLAG_DUPLICATE, 0);
        assert_eq!(marked(&out), ["b", "b"]);
    }

    #[test]
    fn missing_mate_tags_and_wrong_order_are_errors() {
        let no_mc = rec(
            "a",
            0x1 | 0x20 | 0x40,
            0,
            100,
            "50M",
            300,
            30,
            &[("ms", "100")],
        );
        let e = run(vec![no_mc]).unwrap_err().to_string();
        assert!(e.contains("MC"), "{e}");
        // ms fehlt: fällt erst auf, wenn zwei Paare verglichen werden
        let a = rec(
            "a",
            0x1 | 0x20 | 0x40,
            0,
            100,
            "50M",
            300,
            30,
            &[("MC", "50M")],
        );
        let b = rec(
            "b",
            0x1 | 0x20 | 0x40,
            0,
            100,
            "50M",
            300,
            30,
            &[("MC", "50M")],
        );
        let e = run(vec![a, b]).unwrap_err().to_string();
        assert!(e.contains("ms"), "{e}");

        let mut marker = DupMarker::new(DupOptions::default());
        let mut sink = |_: &[u8]| Ok(());
        marker
            .push(&rec("a", 0, 0, 500, "50M", -1, 30, &[]), &mut sink)
            .unwrap();
        let e = marker
            .push(&rec("b", 0, 0, 499, "50M", -1, 30, &[]), &mut sink)
            .unwrap_err()
            .to_string();
        assert!(e.contains("nicht nach Ort sortiert"), "{e}");
    }

    #[test]
    fn mate_ends_follow_the_cigar_text() {
        assert_eq!(text_cigar(b"5S45M"), [(b'S', 5), (b'M', 45)]);
        assert_eq!(text_cigar(b"*"), []);
        // Anfang ab 1 gezählt, abgeschnittene Basen abgezogen; Ende hinter der letzten Base
        assert_eq!(mate_ends(100, &text_cigar(b"5S40M2D3M2S")), (96, 147));
        assert_eq!(mate_ends(100, &text_cigar(b"50M")), (101, 150));
    }

    #[test]
    fn optical_neighbours_follow_the_name() {
        let rule = OpticalRule::default();
        // MGI: gleiches Bildfeld, laufende Nummer 1 daneben
        assert!(rule.neighbours(
            b"E250160857L1C001R00200000021",
            b"E250160857L1C001R00200000022"
        ));
        assert!(rule.neighbours(
            b"E250160857L1C001R00200000021",
            b"E250160857L1C001R00200001021"
        ));
        assert!(!rule.neighbours(
            b"E250160857L1C001R00200000021",
            b"E250160857L1C001R00200001022"
        ));
        // anderes Bildfeld, andere Lane oder Flowcell: nie optisch, auch bei naher Nummer
        assert!(!rule.neighbours(
            b"E250160857L1C001R00200000021",
            b"E250160857L1C001R00300000021"
        ));
        assert!(!rule.neighbours(
            b"E250160857L1C001R00200000021",
            b"E250160857L2C001R00200000021"
        ));
        assert!(!rule.neighbours(
            b"E250160857L1C001R00200000021",
            b"E250160051L1C001R00200000021"
        ));
        // Illumina: Pixelabstand in x und y
        assert!(rule.neighbours(b"A1:1:FC:2:1101:1000:2000", b"A1:1:FC:2:1101:1100:1900"));
        assert!(!rule.neighbours(b"A1:1:FC:2:1101:1000:2000", b"A1:1:FC:2:1101:1101:2000"));
        assert!(!rule.neighbours(b"A1:1:FC:2:1101:1000:2000", b"A1:1:FC:2:1102:1000:2000"));
        // unbekanntes Schema
        assert!(!rule.neighbours(b"read1", b"read2"));
    }

    #[test]
    fn optical_duplicates_are_counted_but_flags_stay_the_same() {
        let mut records = Vec::new();
        records.extend(pair("E250160857L1C001R00200000021", 1000, "50M", 1300, 30));
        records.extend(pair("E250160857L1C001R00200000022", 1000, "50M", 1300, 20));
        records.extend(pair("E250160857L1C002R00200000023", 1000, "50M", 1300, 10));
        let with = run(records.clone()).unwrap();
        let mut sorted = records;
        sorted.sort_by_key(|r| crate::sort_key(&r[4..]));
        let mut out = Vec::new();
        let mut sink = |r: &[u8]| {
            let v = View { b: &r[4..] };
            out.push((String::from_utf8_lossy(v.name()).into_owned(), v.flag()));
            Ok(())
        };
        let mut marker = DupMarker::new(DupOptions {
            optical: None,
            ..DupOptions::default()
        });
        for r in &sorted {
            marker.push(r, &mut sink).unwrap();
        }
        let without = marker.finish(&mut sink).unwrap();
        assert_eq!(with.0, out);
        // Nur das Paar aus demselben Bildfeld ist optisch, je Read gezählt
        assert_eq!((with.1.duplicate_pair, with.1.optical_pair), (4, 2));
        assert_eq!((without.duplicate_pair, without.optical_pair), (4, 0));
    }

    #[test]
    fn optical_duplicates_are_judged_against_the_finally_kept_pair() {
        // C ist am besten und liegt zwischen A und B; A und B sind je 75 Pixel von C entfernt,
        // aber 150 voneinander. Das Ergebnis darf nicht von der Reihenfolge abhängen.
        let a = pair("M:1:FC:1:1101:1000:500", 1000, "50M", 1300, 30);
        let b = pair("M:1:FC:1:1101:1150:500", 1000, "50M", 1300, 20);
        let c = pair("M:1:FC:1:1101:1075:500", 1000, "50M", 1300, 40);
        let count = |order: [&[Vec<u8>; 2]; 3]| {
            let mut firsts: Vec<&Vec<u8>> = order.iter().map(|p| &p[0]).collect();
            firsts.extend(order.iter().map(|p| &p[1]));
            let mut marker = DupMarker::new(DupOptions::default());
            let mut sink = |_: &[u8]| Ok(());
            for r in firsts {
                marker.push(r, &mut sink).unwrap();
            }
            marker.finish(&mut sink).unwrap().optical_pair
        };
        assert_eq!(count([&a, &b, &c]), 4);
        assert_eq!(count([&c, &a, &b]), 4);
        assert_eq!(count([&b, &c, &a]), 4);
    }

    #[test]
    fn optical_duplicates_do_not_count_for_the_library() {
        let base = DupStats {
            paired: 251_952,
            duplicate_pair: 52_830,
            ..DupStats::default()
        };
        let optical = DupStats {
            optical_pair: 20_000,
            ..base
        };
        let (a, b) = (
            base.estimated_library_size().unwrap(),
            optical.estimated_library_size().unwrap(),
        );
        assert!(b > a, "{b} > {a}");
        // Nur optische Duplikate: keine Schätzung möglich
        let all = DupStats {
            optical_pair: 52_830,
            ..base
        };
        assert_eq!(all.estimated_library_size(), None);
    }

    #[test]
    fn library_size_is_estimated_like_picard() {
        let stats = DupStats {
            paired: 251_952,
            duplicate_pair: 52_830,
            ..DupStats::default()
        };
        // Wert von samtools markdup 1.24 für dieselben Zahlen (bench/m5/markdup_vergleich.json)
        assert_eq!(stats.estimated_library_size(), Some(256_745));
        assert_eq!(DupStats::default().estimated_library_size(), None);
    }
}
