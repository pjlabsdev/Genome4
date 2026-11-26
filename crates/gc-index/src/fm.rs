//! FM-Index: BWT mit Zähltabelle in Blöcken und abgetastetem Suffix-Array.
//!
//! Drei Varianten desselben Formats:
//! - `Standard`: Zählstand alle 128 Basen, Suffix-Array jede 8. Zeile (Vorgabe).
//! - `Kompakt`: Zählstand alle 128 Basen, Suffix-Array jede 32. Zeile (wenig RAM).
//! - `Schnell`: Zählstand alle 64 Basen, Suffix-Array jede 8. Zeile (mehr RAM, kürzere Wege).

use std::io::{self, Read, Write};

use gc_log::Level;
use serde_json::json;

use crate::reference::{Contig, Hole, Packed, Reference};
use crate::sais;

const MAGIC: &[u8; 4] = b"GCFM";
const FORMAT: u32 = 2;
/// Bytes je Suffix-Array-Stichprobe (40-Bit-Positionen, reicht für 1 Billion Zeichen).
const SA_BYTES: usize = 5;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Variant {
    /// Vorgabe: Zählblöcke wie `Kompakt`, Stichproben so dicht wie `Schnell`.
    Standard,
    Kompakt,
    Schnell,
}

impl Variant {
    /// (Basen je Zählblock, Abstand der Suffix-Array-Stichproben)
    pub fn parameters(self) -> (u32, u32) {
        match self {
            Variant::Standard => (128, 8),
            Variant::Kompakt => (128, 32),
            Variant::Schnell => (64, 8),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Variant::Standard => "standard",
            Variant::Kompakt => "kompakt",
            Variant::Schnell => "schnell",
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum BuildError {
    #[error("Referenz ist leer")]
    Empty,
    #[error(
        "Referenz hat {bases} Basen; diese Index-Version verarbeitet höchstens {max} Basen \
         (32-Bit-Positionen über beide Stränge)"
    )]
    TooLong { bases: u64, max: u64 },
    #[error("Suffix-Array fehlerhaft in Zeile {row}")]
    SuffixArray { row: usize },
    #[error("Selbstprüfung fehlgeschlagen: {0}")]
    Check(String),
}

/// Ein Treffer im Referenzgenom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Hit {
    /// Index des Contigs.
    pub contig: u32,
    /// Position im Contig, 0-basiert, immer auf dem Vorwärtsstrang (linkes Ende).
    pub pos: u64,
    /// Das Muster passt auf den Gegenstrang.
    pub reverse: bool,
}

/// Ergebnis der Selbstprüfung beim Bau.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CheckReport {
    pub queries: u64,
    pub queries_with_hits: u64,
    pub occurrences_compared: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FmIndex {
    /// Zeilen der BWT: 2 × Referenzlänge + 1.
    rows: u64,
    /// Zeile, deren Suffix bei Textposition 0 beginnt (dort steht das Endzeichen in der BWT).
    primary: u64,
    /// Zeilen vor dem ersten Suffix, das mit Base c beginnt.
    c_before: [u64; 4],
    occ_interval: u32,
    sa_interval: u32,
    /// Beide Abstände sind Zweierpotenzen; geteilt wird durch Schieben um so viele Stellen.
    /// Eine echte Division je Zugriff kostete mehr als das Zählen selbst.
    occ_shift: u32,
    sa_shift: u32,
    /// Je Block: 4 Zählstände, dann `occ_interval / 32` Wörter BWT.
    blocks: Vec<u64>,
    /// Stichproben des Suffix-Arrays, je `SA_BYTES` Bytes little-endian.
    sa_samples: Vec<u8>,
    reference: Reference,
    /// Bereiche aller kurzen Muster, beim Anlegen des Index berechnet (nicht in der Datei).
    pub(crate) short: crate::smem::ShortTable,
    /// Zwischenspeicher für `locate_many` (leer = aus), siehe [`FmIndex::enable_locate_cache`].
    locate_cache: LocateCache,
    /// `true`: `occ4` zählt immer mit Schleifen und Verzweigungen wie früher, siehe
    /// [`FmIndex::set_occ_plain`].
    occ_plain: bool,
}

/// Plätze des Zwischenspeichers. Gehört nicht zum Inhalt des Index: Vergleiche sehen ihn nicht,
/// eine Kopie beginnt leer mit gleich vielen Plätzen.
#[derive(Debug, Default)]
struct LocateCache(Vec<std::sync::atomic::AtomicU64>);

impl Clone for LocateCache {
    fn clone(&self) -> Self {
        LocateCache(
            (0..self.0.len())
                .map(|_| std::sync::atomic::AtomicU64::new(0))
                .collect(),
        )
    }
}

impl PartialEq for LocateCache {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

impl Eq for LocateCache {}

/// Aufteilung eines Platzes im Zwischenspeicher von `locate_many`: Bit 63 belegt, darunter die
/// oberen Bits der Zeile (die unteren sind die Platznummer), darunter 33 Bit Position.
const CACHE_POS_BITS: u32 = 33;
const CACHE_POS_MASK: u64 = (1 << CACHE_POS_BITS) - 1;
const CACHE_TAG_MASK: u64 = (1 << (63 - CACHE_POS_BITS)) - 1;

const REP: [u64; 4] = [
    0,
    0x5555_5555_5555_5555,
    0xAAAA_AAAA_AAAA_AAAA,
    0xFFFF_FFFF_FFFF_FFFF,
];

/// Die unteren Bits aller 32 Basen eines Wortes.
const EVEN: u64 = 0x5555_5555_5555_5555;

/// Maske der unteren Bits der ersten `r` Basen eines Wortes (`r` von 1 bis 31).
#[inline]
fn low_bases(r: u32) -> u64 {
    EVEN & ((1u64 << (2 * r)) - 1)
}

/// Je Base des Wortes ein Bit (das untere ihrer zwei), gesetzt, wo die Base `c` steht.
#[inline]
fn matches_of(w: u64, c: u8) -> u64 {
    let x = w ^ REP[c as usize];
    !(x | x >> 1) & EVEN
}

/// Zählen über mehrere Wörter, ohne jedes Wort einzeln fertig zu zählen: `pairs` verdichtet ein
/// Wort nur bis zu Summen je Halbbyte, diese Summen werden addiert, und erst am Ende zählt
/// `fold` einmal ganz zusammen. Das braucht weniger Rechenschritte als eine volle Zählung je
/// Wort und keinen besonderen Befehl des Prozessors. Gilt für Wörter, in denen nur die Bits an
/// geraden Stellen gesetzt sein können (eines je Base): Je Halbbyte kommen dann höchstens 2
/// dazu, sieben Wörter passen in eine Summe. Ein Block hat deshalb höchstens
/// `MAX_OCC_INTERVAL` Basen.
#[inline]
fn pairs(bits: u64) -> u64 {
    const P: u64 = 0x3333_3333_3333_3333;
    (bits & P) + ((bits >> 2) & P)
}

#[inline]
fn fold(sums: u64) -> u64 {
    const B: u64 = 0x0F0F_0F0F_0F0F_0F0F;
    ((sums & B) + ((sums >> 4) & B)).wrapping_mul(0x0101_0101_0101_0101) >> 56
}

/// Gesetzte Bits zweier Wörter. Wo der Prozessor einen Befehl zum Zählen hat (Apple M und
/// alle 64-Bit-ARM, x86 mit `popcnt`), zählt dieser; sonst die Halbbyte-Summen wie oben.
#[inline]
fn count2(a: u64, b: u64) -> u64 {
    if cfg!(any(target_arch = "aarch64", target_feature = "popcnt")) {
        u64::from(a.count_ones() + b.count_ones())
    } else {
        fold(pairs(a) + pairs(b))
    }
}

/// Größter Abstand der Zählstände: sieben Wörter BWT je Block, als Zweierpotenz also vier.
pub(crate) const MAX_OCC_INTERVAL: u32 = 128;

/// Je Base C, G und T in den ersten `within` Basen der Wörter (A ist der Rest).
#[inline]
fn count_cgt(words: &[u64], within: u32) -> [u64; 3] {
    let (full, rest) = ((within / 32) as usize, within % 32);
    let (mut c, mut g, mut t) = (0u64, 0u64, 0u64);
    let mut add = |w: u64, keep: u64| {
        let (low, high) = (w & keep, (w >> 1) & keep);
        let both = low & high;
        c += pairs(low ^ both);
        g += pairs(high ^ both);
        t += pairs(both);
    };
    for &w in &words[..full] {
        add(w, EVEN);
    }
    if rest != 0 {
        add(words[full], low_bases(rest));
    }
    [fold(c), fold(g), fold(t)]
}

/// Je Base C, G und T in den Basen ab `within` bis zum Ende der Wörter (A ist der Rest).
#[inline]
fn count_cgt_tail(words: &[u64], within: u32) -> [u64; 3] {
    let (first, rest) = ((within / 32) as usize, within % 32);
    let (mut c, mut g, mut t) = (0u64, 0u64, 0u64);
    let mut add = |w: u64, keep: u64| {
        let (low, high) = (w & keep, (w >> 1) & keep);
        let both = low & high;
        c += pairs(low ^ both);
        g += pairs(high ^ both);
        t += pairs(both);
    };
    let mut from = first;
    if rest != 0 {
        add(words[first], EVEN & !((1u64 << (2 * rest)) - 1));
        from += 1;
    }
    for &w in &words[from..] {
        add(w, EVEN);
    }
    [fold(c), fold(g), fold(t)]
}

/// Base `c` in den ersten `within` Basen der Wörter.
#[inline]
fn count_one(words: &[u64], c: u8, within: u32) -> u64 {
    let (full, rest) = ((within / 32) as usize, within % 32);
    let mut sum = 0;
    for &w in &words[..full] {
        sum += pairs(matches_of(w, c));
    }
    if rest != 0 {
        sum += pairs(matches_of(words[full], c) & low_bases(rest));
    }
    fold(sum)
}

/// Zeichen vor Textposition 0: das Endzeichen.
pub(crate) const SENTINEL: u8 = 4;

/// Nimmt die BWT Zeile für Zeile entgegen und baut daraus die Tabellen des Index.
pub(crate) struct Builder {
    rows: u64,
    row: u64,
    primary: u64,
    occ_interval: u32,
    sa_interval: u32,
    stride: usize,
    blocks: Vec<u64>,
    counts: [u64; 4],
    sa_samples: Vec<u8>,
}

impl Builder {
    pub(crate) fn new(rows: u64, occ_interval: u32, sa_interval: u32) -> Self {
        assert!(
            (32..=MAX_OCC_INTERVAL).contains(&occ_interval)
                && occ_interval.is_power_of_two()
                && sa_interval.is_power_of_two(),
            "Abstände müssen Zweierpotenzen sein, der der Zählstände höchstens 128"
        );
        let stride = 4 + (occ_interval / 32) as usize;
        let n_blocks = (rows / u64::from(occ_interval)) as usize + 1;
        Builder {
            rows,
            row: 0,
            primary: 0,
            occ_interval,
            sa_interval,
            stride,
            blocks: vec![0u64; n_blocks * stride],
            counts: [0; 4],
            sa_samples: Vec::with_capacity(
                rows.div_ceil(u64::from(sa_interval)) as usize * SA_BYTES,
            ),
        }
    }

    /// Nächste Zeile: `prev` ist die Base vor dem Suffix (oder `SENTINEL`), `pos` seine Textposition.
    #[inline]
    pub(crate) fn push(&mut self, prev: u8, pos: u64) {
        let interval = u64::from(self.occ_interval);
        let within = (self.row % interval) as usize;
        let base = (self.row / interval) as usize * self.stride;
        if within == 0 {
            self.blocks[base..base + 4].copy_from_slice(&self.counts);
        }
        if prev == SENTINEL {
            self.primary = self.row; // als Code 0 gespeichert, beim Zählen abgezogen
        } else {
            self.counts[prev as usize] += 1;
            self.blocks[base + 4 + within / 32] |= u64::from(prev) << (2 * (within % 32));
        }
        if self.row.is_multiple_of(u64::from(self.sa_interval)) {
            self.sa_samples
                .extend_from_slice(&pos.to_le_bytes()[..SA_BYTES]);
        }
        self.row += 1;
    }

    pub(crate) fn finish(mut self, reference: Reference) -> FmIndex {
        assert_eq!(self.row, self.rows, "Zeilenzahl der BWT stimmt nicht");
        if self.rows.is_multiple_of(u64::from(self.occ_interval)) {
            let base = self.blocks.len() - self.stride;
            self.blocks[base..base + 4].copy_from_slice(&self.counts);
        }
        let mut c_before = [0u64; 4];
        let mut sum = 1; // das Endzeichen steht in Zeile 0
        for (before, count) in c_before.iter_mut().zip(self.counts) {
            *before = sum;
            sum += count;
        }
        let mut index = FmIndex {
            rows: self.rows,
            primary: self.primary,
            c_before,
            occ_interval: self.occ_interval,
            sa_interval: self.sa_interval,
            occ_shift: self.occ_interval.trailing_zeros(),
            sa_shift: self.sa_interval.trailing_zeros(),
            blocks: self.blocks,
            sa_samples: self.sa_samples,
            reference,
            short: crate::smem::ShortTable::default(),
            locate_cache: LocateCache::default(),
            occ_plain: false,
        };
        index.short = crate::smem::ShortTable::build(&index);
        index
    }
}

pub(crate) struct Rng(pub(crate) u64);

impl Rng {
    pub(crate) fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

/// Zeilenbereich aller Suffixe, die mit `pattern` beginnen, durch binäre Suche im vollen
/// Suffix-Array. Unabhängiger Vergleichsweg für die Selbstprüfung.
fn sa_interval_by_binary_search(text: &[u8], sa: &[u32], pattern: &[u8]) -> (usize, usize) {
    let cmp_prefix = |row: usize| -> std::cmp::Ordering {
        let start = sa[row] as usize;
        let end = (start + pattern.len()).min(text.len());
        let suffix = &text[start..end];
        // Zeichen im Text sind Code + 1.
        for (i, &p) in pattern.iter().enumerate() {
            match suffix.get(i) {
                None => return std::cmp::Ordering::Less,
                Some(&t) => match t.cmp(&(p + 1)) {
                    std::cmp::Ordering::Equal => {}
                    other => return other,
                },
            }
        }
        std::cmp::Ordering::Equal
    };
    let lower = {
        let (mut lo, mut hi) = (0, sa.len());
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if cmp_prefix(mid) == std::cmp::Ordering::Less {
                lo = mid + 1
            } else {
                hi = mid
            }
        }
        lo
    };
    let upper = {
        let (mut lo, mut hi) = (lower, sa.len());
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if cmp_prefix(mid) == std::cmp::Ordering::Greater {
                hi = mid
            } else {
                lo = mid + 1
            }
        }
        lo
    };
    (lower, upper)
}

impl FmIndex {
    /// Baut den Index. `check_queries > 0` prüft zusätzlich das Suffix-Array vollständig und
    /// vergleicht so viele Zufallsabfragen mit der binären Suche im vollen Suffix-Array.
    pub fn build(
        reference: Reference,
        variant: Variant,
        check_queries: u64,
    ) -> Result<(FmIndex, Option<CheckReport>), BuildError> {
        let l = reference.len();
        if l == 0 {
            return Err(BuildError::Empty);
        }
        let max = (sais::MAX_LEN as u64 - 1) / 2;
        if l > max {
            return Err(BuildError::TooLong { bases: l, max });
        }
        let (occ_interval, sa_interval) = variant.parameters();
        gc_log::info(
            "INDEX.BAU.PARAMETER",
            "Index-Variante",
            json!({
                "variante": variant.name(), "basen_je_zaehlblock": occ_interval,
                "sa_abstand": sa_interval, "referenz_basen": l,
                "contigs": reference.contigs.len(), "n_strecken": reference.holes.len(),
                "n_basen": reference.holes.iter().map(|h| h.len).sum::<u64>(),
            }),
        );

        let st = gc_log::stage("index.text", json!({}));
        let text = reference.index_text();
        st.finish(json!({ "zeichen": text.len() }));

        let st = gc_log::stage("index.suffix_array", json!({ "verfahren": "SA-IS" }));
        let sa = sais::suffix_array(&text, 5);
        st.finish(json!({ "zeilen": sa.len() }));

        if check_queries > 0 {
            let st = gc_log::stage("index.suffix_array.pruefung", json!({}));
            sais::verify(&text, &sa).map_err(|row| BuildError::SuffixArray { row })?;
            st.finish(json!({ "ergebnis": "richtig sortiert, jede Position genau einmal" }));
        }

        let st = gc_log::stage("index.bwt", json!({}));
        let mut builder = Builder::new(sa.len() as u64, occ_interval, sa_interval);
        for &p in &sa {
            let prev = if p == 0 {
                SENTINEL
            } else {
                text[p as usize - 1] - 1
            };
            builder.push(prev, u64::from(p));
        }
        let index = builder.finish(reference);
        st.finish(json!({
            "zeilen": index.rows, "primaerzeile": index.primary,
            "sa_stichproben": index.sa_samples.len() / SA_BYTES,
        }));

        let report = if check_queries > 0 {
            let st = gc_log::stage("index.selbstpruefung", json!({ "abfragen": check_queries }));
            let report = index
                .check_against_suffix_array(&text, &sa, check_queries)
                .map_err(BuildError::Check)?;
            st.finish(json!({
                "abfragen": report.queries, "abfragen_mit_treffern": report.queries_with_hits,
                "verglichene_fundstellen": report.occurrences_compared, "abweichungen": 0,
            }));
            Some(report)
        } else {
            None
        };
        Ok((index, report))
    }

    fn check_against_suffix_array(
        &self,
        text: &[u8],
        sa: &[u32],
        queries: u64,
    ) -> Result<CheckReport, String> {
        let mut rng = Rng(0xC0FF_EE12_3456_789B);
        let body = text.len() as u64 - 1;
        let mut report = CheckReport::default();
        let mut pattern = Vec::new();
        for q in 0..queries {
            // Länge 1 bis 150; Herkunft: Stück aus dem Text, davon ein Drittel mit einer Änderung,
            // jede sechzehnte Abfrage reiner Zufall.
            let len = 1 + rng.next() % 150;
            pattern.clear();
            if q % 16 == 15 || body < len {
                pattern.extend((0..len.min(40)).map(|_| (rng.next() & 3) as u8));
            } else {
                let start = rng.next() % (body - len + 1);
                pattern.extend(
                    text[start as usize..(start + len) as usize]
                        .iter()
                        .map(|c| c - 1),
                );
                if q % 3 == 2 {
                    let at = (rng.next() % len) as usize;
                    pattern[at] = (pattern[at] + 1 + (rng.next() % 3) as u8) & 3;
                }
            }
            let (lo, hi) = self.interval(&pattern);
            let (elo, ehi) = sa_interval_by_binary_search(text, sa, &pattern);
            if hi - lo != (ehi - elo) as u64 || (hi > lo && lo != elo as u64) {
                return Err(format!(
                    "Abfrage {q} (Länge {}): FM-Index Zeilen {lo}..{hi}, Suffix-Array {elo}..{ehi}",
                    pattern.len()
                ));
            }
            report.queries += 1;
            if hi > lo {
                report.queries_with_hits += 1;
            }
            // Fundstellen einzeln vergleichen (bei sehr häufigen Mustern die ersten 200).
            for row in lo..hi.min(lo + 200) {
                let got = self.locate(row);
                let want = u64::from(sa[row as usize]);
                if got != want {
                    return Err(format!("Zeile {row}: geortet {got}, Suffix-Array {want}"));
                }
                report.occurrences_compared += 1;
            }
            if gc_log::enabled(Level::Debug) && q % 100_000 == 0 {
                gc_log::debug(
                    "INDEX.PRUEFUNG.STAND",
                    "Selbstprüfung läuft",
                    json!({ "abfrage": q }),
                );
            }
        }
        Ok(report)
    }

    pub fn reference(&self) -> &Reference {
        &self.reference
    }

    pub fn rows(&self) -> u64 {
        self.rows
    }

    pub fn occ_interval(&self) -> u32 {
        self.occ_interval
    }

    pub fn sa_interval(&self) -> u32 {
        self.sa_interval
    }

    /// Speicherbedarf der Tabellen in Bytes (ohne Contig-Namen).
    pub fn table_bytes(&self) -> u64 {
        (self.blocks.len() * 8 + self.sa_samples.len() + self.reference.forward.words.len() * 8)
            as u64
    }

    #[inline]
    fn sa_sample(&self, i: usize) -> u64 {
        let b = &self.sa_samples[i * SA_BYTES..(i + 1) * SA_BYTES];
        u64::from(b[0])
            | u64::from(b[1]) << 8
            | u64::from(b[2]) << 16
            | u64::from(b[3]) << 24
            | u64::from(b[4]) << 32
    }

    /// Erste Zeile, deren Suffix mit Base `c` beginnt.
    #[inline]
    pub fn first_row(&self, c: u8) -> u64 {
        self.c_before[c as usize]
    }

    /// Anzahl der Base `c` in den BWT-Zeilen `0..i`.
    #[inline]
    pub fn occ(&self, c: u8, i: u64) -> u64 {
        let interval = u64::from(self.occ_interval);
        let stride = 4 + (self.occ_interval / 32) as usize;
        let base = (i >> self.occ_shift) as usize * stride;
        let within = (i & (interval - 1)) as u32;
        let block = &self.blocks[base..base + stride];
        let n = block[c as usize] + count_one(&block[4..], c, within);
        // Das Endzeichen ist als Code 0 gespeichert, zählt aber nicht als A.
        let block_start = i - u64::from(within);
        if c == 0 && self.primary >= block_start && self.primary < i {
            n - 1
        } else {
            n
        }
    }

    /// Anzahl jeder der vier Basen in den BWT-Zeilen `0..i`, in einem Durchgang durch den
    /// Zählblock (bwa: `bwt_occ4`). Dasselbe wie viermal `occ`.
    #[inline]
    pub fn occ4(&self, i: u64) -> [u64; 4] {
        if self.occ_interval == 128 && !self.occ_plain {
            let within = i & 127;
            let base = (i >> 7) as usize * 8;
            if within <= 64 || base + 16 <= self.blocks.len() {
                return self.occ4_128(i);
            }
        }
        self.occ4_plain(i)
    }

    /// Schaltet das Zählen ohne Verzweigung ab (`true`) oder wieder ein. Nur zum Vergleichen:
    /// Beide Wege geben dieselben Zahlen.
    pub fn set_occ_plain(&mut self, plain: bool) {
        self.occ_plain = plain;
    }

    /// `occ4` für Blöcke zu 128 Basen ohne Schleife und ohne Sprung nach der Blockhälfte.
    /// Vorn (bis 64 Basen in den Block) wird vom Zählstand dieses Blocks vorwärts gezählt,
    /// hinten vom Zählstand des nächsten Blocks rückwärts; welche Hälfte, entscheidet sich für
    /// jeden Zugriff anders und war für den Prozessor kaum vorherzusagen. Hier werden in beiden
    /// Fällen genau zwei Wörter mit einer Maske gezählt, die Wahl der Hälfte ist eine Auswahl
    /// von Werten statt eines Sprungs. Gleiche Zahlen wie `occ4_plain`. Im letzten Block gibt
    /// es keinen nächsten; dort ruft `occ4` für die hintere Hälfte `occ4_plain`.
    #[inline]
    fn occ4_128(&self, i: u64) -> [u64; 4] {
        let base = (i >> 7) as usize * 8;
        let within = (i & 127) as u32;
        let back = within > 64;
        let counts = base + 8 * usize::from(back);
        let words = base + 4 + 2 * usize::from(back);
        // Basen 0..x der beiden Wörter: vorn die gezählten, hinten die nicht gezählten.
        let x = within - 64 * u32::from(back);
        let below = u128::MAX.checked_shr(128 - 2 * x).unwrap_or(0);
        let flip = 0u128.wrapping_sub(u128::from(back));
        let mask = (below ^ flip) & (u128::from(EVEN) << 64 | u128::from(EVEN));
        let (m0, m1) = (mask as u64, (mask >> 64) as u64);
        let w: &[u64; 2] = self.blocks[words..words + 2].try_into().unwrap();
        let (w0, w1) = (w[0], w[1]);
        let (l0, h0) = (w0 & m0, (w0 >> 1) & m0);
        let (l1, h1) = (w1 & m1, (w1 >> 1) & m1);
        let (b0, b1) = (l0 & h0, l1 & h1);
        let c = count2(l0 ^ b0, l1 ^ b1);
        let g = count2(h0 ^ b0, h1 ^ b1);
        let t = count2(b0, b1);
        let n = if back { 64 - x } else { x };
        // Das Endzeichen ist als Code 0 gespeichert, zählt aber nicht als A.
        let block_start = i - u64::from(within);
        let (lo, hi) = if back {
            (i, block_start + 128)
        } else {
            (block_start, i)
        };
        let a = u64::from(n) - c - g - t - u64::from(self.primary >= lo && self.primary < hi);
        // Vorn wird addiert, hinten abgezogen: `sign` ist 0 oder alle Bits.
        let sign = 0u64.wrapping_sub(u64::from(back));
        let have: &[u64; 4] = self.blocks[counts..counts + 4].try_into().unwrap();
        let at = |j: usize, v: u64| have[j].wrapping_add((v ^ sign).wrapping_sub(sign));
        [at(0, a), at(1, c), at(2, g), at(3, t)]
    }

    /// `occ4` mit Schleifen über die Wörter des Blocks, für jeden Abstand der Zählstände.
    #[inline]
    fn occ4_plain(&self, i: u64) -> [u64; 4] {
        let interval = u64::from(self.occ_interval);
        let stride = 4 + (self.occ_interval / 32) as usize;
        let base = (i >> self.occ_shift) as usize * stride;
        let within = (i & (interval - 1)) as u32;
        let block = &self.blocks[base..base + stride];
        // Liegt `i` in der hinteren Hälfte des Blocks, wird vom Zählstand des nächsten Blocks
        // zurückgezählt: Dann sind es höchstens halb so viele Basen. Gleiches Ergebnis.
        if within > (interval / 2) as u32 && base + 2 * stride <= self.blocks.len() {
            let next = &self.blocks[base + stride..base + stride + 4];
            let [c, g, t] = count_cgt_tail(&block[4..], within);
            let block_end = i - u64::from(within) + interval;
            let mut a = (block_end - i) - c - g - t;
            // Das Endzeichen ist als Code 0 gespeichert, zählt aber nicht als A.
            if self.primary >= i && self.primary < block_end {
                a -= 1;
            }
            return [next[0] - a, next[1] - c, next[2] - g, next[3] - t];
        }
        let [c, g, t] = count_cgt(&block[4..], within);
        let mut a = u64::from(within) - c - g - t;
        // Das Endzeichen ist als Code 0 gespeichert, zählt aber nicht als A.
        let block_start = i - u64::from(within);
        if self.primary >= block_start && self.primary < i {
            a -= 1;
        }
        [block[0] + a, block[1] + c, block[2] + g, block[3] + t]
    }

    /// `occ4` für beide Enden eines Zeilenbereichs `[k, l)` (bwa: `bwt_2occ4`). Bei der
    /// Seed-Suche sind die Bereiche nach wenigen Basen eng, beide Enden liegen dann fast immer
    /// im selben Wort der BWT. In dem Fall wird nur bis `k` gezählt und das kurze Stück bis `l`
    /// dazu, statt den Block zweimal von vorn durchzugehen.
    #[inline]
    pub fn occ4_pair(&self, k: u64, l: u64) -> ([u64; 4], [u64; 4]) {
        let lo = self.occ4(k);
        if k >> 5 != l >> 5 || l < k {
            return (lo, self.occ4(l));
        }
        // Gleiches Wort, also auch derselbe Block.
        let interval = u64::from(self.occ_interval);
        let stride = 4 + (self.occ_interval / 32) as usize;
        let base = (k >> self.occ_shift) as usize * stride;
        let within = (k & (interval - 1)) as usize;
        let (from, to) = ((k & 31) as u32, (l & 31) as u32);
        let keep = EVEN & ((1u64 << (2 * to)) - 1) & !((1u64 << (2 * from)) - 1);
        let w = self.blocks[base + 4 + within / 32];
        let (low, high) = (w & keep, (w >> 1) & keep);
        let both = low & high;
        let (c, g, t) = (
            fold(pairs(low ^ both)),
            fold(pairs(high ^ both)),
            fold(pairs(both)),
        );
        let mut a = l - k - c - g - t;
        // Das Endzeichen ist als Code 0 gespeichert, zählt aber nicht als A.
        if self.primary >= k && self.primary < l {
            a -= 1;
        }
        (lo, [lo[0] + a, lo[1] + c, lo[2] + g, lo[3] + t])
    }

    /// Holt den Zählblock der Zeile `i` schon in den Cache, damit ein folgendes `occ4` an
    /// dieser Zeile nicht auf den Speicher warten muss. Ändert kein Ergebnis.
    #[inline]
    pub fn prefetch_row(&self, i: u64) {
        let stride = 4 + (self.occ_interval / 32) as usize;
        let base = (i >> self.occ_shift) as usize * stride;
        // Ein Block hat bis zu 64 Bytes und liegt nicht ausgerichtet: Anfang und Ende holen.
        gc_sys::prefetch(&self.blocks[base]);
        gc_sys::prefetch(&self.blocks[base + stride - 1]);
        // In der hinteren Hälfte zählt `occ4` vom Zählstand des nächsten Blocks zurück.
        let within = i & (u64::from(self.occ_interval) - 1);
        if within > u64::from(self.occ_interval / 2) && base + 2 * stride <= self.blocks.len() {
            gc_sys::prefetch(&self.blocks[base + stride]);
        }
    }

    #[inline]
    fn bwt(&self, row: u64) -> u8 {
        let interval = u64::from(self.occ_interval);
        let stride = 4 + (self.occ_interval / 32) as usize;
        let base = (row >> self.occ_shift) as usize * stride;
        let within = row & (interval - 1);
        (self.blocks[base + 4 + (within / 32) as usize] >> (2 * (within % 32)) & 3) as u8
    }

    /// Zeilenbereich `[lo, hi)` aller Suffixe, die mit `pattern` (Codes 0..3) beginnen.
    pub fn interval(&self, pattern: &[u8]) -> (u64, u64) {
        let (mut lo, mut hi) = (0, self.rows);
        for &c in pattern.iter().rev() {
            lo = self.c_before[c as usize] + self.occ(c, lo);
            hi = self.c_before[c as usize] + self.occ(c, hi);
            if lo >= hi {
                return (lo, lo);
            }
        }
        (lo, hi)
    }

    /// Anzahl der Fundstellen auf beiden Strängen, ohne Filter.
    pub fn count(&self, pattern: &[u8]) -> u64 {
        let (lo, hi) = self.interval(pattern);
        hi - lo
    }

    /// Textposition des Suffixes in dieser Zeile.
    pub fn locate(&self, mut row: u64) -> u64 {
        let mut steps = 0;
        let interval = u64::from(self.sa_interval);
        loop {
            if row == self.primary {
                return steps;
            }
            if row & (interval - 1) == 0 {
                return self.sa_sample((row >> self.sa_shift) as usize) + steps;
            }
            let c = self.bwt(row);
            row = self.c_before[c as usize] + self.occ(c, row);
            steps += 1;
        }
    }

    /// Schaltet einen Zwischenspeicher für [`FmIndex::locate_many`] ein: `2^bits` Plätze zu je 8
    /// Bytes, `bits = 0` schaltet ihn aus. Reads aus Wiederholungen (Alu, L1) finden ihre Seeds
    /// immer wieder in denselben Zeilen des Index; jede Zeile kostet ohne Speicher bis zu
    /// `sa_interval` Schritte durch die BWT, fast jeder davon außerhalb des Prozessor-Caches. Ein
    /// Platz hält Zeile und Position zusammen in einem Wort; ein Treffer liefert genau den Wert,
    /// den der Weg durch die BWT ergäbe. Mehrere Kerne lesen und schreiben ohne Sperre.
    /// Geht nur, wenn Zeilen und Positionen in 33 Bit passen (bis 4,2 Mrd. Basen Referenz);
    /// sonst bleibt er aus.
    pub fn enable_locate_cache(&mut self, bits: u32) {
        self.locate_cache = LocateCache(
            if bits == 0 || bits > 30 || self.rows >= 1 << CACHE_POS_BITS {
                Vec::new()
            } else {
                (0..1usize << bits)
                    .map(|_| std::sync::atomic::AtomicU64::new(0))
                    .collect()
            },
        );
    }

    fn cache_get(&self, row: u64) -> Option<u64> {
        let n = self.locate_cache.0.len() as u64;
        if n == 0 {
            return None;
        }
        let e = self.locate_cache.0[(row & (n - 1)) as usize]
            .load(std::sync::atomic::Ordering::Relaxed);
        let tag = row >> n.trailing_zeros();
        (e >> 63 == 1 && (e >> CACHE_POS_BITS) & CACHE_TAG_MASK == tag)
            .then_some(e & CACHE_POS_MASK)
    }

    fn cache_put(&self, row: u64, pos: u64) {
        let n = self.locate_cache.0.len() as u64;
        if n == 0 {
            return;
        }
        let tag = row >> n.trailing_zeros();
        let e = 1 << 63 | tag << CACHE_POS_BITS | pos;
        self.locate_cache.0[(row & (n - 1)) as usize]
            .store(e, std::sync::atomic::Ordering::Relaxed);
    }

    /// Textpositionen vieler Zeilen auf einmal, angehängt an `out` in der Reihenfolge von `rows`.
    /// Dasselbe Ergebnis wie `locate` je Zeile, aber die Wege mehrerer Zeilen laufen im
    /// Gleichschritt: In jeder Runde wird erst für alle offenen Zeilen der Block der BWT gelesen,
    /// dann für alle gezählt. Die Zugriffe einer Runde hängen nicht voneinander ab, der Prozessor
    /// kann also auf mehrere Speicherzeilen gleichzeitig warten statt nacheinander. Bei einem
    /// ganzen Genom liegt fast jeder Schritt außerhalb des Zwischenspeichers, dort zählt das.
    pub fn locate_many(&self, rows: &[u64], out: &mut Vec<u64>) {
        const LANES: usize = 32;
        let start = out.len();
        out.resize(start + rows.len(), 0);
        let out = &mut out[start..];
        let interval = u64::from(self.sa_interval);
        let mut row = [0u64; LANES];
        let mut first = [0u64; LANES];
        let mut steps = [0u64; LANES];
        let mut slot = [0usize; LANES];
        let mut code = [0u8; LANES];
        let (mut live, mut next) = (0, 0);
        loop {
            while live < LANES && next < rows.len() {
                if let Some(pos) = self.cache_get(rows[next]) {
                    out[next] = pos;
                } else {
                    (row[live], first[live], steps[live], slot[live]) =
                        (rows[next], rows[next], 0, next);
                    live += 1;
                }
                next += 1;
            }
            if live == 0 {
                break;
            }
            // Fertige Zeilen abschließen; die letzte offene rückt auf den frei gewordenen Platz.
            let mut i = 0;
            while i < live {
                let r = row[i];
                let done = if r == self.primary {
                    Some(steps[i])
                } else if r & (interval - 1) == 0 {
                    Some(self.sa_sample((r >> self.sa_shift) as usize) + steps[i])
                } else {
                    None
                };
                match done {
                    Some(pos) => {
                        out[slot[i]] = pos;
                        self.cache_put(first[i], pos);
                        live -= 1;
                        (row[i], first[i], steps[i], slot[i]) =
                            (row[live], first[live], steps[live], slot[live]);
                    }
                    None => i += 1,
                }
            }
            for i in 0..live {
                code[i] = self.bwt(row[i]);
            }
            for i in 0..live {
                row[i] = self.c_before[code[i] as usize] + self.occ(code[i], row[i]);
                steps[i] += 1;
            }
        }
    }

    /// Basen der Textpositionen `[from, to)` (Vorwärtsstrang, dann reverses Komplement). Stellen,
    /// an denen die Referenz N hat, kommen als Code 4 (unbekannt), nicht als die Füllbase, die
    /// der Index dort für die Suche führt.
    pub fn text_range(&self, from: u64, to: u64) -> Vec<u8> {
        let mut bases: Vec<u8> = (from..to).map(|p| self.text_base(p)).collect();
        let l = self.reference.len();
        // Der Bereich kann auf beiden Strängen liegen; je Strang das Stück des Vorwärtsstrangs.
        let forward = (from.min(l), to.min(l));
        let reverse = (2 * l - to.max(l), 2 * l - from.max(l));
        let holes = &self.reference.holes;
        for (strand, (a, b)) in [forward, reverse].into_iter().enumerate() {
            if a >= b {
                continue;
            }
            let first = holes.partition_point(|h| h.offset + h.len <= a);
            for h in holes[first..].iter().take_while(|h| h.offset < b) {
                for pos in h.offset.max(a)..(h.offset + h.len).min(b) {
                    let text = if strand == 0 { pos } else { 2 * l - 1 - pos };
                    bases[(text - from) as usize] = 4;
                }
            }
        }
        bases
    }

    /// Base an Textposition `p` (Vorwärtsstrang, dann reverses Komplement), aus der Referenz gelesen.
    #[inline]
    pub fn text_base(&self, p: u64) -> u8 {
        let l = self.reference.len();
        if p < l {
            self.reference.forward.get(p)
        } else {
            3 - self.reference.forward.get(2 * l - 1 - p)
        }
    }

    /// Vollprüfung ohne Suffix-Array: geht die BWT einmal rückwärts durch. Jeder Schritt muss die
    /// nächste Base der Referenz liefern (von hinten nach vorn, beide Stränge), jede Stichprobe
    /// des Suffix-Arrays muss die erreichte Textposition nennen, und nach genau so vielen Schritten,
    /// wie der Text Zeichen hat, muss die Primärzeile erreicht sein. `progress` wird gelegentlich
    /// mit der Zahl der geprüften Zeichen gerufen.
    pub fn verify_by_inversion(&self, mut progress: impl FnMut(u64)) -> Result<u64, String> {
        let n = self.rows - 1;
        let sa_interval = u64::from(self.sa_interval);
        let (mut row, mut p) = (0u64, n);
        if self.sa_sample(0) != n {
            return Err(format!(
                "Zeile 0: Stichprobe {} statt {n}",
                self.sa_sample(0)
            ));
        }
        while p > 0 {
            if row == self.primary {
                return Err(format!(
                    "Primärzeile zu früh erreicht, noch {p} Zeichen offen"
                ));
            }
            let c = self.bwt(row);
            p -= 1;
            if c != self.text_base(p) {
                return Err(format!(
                    "Textposition {p}: BWT liefert Base {c}, Referenz {}",
                    self.text_base(p)
                ));
            }
            row = self.c_before[c as usize] + self.occ(c, row);
            if row.is_multiple_of(sa_interval) && self.sa_sample((row / sa_interval) as usize) != p
            {
                return Err(format!(
                    "Zeile {row}: Stichprobe {} statt {p}",
                    self.sa_sample((row / sa_interval) as usize)
                ));
            }
            if p.is_multiple_of(1 << 24) {
                progress(n - p);
            }
        }
        if row != self.primary {
            return Err(format!(
                "am Textanfang in Zeile {row}, Primärzeile ist {}",
                self.primary
            ));
        }
        Ok(n)
    }

    /// Stichprobenprüfung ohne Suffix-Array: Stücke aus der Referenz müssen an ihrer Herkunft
    /// gefunden werden, und jede geortete Fundstelle muss in der Referenz tatsächlich das Muster tragen.
    pub fn check_against_reference(&self, queries: u64) -> Result<CheckReport, String> {
        let n = self.rows - 1;
        let mut rng = Rng(0x5EED_1234_ABCD_EF01);
        let mut report = CheckReport::default();
        let mut pattern = Vec::new();
        for q in 0..queries {
            let len = (12 + rng.next() % 139).min(n);
            let start = rng.next() % (n - len + 1);
            pattern.clear();
            pattern.extend((start..start + len).map(|p| self.text_base(p)));
            if q % 3 == 2 {
                let at = (rng.next() % len) as usize;
                pattern[at] = (pattern[at] + 1 + (rng.next() % 3) as u8) & 3;
            }
            let (lo, hi) = self.interval(&pattern);
            report.queries += 1;
            if hi > lo {
                report.queries_with_hits += 1;
            }
            let mut found_origin = q % 3 == 2;
            for row in lo..hi.min(lo + 200) {
                let at = self.locate(row);
                if at + len > n || (0..len).any(|i| self.text_base(at + i) != pattern[i as usize]) {
                    return Err(format!(
                        "Abfrage {q}: Zeile {row} ortet {at}, dort steht das Muster nicht"
                    ));
                }
                found_origin |= at == start;
                report.occurrences_compared += 1;
            }
            if !found_origin && hi - lo <= 200 {
                return Err(format!(
                    "Abfrage {q}: Herkunft {start} (Länge {len}) nicht unter den {} Fundstellen",
                    hi - lo
                ));
            }
        }
        Ok(report)
    }

    /// Bereich des Contigs im Indextext auf dem Strang, in dem `text_pos` liegt: `[von, bis)`.
    /// `None`, wenn die Position hinter dem Text liegt.
    pub fn contig_text_range(&self, text_pos: u64) -> Option<(u32, u64, u64)> {
        let l = self.reference.len();
        if text_pos >= 2 * l {
            return None;
        }
        let forward = if text_pos < l {
            text_pos
        } else {
            2 * l - 1 - text_pos
        };
        let contig = self.reference.contig_at(forward)?;
        let c = &self.reference.contigs[contig];
        let (beg, end) = if text_pos < l {
            (c.offset, c.offset + c.len)
        } else {
            (2 * l - c.offset - c.len, 2 * l - c.offset)
        };
        Some((contig as u32, beg, end))
    }

    /// Wandelt eine Textposition in einen Treffer im Genom um. `None`, wenn die Fundstelle über
    /// die Grenze zwischen den Strängen oder zwischen zwei Contigs reicht oder eine N-Strecke berührt.
    pub fn hit_at(&self, text_pos: u64, len: u64) -> Option<Hit> {
        let l = self.reference.len();
        let (start, reverse) = if text_pos + len <= l {
            (text_pos, false)
        } else if text_pos >= l && text_pos + len <= 2 * l {
            (2 * l - text_pos - len, true)
        } else {
            return None;
        };
        let contig = self.reference.contig_at(start)?;
        let c = &self.reference.contigs[contig];
        if len == 0 || start + len > c.offset + c.len || self.reference.overlaps_hole(start, len) {
            return None;
        }
        Some(Hit {
            contig: contig as u32,
            pos: start - c.offset,
            reverse,
        })
    }

    /// Alle gültigen Treffer eines Musters, sortiert. `limit` begrenzt die georteten Zeilen.
    pub fn hits(&self, pattern: &[u8], limit: u64) -> Vec<Hit> {
        let (lo, hi) = self.interval(pattern);
        let mut out: Vec<Hit> = (lo..hi.min(lo.saturating_add(limit)))
            .filter_map(|row| self.hit_at(self.locate(row), pattern.len() as u64))
            .collect();
        out.sort_unstable();
        out
    }

    // ---- Datei ----

    /// Schreibt den Index. Alle Zahlen little-endian, am Ende eine CRC-32 über alles davor.
    pub fn write<W: Write>(&self, out: W) -> io::Result<()> {
        let mut w = CrcWriter {
            inner: out,
            crc: flate2::Crc::new(),
        };
        w.write_all(MAGIC)?;
        w.write_all(&FORMAT.to_le_bytes())?;
        for v in [
            self.rows,
            self.primary,
            u64::from(self.occ_interval),
            u64::from(self.sa_interval),
        ] {
            w.write_all(&v.to_le_bytes())?;
        }
        for v in self.c_before {
            w.write_all(&v.to_le_bytes())?;
        }
        let r = &self.reference;
        w.write_all(&(r.contigs.len() as u64).to_le_bytes())?;
        for c in &r.contigs {
            w.write_all(&(c.name.len() as u64).to_le_bytes())?;
            w.write_all(c.name.as_bytes())?;
            w.write_all(&c.len.to_le_bytes())?;
            w.write_all(&c.offset.to_le_bytes())?;
        }
        w.write_all(&(r.holes.len() as u64).to_le_bytes())?;
        for h in &r.holes {
            w.write_all(&h.offset.to_le_bytes())?;
            w.write_all(&h.len.to_le_bytes())?;
        }
        w.write_all(&r.forward.len.to_le_bytes())?;
        write_u64s(&mut w, &r.forward.words)?;
        write_u64s(&mut w, &self.blocks)?;
        w.write_all(&((self.sa_samples.len() / SA_BYTES) as u64).to_le_bytes())?;
        for chunk in self.sa_samples.chunks(1 << 16) {
            w.write_all(chunk)?;
        }
        let crc = w.crc.sum();
        w.inner.write_all(&crc.to_le_bytes())?;
        w.inner.flush()
    }

    /// Liest nur den Kopf einer Indexdatei: Zeilenzahl und Abstand der Suffix-Array-Stichproben.
    /// Daraus folgt der Speicher der Stichproben (Zeilen / Abstand × 5 Bytes), ohne die Datei
    /// zu laden.
    pub fn read_sample_layout<R: Read>(mut input: R) -> io::Result<(u64, u32)> {
        let mut head = [0u8; 40];
        input.read_exact(&mut head)?;
        if &head[..4] != MAGIC || u32::from_le_bytes(head[4..8].try_into().unwrap()) != FORMAT {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Indexdatei: keine Genome-Indexdatei dieser Version",
            ));
        }
        let rows = u64::from_le_bytes(head[8..16].try_into().unwrap());
        let sa = u64::from_le_bytes(head[32..40].try_into().unwrap());
        if sa == 0 || sa > 1024 || !sa.is_power_of_two() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Indexdatei: unzulässige Abstände",
            ));
        }
        Ok((rows, sa as u32))
    }

    /// Bytes je Suffix-Array-Stichprobe.
    pub const SAMPLE_BYTES: u64 = SA_BYTES as u64;

    /// Liest einen Index und prüft Format, Größen und Prüfsumme.
    pub fn read<R: Read>(input: R) -> io::Result<FmIndex> {
        Self::read_thinned(input, 1)
    }

    /// Wie [`FmIndex::read`], behält aber nur jede `keep_every`-te Stichprobe des
    /// Suffix-Arrays (eine Zweierpotenz). Der Abstand der Stichproben wächst um diesen Faktor,
    /// der Speicher für sie sinkt um ihn. Alle Ergebnisse bleiben gleich, nur das Auffinden
    /// der Textpositionen braucht im Mittel mehr Schritte. Ist der Abstand in der Datei schon
    /// 1024, bleibt es dabei.
    pub fn read_thinned<R: Read>(input: R, keep_every: u32) -> io::Result<FmIndex> {
        let bad = |m: &str| io::Error::new(io::ErrorKind::InvalidData, format!("Indexdatei: {m}"));
        let mut r = CrcReader {
            inner: input,
            crc: flate2::Crc::new(),
        };
        let mut magic = [0u8; 4];
        r.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(bad("keine Genome-Indexdatei (Kennung fehlt)"));
        }
        let format = read_u32(&mut r)?;
        if format != FORMAT {
            return Err(bad(&format!("Formatversion {format}, erwartet {FORMAT}")));
        }
        let rows = read_u64(&mut r)?;
        let primary = read_u64(&mut r)?;
        let occ_interval = read_u64(&mut r)?;
        let sa_interval = read_u64(&mut r)?;
        if !(32..=u64::from(MAX_OCC_INTERVAL)).contains(&occ_interval)
            || !occ_interval.is_power_of_two()
            || !(1..=1024).contains(&sa_interval)
            || !sa_interval.is_power_of_two()
        {
            return Err(bad("unzulässige Abstände"));
        }
        if rows == 0 || rows >= 1 << 40 || primary >= rows {
            return Err(bad("unzulässige Zeilenzahl"));
        }
        let mut c_before = [0u64; 4];
        for c in &mut c_before {
            *c = read_u64(&mut r)?;
        }
        let n_contigs = read_u64(&mut r)?;
        if n_contigs > 1 << 24 {
            return Err(bad("zu viele Contigs"));
        }
        let mut contigs = Vec::with_capacity(n_contigs as usize);
        for _ in 0..n_contigs {
            let name_len = read_u64(&mut r)?;
            if name_len > 1 << 16 {
                return Err(bad("Contig-Name zu lang"));
            }
            let mut name = vec![0u8; name_len as usize];
            r.read_exact(&mut name)?;
            let name = String::from_utf8(name).map_err(|_| bad("Contig-Name ist kein UTF-8"))?;
            contigs.push(Contig {
                name,
                len: read_u64(&mut r)?,
                offset: read_u64(&mut r)?,
            });
        }
        let n_holes = read_u64(&mut r)?;
        if n_holes > rows {
            return Err(bad("zu viele N-Strecken"));
        }
        let mut holes = Vec::with_capacity(n_holes as usize);
        for _ in 0..n_holes {
            holes.push(Hole {
                offset: read_u64(&mut r)?,
                len: read_u64(&mut r)?,
            });
        }
        let l = read_u64(&mut r)?;
        if 2 * l + 1 != rows {
            return Err(bad("Referenzlänge passt nicht zur Zeilenzahl"));
        }
        let words = read_u64s(&mut r, l.div_ceil(32))?;
        let stride = 4 + occ_interval / 32;
        let blocks = read_u64s(&mut r, (rows / occ_interval + 1) * stride)?;
        let n_samples = read_u64(&mut r)?;
        if n_samples != rows.div_ceil(sa_interval) {
            return Err(bad("Zahl der Suffix-Array-Stichproben passt nicht"));
        }
        if keep_every == 0 || !keep_every.is_power_of_two() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Ausdünnung der Stichproben muss eine Zweierpotenz sein",
            ));
        }
        let keep = u64::from(keep_every).min(1024 / sa_interval).max(1);
        let sa_interval = sa_interval * keep;
        let mut sa_samples = vec![0u8; rows.div_ceil(sa_interval) as usize * SA_BYTES];
        if keep == 1 {
            for chunk in sa_samples.chunks_mut(1 << 16) {
                r.read_exact(chunk)?;
            }
        } else {
            // Stückweise lesen, damit nie alle Stichproben zugleich im Speicher liegen.
            let per_chunk = (keep as usize) * 8192;
            let mut buf = vec![0u8; per_chunk * SA_BYTES];
            let (mut left, mut seen, mut out) = (n_samples as usize, 0usize, 0usize);
            while left > 0 {
                let take = left.min(per_chunk);
                r.read_exact(&mut buf[..take * SA_BYTES])?;
                for i in 0..take {
                    if (seen + i) % keep as usize == 0 {
                        sa_samples[out..out + SA_BYTES]
                            .copy_from_slice(&buf[i * SA_BYTES..(i + 1) * SA_BYTES]);
                        out += SA_BYTES;
                    }
                }
                seen += take;
                left -= take;
            }
            debug_assert_eq!(out, sa_samples.len());
        }
        let computed = r.crc.sum();
        let mut stored = [0u8; 4];
        r.inner.read_exact(&mut stored)?;
        if u32::from_le_bytes(stored) != computed {
            return Err(bad(
                "Prüfsumme stimmt nicht, Datei ist beschädigt oder unvollständig",
            ));
        }
        let mut reference = Reference::new();
        reference.contigs = contigs;
        reference.holes = holes;
        reference.forward = Packed { len: l, words };
        let mut index = FmIndex {
            rows,
            primary,
            c_before,
            occ_interval: occ_interval as u32,
            sa_interval: sa_interval as u32,
            occ_shift: occ_interval.trailing_zeros(),
            sa_shift: sa_interval.trailing_zeros(),
            blocks,
            sa_samples,
            reference,
            short: crate::smem::ShortTable::default(),
            locate_cache: LocateCache::default(),
            occ_plain: false,
        };
        index.short = crate::smem::ShortTable::build(&index);
        Ok(index)
    }
}

struct CrcWriter<W> {
    inner: W,
    crc: flate2::Crc,
}

impl<W: Write> CrcWriter<W> {
    fn write_all(&mut self, data: &[u8]) -> io::Result<()> {
        self.crc.update(data);
        self.inner.write_all(data)
    }
}

struct CrcReader<R> {
    inner: R,
    crc: flate2::Crc,
}

impl<R: Read> CrcReader<R> {
    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        self.inner.read_exact(buf)?;
        self.crc.update(buf);
        Ok(())
    }
}

fn write_u64s<W: Write>(w: &mut CrcWriter<W>, values: &[u64]) -> io::Result<()> {
    w.write_all(&(values.len() as u64).to_le_bytes())?;
    let mut buf = Vec::with_capacity(1 << 16);
    for chunk in values.chunks(1 << 13) {
        buf.clear();
        for v in chunk {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        w.write_all(&buf)?;
    }
    Ok(())
}

fn read_u32<R: Read>(r: &mut CrcReader<R>) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64<R: Read>(r: &mut CrcReader<R>) -> io::Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

fn read_u64s<R: Read>(r: &mut CrcReader<R>, expected: u64) -> io::Result<Vec<u64>> {
    let n = read_u64(r)?;
    if n != expected {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Indexdatei: Tabelle hat {n} Einträge, erwartet {expected}"),
        ));
    }
    let mut out = Vec::with_capacity(n as usize);
    let mut buf = vec![0u8; 1 << 16];
    let mut left = n as usize;
    while left > 0 {
        let take = left.min(buf.len() / 8);
        r.read_exact(&mut buf[..take * 8])?;
        out.extend(
            buf[..take * 8]
                .as_chunks::<8>()
                .0
                .iter()
                .map(|b| u64::from_le_bytes(*b)),
        );
        left -= take;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::encode;

    fn codes(s: &[u8]) -> Vec<u8> {
        s.iter().map(|&b| encode(b).unwrap()).collect()
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|c| 3 - c).collect()
    }

    /// Unabhängige Wahrheit: jedes Fenster jedes Contigs auf beiden Strängen vergleichen.
    fn naive_hits(contigs: &[(&str, Vec<u8>)], pattern: &[u8]) -> Vec<Hit> {
        let rc = revcomp(pattern);
        let mut out = Vec::new();
        for (ci, (_, seq)) in contigs.iter().enumerate() {
            let enc: Vec<Option<u8>> = seq.iter().map(|&b| encode(b)).collect();
            if enc.len() < pattern.len() {
                continue;
            }
            for pos in 0..=enc.len() - pattern.len() {
                let window = &enc[pos..pos + pattern.len()];
                if window.iter().zip(pattern).all(|(a, b)| *a == Some(*b)) {
                    out.push(Hit {
                        contig: ci as u32,
                        pos: pos as u64,
                        reverse: false,
                    });
                }
                if window.iter().zip(&rc).all(|(a, b)| *a == Some(*b)) {
                    out.push(Hit {
                        contig: ci as u32,
                        pos: pos as u64,
                        reverse: true,
                    });
                }
            }
        }
        out.sort_unstable();
        out
    }

    fn random_genome(seed: u64) -> Vec<(&'static str, Vec<u8>)> {
        let mut rng = Rng(seed);
        let mut make = |len: usize| -> Vec<u8> {
            let mut s: Vec<u8> = (0..len)
                .map(|_| b"ACGT"[(rng.next() & 3) as usize])
                .collect();
            // Wiederholungen, Homopolymer, N-Strecken, Kleinbuchstaben
            if len > 400 {
                let copy = s[50..170].to_vec();
                s[200..320].copy_from_slice(&copy);
                for b in &mut s[330..360] {
                    *b = b'A';
                }
                for b in &mut s[370..380] {
                    *b = b'N';
                }
                s[385] = b'n';
                s[390] = b'c';
            }
            s
        };
        vec![
            ("eins", make(1500)),
            ("zwei", make(700)),
            ("drei", make(33)),
            ("vier", make(1024)),
        ]
    }

    fn build(contigs: &[(&str, Vec<u8>)], variant: Variant) -> FmIndex {
        let mut r = Reference::new();
        for (name, seq) in contigs {
            r.push_contig(name, seq);
        }
        let (index, report) = FmIndex::build(r, variant, 3000).expect("Bau mit Selbstprüfung");
        assert_eq!(report.unwrap().queries, 3000);
        index
    }

    #[test]
    fn hits_equal_naive_search_in_both_variants() {
        let contigs = random_genome(7);
        for variant in [Variant::Standard, Variant::Kompakt, Variant::Schnell] {
            let index = build(&contigs, variant);
            let mut rng = Rng(99);
            let mut with_hits = 0;
            for q in 0..1500 {
                let (_, seq) = &contigs[(rng.next() % 4) as usize];
                let len = 1 + (rng.next() % 40) as usize;
                let pattern: Vec<u8> = if q % 5 == 0 || seq.len() < len {
                    (0..len.min(8)).map(|_| (rng.next() & 3) as u8).collect()
                } else {
                    let start = (rng.next() as usize) % (seq.len() - len + 1);
                    let piece = &seq[start..start + len];
                    if piece.iter().any(|&b| encode(b).is_none()) {
                        continue;
                    }
                    let p = codes(piece);
                    if q % 2 == 0 { revcomp(&p) } else { p }
                };
                let got = index.hits(&pattern, u64::MAX);
                let want = naive_hits(&contigs, &pattern);
                assert_eq!(got, want, "{variant:?} Abfrage {q} Muster {pattern:?}");
                with_hits += usize::from(!want.is_empty());
            }
            assert!(with_hits > 1000);
        }
    }

    #[test]
    fn text_range_marks_n_as_unknown_on_both_strands() {
        let contigs = vec![("a", b"ACGTNNACGTAC".to_vec()), ("b", b"NTTGCA".to_vec())];
        let index = build(&contigs, Variant::Kompakt);
        let l = index.reference().len();
        let forward: Vec<u8> = contigs.iter().flat_map(|c| c.1.clone()).collect();
        let code = |b: u8| match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        };
        let mut text: Vec<u8> = forward.iter().map(|&b| code(b)).collect();
        let rc: Vec<u8> = text
            .iter()
            .rev()
            .map(|&c| if c < 4 { 3 - c } else { 4 })
            .collect();
        text.extend(rc);
        for from in 0..=2 * l {
            for to in from..=2 * l {
                assert_eq!(
                    index.text_range(from, to),
                    text[from as usize..to as usize],
                    "[{from}, {to})"
                );
            }
        }
    }

    #[test]
    fn occ_matches_counting_and_edges() {
        let contigs = random_genome(3);
        let index = build(&contigs, Variant::Kompakt);
        // occ gegen direktes Zählen der BWT, an jeder Zeile
        let mut counts = [0u64; 4];
        for row in 0..=index.rows() {
            for c in 0..4u8 {
                assert_eq!(index.occ(c, row), counts[c as usize], "occ({c}, {row})");
            }
            assert_eq!(index.occ4(row), counts, "occ4({row})");
            if row < index.rows() && row != index.primary {
                counts[index.bwt(row) as usize] += 1;
            }
        }
        assert_eq!(counts.iter().sum::<u64>(), index.rows() - 1);
        // occ4 in allen Varianten (es zählt je nach Lage vom Anfang oder vom Ende des Blocks)
        // und auf beiden Wegen (ohne Verzweigung und mit Schleifen)
        for (variant, plain) in [Variant::Standard, Variant::Kompakt, Variant::Schnell]
            .into_iter()
            .flat_map(|v| [(v, false), (v, true)])
        {
            let mut index = build(&contigs, variant);
            index.set_occ_plain(plain);
            let mut counts = [0u64; 4];
            for row in 0..=index.rows() {
                assert_eq!(index.occ4(row), counts, "{variant:?} {plain} occ4({row})");
                if row < index.rows() && row != index.primary {
                    counts[index.bwt(row) as usize] += 1;
                }
            }
        }
        // beide Enden eines Bereichs auf einmal, in beiden Varianten, auch über das Endzeichen
        for variant in [Variant::Standard, Variant::Kompakt, Variant::Schnell] {
            let index = build(&contigs, variant);
            for k in 0..=index.rows() {
                for l in (k..=index.rows()).take(140) {
                    assert_eq!(
                        index.occ4_pair(k, l),
                        (index.occ4(k), index.occ4(l)),
                        "occ4_pair({k}, {l})"
                    );
                }
            }
        }
        // leeres Muster trifft jede Zeile; ein Muster länger als das Genom trifft nie
        assert_eq!(index.count(&[]), index.rows());
        assert_eq!(index.count(&vec![0u8; 10_000]), 0);
        assert!(index.hits(&[], 10).is_empty());
    }

    #[test]
    fn locate_many_equals_locate() {
        for variant in [Variant::Standard, Variant::Kompakt, Variant::Schnell] {
            let index = build(&random_genome(5), variant);
            let mut rng = Rng(9);
            // Alle Zeilen der Reihe nach, dann zufällige Mengen mit Wiederholungen, auch leere.
            let all: Vec<u64> = (0..index.rows()).collect();
            let mut sets = vec![all, Vec::new()];
            for n in [1usize, 2, 15, 16, 17, 100] {
                sets.push((0..n).map(|_| rng.next() % index.rows()).collect());
            }
            sets.push(vec![index.primary; 40]);
            for rows in sets {
                let mut got = vec![7u64]; // vorhandener Inhalt bleibt stehen
                index.locate_many(&rows, &mut got);
                let mut want = vec![7u64];
                want.extend(rows.iter().map(|&r| index.locate(r)));
                assert_eq!(got, want);
            }
        }
    }

    #[test]
    fn the_locate_cache_gives_the_same_positions() {
        for variant in [Variant::Kompakt, Variant::Schnell] {
            let plain = build(&random_genome(5), variant);
            // Kleiner Speicher mit 2^4 Plätzen: viele Zeilen teilen sich einen Platz.
            for bits in [4, 12] {
                let mut index = plain.clone();
                index.enable_locate_cache(bits);
                let mut rng = Rng(11);
                for round in 0..3 {
                    let rows: Vec<u64> = (0..500).map(|_| rng.next() % index.rows()).collect();
                    // Zweimal dieselben Zeilen: beim zweiten Mal kommen viele aus dem Speicher.
                    for _ in 0..2 {
                        let mut got = Vec::new();
                        index.locate_many(&rows, &mut got);
                        let want: Vec<u64> = rows.iter().map(|&r| plain.locate(r)).collect();
                        assert_eq!(got, want, "{variant:?}, 2^{bits} Plätze, Runde {round}");
                    }
                }
            }
        }
    }

    #[test]
    fn thinned_samples_find_the_same_positions() {
        let contigs = random_genome(12);
        let index = build(&contigs, Variant::Schnell);
        let mut bytes = Vec::new();
        index.write(&mut bytes).unwrap();
        for keep in [1u32, 2, 4, 16, 4096] {
            let thin = FmIndex::read_thinned(&bytes[..], keep).unwrap();
            let want = (index.sa_interval() * keep).min(1024.max(index.sa_interval()));
            assert_eq!(thin.sa_interval(), want);
            let rows: Vec<u64> = (0..index.rows()).collect();
            let (mut a, mut b) = (Vec::new(), Vec::new());
            index.locate_many(&rows, &mut a);
            thin.locate_many(&rows, &mut b);
            assert_eq!(a, b, "keep {keep}");
        }
        assert!(FmIndex::read_thinned(&bytes[..], 3).is_err());
        let (rows, sa) = FmIndex::read_sample_layout(&bytes[..]).unwrap();
        assert_eq!((rows, sa), (index.rows(), index.sa_interval()));
    }

    #[test]
    fn file_roundtrip_and_damage() {
        let contigs = random_genome(11);
        let index = build(&contigs, Variant::Schnell);
        let mut bytes = Vec::new();
        index.write(&mut bytes).unwrap();
        let back = FmIndex::read(&bytes[..]).unwrap();
        assert_eq!(back, index);
        let pattern = codes(&contigs[0].1[10..40]);
        assert_eq!(back.hits(&pattern, 100), index.hits(&pattern, 100));

        // ein gekipptes Bit irgendwo in den Tabellen fällt auf
        let mut damaged = bytes.clone();
        let at = damaged.len() / 2;
        damaged[at] ^= 0x10;
        let err = FmIndex::read(&damaged[..]).unwrap_err().to_string();
        assert!(
            err.contains("Prüfsumme") || err.contains("Indexdatei"),
            "{err}"
        );
        // abgeschnittene Datei
        assert!(FmIndex::read(&bytes[..bytes.len() - 9]).is_err());
        // falsche Kennung
        assert!(
            FmIndex::read(&b"BAM\x01aaaaaaaaaaaaaaaa"[..])
                .unwrap_err()
                .to_string()
                .contains("Kennung")
        );
    }

    #[test]
    fn inversion_check_passes_and_detects_damage() {
        let contigs = random_genome(21);
        for variant in [Variant::Standard, Variant::Kompakt, Variant::Schnell] {
            let mut index = build(&contigs, variant);
            assert_eq!(index.verify_by_inversion(|_| {}), Ok(index.rows() - 1));
            assert!(index.check_against_reference(500).is_ok());
            // eine Base der BWT ändern
            let stride = 4 + (index.occ_interval() / 32) as usize;
            index.blocks[3 * stride + 4] ^= 1 << 10;
            assert!(index.verify_by_inversion(|_| {}).is_err());
            index.blocks[3 * stride + 4] ^= 1 << 10;
            // eine Stichprobe ändern
            index.sa_samples[SA_BYTES * 17] ^= 1;
            assert!(index.verify_by_inversion(|_| {}).is_err());
        }
    }

    #[test]
    fn build_errors() {
        assert!(matches!(
            FmIndex::build(Reference::new(), Variant::Kompakt, 0),
            Err(BuildError::Empty)
        ));
    }
}
