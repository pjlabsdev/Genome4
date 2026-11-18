//! FM-Index: BWT mit Zähltabelle in Blöcken und abgetastetem Suffix-Array.
//!
//! Zwei Varianten desselben Formats:
//! - `Kompakt`: Zählstand alle 128 Basen, Suffix-Array jede 32. Zeile (wenig RAM).
//! - `Schnell`: Zählstand alle 64 Basen, Suffix-Array jede 8. Zeile (mehr RAM, kürzere Wege).

use std::io::{self, Read, Write};

use gc_log::Level;
use serde_json::json;

use crate::reference::{Contig, Hole, Packed, Reference};
use crate::sais;

const MAGIC: &[u8; 4] = b"GCFM";
const FORMAT: u32 = 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Variant {
    Kompakt,
    Schnell,
}

impl Variant {
    /// (Basen je Zählblock, Abstand der Suffix-Array-Stichproben)
    pub fn parameters(self) -> (u32, u32) {
        match self {
            Variant::Kompakt => (128, 32),
            Variant::Schnell => (64, 8),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
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
    /// Je Block: 4 Zählstände, dann `occ_interval / 32` Wörter BWT.
    blocks: Vec<u64>,
    sa_samples: Vec<u32>,
    reference: Reference,
}

const REP: [u64; 4] = [
    0,
    0x5555_5555_5555_5555,
    0xAAAA_AAAA_AAAA_AAAA,
    0xFFFF_FFFF_FFFF_FFFF,
];

/// Anzahl der Base `c` unter den ersten `r` Basen des Wortes (`r` höchstens 32).
#[inline]
fn count_in_word(w: u64, c: u8, r: u32) -> u64 {
    if r == 0 {
        return 0;
    }
    let x = w ^ REP[c as usize];
    let mut m = !(x | x >> 1) & 0x5555_5555_5555_5555;
    if r < 32 {
        m &= (1u64 << (2 * r)) - 1;
    }
    u64::from(m.count_ones())
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
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
        let rows = sa.len() as u64;
        let words_per_block = (occ_interval / 32) as usize;
        let stride = 4 + words_per_block;
        let n_blocks = (rows / u64::from(occ_interval)) as usize + 1;
        let mut blocks = vec![0u64; n_blocks * stride];
        let mut counts = [0u64; 4];
        let mut primary = 0u64;
        for (i, &p) in sa.iter().enumerate() {
            let within = i % occ_interval as usize;
            let base = i / occ_interval as usize * stride;
            if within == 0 {
                blocks[base..base + 4].copy_from_slice(&counts);
            }
            if p == 0 {
                primary = i as u64; // Endzeichen: als Code 0 gespeichert, beim Zählen abgezogen
                continue;
            }
            let code = text[p as usize - 1] - 1;
            counts[code as usize] += 1;
            blocks[base + 4 + within / 32] |= u64::from(code) << (2 * (within % 32));
        }
        if rows.is_multiple_of(u64::from(occ_interval)) {
            let base = (n_blocks - 1) * stride;
            blocks[base..base + 4].copy_from_slice(&counts);
        }
        let mut c_before = [0u64; 4];
        let mut sum = 1; // das Endzeichen steht in Zeile 0
        for c in 0..4 {
            c_before[c] = sum;
            sum += counts[c];
        }
        let sa_samples: Vec<u32> = sa.iter().step_by(sa_interval as usize).copied().collect();
        st.finish(json!({
            "zeilen": rows, "primaerzeile": primary,
            "a": counts[0], "c": counts[1], "g": counts[2], "t": counts[3],
            "sa_stichproben": sa_samples.len(),
        }));

        let index = FmIndex {
            rows,
            primary,
            c_before,
            occ_interval,
            sa_interval,
            blocks,
            sa_samples,
            reference,
        };

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
        (self.blocks.len() * 8 + self.sa_samples.len() * 4 + self.reference.forward.words.len() * 8)
            as u64
    }

    /// Anzahl der Base `c` in den BWT-Zeilen `0..i`.
    #[inline]
    pub fn occ(&self, c: u8, i: u64) -> u64 {
        let interval = u64::from(self.occ_interval);
        let stride = 4 + (self.occ_interval / 32) as usize;
        let base = (i / interval) as usize * stride;
        let within = (i % interval) as u32;
        let mut n = self.blocks[base + c as usize];
        let full = (within / 32) as usize;
        for w in 0..full {
            n += count_in_word(self.blocks[base + 4 + w], c, 32);
        }
        if !within.is_multiple_of(32) {
            n += count_in_word(self.blocks[base + 4 + full], c, within % 32);
        }
        // Das Endzeichen ist als Code 0 gespeichert, zählt aber nicht als A.
        let block_start = i - u64::from(within);
        if c == 0 && self.primary >= block_start && self.primary < i {
            n -= 1;
        }
        n
    }

    #[inline]
    fn bwt(&self, row: u64) -> u8 {
        let interval = u64::from(self.occ_interval);
        let stride = 4 + (self.occ_interval / 32) as usize;
        let base = (row / interval) as usize * stride;
        let within = row % interval;
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
            if row.is_multiple_of(interval) {
                return u64::from(self.sa_samples[(row / interval) as usize]) + steps;
            }
            let c = self.bwt(row);
            row = self.c_before[c as usize] + self.occ(c, row);
            steps += 1;
        }
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
        w.write_all(&(self.sa_samples.len() as u64).to_le_bytes())?;
        let mut buf = Vec::with_capacity(1 << 16);
        for chunk in self.sa_samples.chunks(1 << 14) {
            buf.clear();
            for v in chunk {
                buf.extend_from_slice(&v.to_le_bytes());
            }
            w.write_all(&buf)?;
        }
        let crc = w.crc.sum();
        w.inner.write_all(&crc.to_le_bytes())?;
        w.inner.flush()
    }

    /// Liest einen Index und prüft Format, Größen und Prüfsumme.
    pub fn read<R: Read>(input: R) -> io::Result<FmIndex> {
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
        if !(32..=1024).contains(&occ_interval)
            || occ_interval % 32 != 0
            || !(1..=1024).contains(&sa_interval)
        {
            return Err(bad("unzulässige Abstände"));
        }
        if rows == 0 || rows > sais::MAX_LEN as u64 || primary >= rows {
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
        let mut sa_samples = Vec::with_capacity(n_samples as usize);
        let mut buf = vec![0u8; 1 << 16];
        let mut left = n_samples as usize;
        while left > 0 {
            let take = left.min(buf.len() / 4);
            r.read_exact(&mut buf[..take * 4])?;
            sa_samples.extend(
                buf[..take * 4]
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .map(|b| u32::from_le_bytes(*b)),
            );
            left -= take;
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
        Ok(FmIndex {
            rows,
            primary,
            c_before,
            occ_interval: occ_interval as u32,
            sa_interval: sa_interval as u32,
            blocks,
            sa_samples,
            reference,
        })
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
        for variant in [Variant::Kompakt, Variant::Schnell] {
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
    fn occ_matches_counting_and_edges() {
        let contigs = random_genome(3);
        let index = build(&contigs, Variant::Kompakt);
        // occ gegen direktes Zählen der BWT, an jeder Zeile
        let mut counts = [0u64; 4];
        for row in 0..=index.rows() {
            for c in 0..4u8 {
                assert_eq!(index.occ(c, row), counts[c as usize], "occ({c}, {row})");
            }
            if row < index.rows() && row != index.primary {
                counts[index.bwt(row) as usize] += 1;
            }
        }
        assert_eq!(counts.iter().sum::<u64>(), index.rows() - 1);
        // leeres Muster trifft jede Zeile; ein Muster länger als das Genom trifft nie
        assert_eq!(index.count(&[]), index.rows());
        assert_eq!(index.count(&vec![0u8; 10_000]), 0);
        assert!(index.hits(&[], 10).is_empty());
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
    fn build_errors() {
        assert!(matches!(
            FmIndex::build(Reference::new(), Variant::Kompakt, 0),
            Err(BuildError::Empty)
        ));
    }
}
