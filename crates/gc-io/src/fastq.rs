//! FASTQ lesen. Erkennt gzip/BGZF am Dateianfang, prüft jeden Datensatz streng
//! und meldet Fehler mit Datei und Zeilennummer, statt nach Stunden abzubrechen.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

/// Ein FASTQ-Datensatz. `name` ist der Name ohne `@` und ohne Kommentar.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Record {
    pub name: Vec<u8>,
    pub comment: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

impl Record {
    /// Name ohne Paarkennung (`/1`, `/2`), damit beide Lesungen eines Paars vergleichbar sind.
    pub fn pair_key(&self) -> &[u8] {
        strip_mate_suffix(&self.name)
    }
}

/// Entfernt `/1` oder `/2` am Namensende.
pub fn strip_mate_suffix(name: &[u8]) -> &[u8] {
    match name {
        [rest @ .., b'/', b'1' | b'2'] => rest,
        _ => name,
    }
}

#[derive(Debug, thiserror::Error)]
pub enum FastqError {
    #[error("{path}: {source}")]
    Io { path: PathBuf, source: io::Error },
    #[error("{path}, Zeile {line}: {message}")]
    Format {
        path: PathBuf,
        line: u64,
        message: String,
    },
    #[error("Paar {index}: Namen passen nicht zusammen ({left} / {right})")]
    PairMismatch {
        index: u64,
        left: String,
        right: String,
    },
    #[error("{longer} hat mehr Lesungen als {shorter} (abgeschnittene oder falsche Datei?)")]
    PairLength { longer: PathBuf, shorter: PathBuf },
}

/// Liest FASTQ-Datensätze aus einer Datei (roh, gzip oder BGZF).
pub struct Reader<R> {
    inner: R,
    path: PathBuf,
    line: u64,
    buf: Vec<u8>,
}

/// Öffnet eine FASTQ-Datei und erkennt gzip am Magic `1f 8b`.
pub fn open(path: impl AsRef<Path>) -> Result<Reader<Box<dyn BufRead + Send>>, FastqError> {
    let path = path.as_ref().to_path_buf();
    let io_err = |source| FastqError::Io {
        path: path.clone(),
        source,
    };
    let mut file = BufReader::with_capacity(1 << 20, File::open(&path).map_err(io_err)?);
    let gz = {
        let head = file.fill_buf().map_err(io_err)?;
        head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
    };
    let inner: Box<dyn BufRead + Send> = if gz {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(file)
    };
    Ok(Reader::new(inner, path))
}

impl<R: BufRead> Reader<R> {
    pub fn new(inner: R, path: impl Into<PathBuf>) -> Self {
        Self {
            inner,
            path: path.into(),
            line: 0,
            buf: Vec::with_capacity(512),
        }
    }

    fn format(&self, message: impl Into<String>) -> FastqError {
        FastqError::Format {
            path: self.path.clone(),
            line: self.line,
            message: message.into(),
        }
    }

    fn next_line(&mut self) -> Result<bool, FastqError> {
        self.buf.clear();
        let n = self
            .inner
            .read_until(b'\n', &mut self.buf)
            .map_err(|source| FastqError::Io {
                path: self.path.clone(),
                source,
            })?;
        if n == 0 {
            return Ok(false);
        }
        self.line += 1;
        if self.buf.last() == Some(&b'\n') {
            self.buf.pop();
            if self.buf.last() == Some(&b'\r') {
                self.buf.pop();
            }
        }
        Ok(true)
    }

    /// Liest den nächsten Datensatz in `record`. `Ok(false)` am Dateiende.
    pub fn read_record(&mut self, record: &mut Record) -> Result<bool, FastqError> {
        // Leerzeilen am Dateiende sind erlaubt, mitten im Datensatz nicht.
        loop {
            if !self.next_line()? {
                return Ok(false);
            }
            if !self.buf.is_empty() {
                break;
            }
        }
        let header = &self.buf;
        if header.first() != Some(&b'@') {
            return Err(self.format("Kopfzeile beginnt nicht mit '@'"));
        }
        let body = &header[1..];
        let split = body.iter().position(|&b| b == b' ' || b == b'\t');
        record.name.clear();
        record.comment.clear();
        match split {
            Some(i) => {
                record.name.extend_from_slice(&body[..i]);
                record.comment.extend_from_slice(&body[i + 1..]);
            }
            None => record.name.extend_from_slice(body),
        }
        if record.name.is_empty() {
            return Err(self.format("leerer Lesungsname"));
        }

        if !self.next_line()? {
            return Err(self.format("Datei endet nach der Kopfzeile (abgeschnitten?)"));
        }
        record.seq.clear();
        record.seq.extend_from_slice(&self.buf);
        for b in record.seq.iter_mut() {
            *b = b.to_ascii_uppercase();
            if !matches!(*b, b'A' | b'C' | b'G' | b'T' | b'N') {
                if *b == b'.' {
                    *b = b'N';
                } else {
                    let c = *b as char;
                    return Err(self.format(format!("unerlaubtes Zeichen '{c}' in der Sequenz")));
                }
            }
        }

        if !self.next_line()? || self.buf.first() != Some(&b'+') {
            return Err(self.format("Trennzeile '+' fehlt"));
        }
        if !self.next_line()? {
            return Err(self.format("Qualitätszeile fehlt (abgeschnitten?)"));
        }
        if self.buf.len() != record.seq.len() {
            let (q, s) = (self.buf.len(), record.seq.len());
            return Err(self.format(format!("Qualität hat {q} Zeichen, Sequenz {s}")));
        }
        if let Some(&bad) = self.buf.iter().find(|&&q| !(33..=126).contains(&q)) {
            return Err(self.format(format!("ungültiger Qualitätswert {bad}")));
        }
        record.qual.clear();
        record.qual.extend_from_slice(&self.buf);
        Ok(true)
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
}

/// Liest zwei FASTQ-Dateien im Gleichschritt und prüft, dass die Namen zusammenpassen.
pub struct PairedReader<R1, R2> {
    r1: Reader<R1>,
    r2: Reader<R2>,
    index: u64,
}

impl<R1: BufRead, R2: BufRead> PairedReader<R1, R2> {
    pub fn new(r1: Reader<R1>, r2: Reader<R2>) -> Self {
        Self { r1, r2, index: 0 }
    }

    pub fn read_pair(&mut self, a: &mut Record, b: &mut Record) -> Result<bool, FastqError> {
        let ha = self.r1.read_record(a)?;
        let hb = self.r2.read_record(b)?;
        match (ha, hb) {
            (false, false) => Ok(false),
            (true, false) => Err(FastqError::PairLength {
                longer: self.r1.path().into(),
                shorter: self.r2.path().into(),
            }),
            (false, true) => Err(FastqError::PairLength {
                longer: self.r2.path().into(),
                shorter: self.r1.path().into(),
            }),
            (true, true) => {
                self.index += 1;
                if a.pair_key() != b.pair_key() {
                    return Err(FastqError::PairMismatch {
                        index: self.index,
                        left: String::from_utf8_lossy(&a.name).into_owned(),
                        right: String::from_utf8_lossy(&b.name).into_owned(),
                    });
                }
                Ok(true)
            }
        }
    }
}

/// Kennzahlen eines Prüflaufs über eine FASTQ-Datei oder ein Paar.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Stats {
    pub reads: u64,
    pub bases: u64,
    pub min_len: usize,
    pub max_len: usize,
    pub q30_bases: u64,
    pub n_bases: u64,
    pub gc_bases: u64,
}

impl Stats {
    pub fn add(&mut self, r: &Record) {
        let len = r.seq.len();
        if self.reads == 0 || len < self.min_len {
            self.min_len = len;
        }
        self.max_len = self.max_len.max(len);
        self.reads += 1;
        self.bases += len as u64;
        self.q30_bases += r.qual.iter().filter(|&&q| q >= 33 + 30).count() as u64;
        for &b in &r.seq {
            match b {
                b'G' | b'C' => self.gc_bases += 1,
                b'N' => self.n_bases += 1,
                _ => {}
            }
        }
    }
}

/// Erlaubt `Reader` direkt über einen beliebigen `Read` (z. B. Tests).
pub fn from_read<R: Read>(r: R, path: impl Into<PathBuf>) -> Reader<BufReader<R>> {
    Reader::new(BufReader::new(r), path)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reader(s: &str) -> Reader<&[u8]> {
        Reader::new(s.as_bytes(), "test.fq")
    }

    #[test]
    fn reads_records_and_comment() {
        let mut r = reader("@r1/1 BC:Z:ACGT\nacgtn\n+\nIIIII\n@r2\nACG\n+r2\n###\n\n");
        let mut rec = Record::default();
        assert!(r.read_record(&mut rec).unwrap());
        assert_eq!(rec.name, b"r1/1");
        assert_eq!(rec.pair_key(), b"r1");
        assert_eq!(rec.comment, b"BC:Z:ACGT");
        assert_eq!(rec.seq, b"ACGTN");
        assert!(r.read_record(&mut rec).unwrap());
        assert_eq!(rec.name, b"r2");
        assert!(!r.read_record(&mut rec).unwrap());
    }

    #[test]
    fn rejects_truncated_and_mismatched_quality() {
        let mut rec = Record::default();
        let e = reader("@r1\nACGT\n+\nIII\n")
            .read_record(&mut rec)
            .unwrap_err();
        assert!(e.to_string().contains("Zeile 4"), "{e}");
        let e = reader("@r1\nACGT\n").read_record(&mut rec).unwrap_err();
        assert!(e.to_string().contains("'+'"), "{e}");
        let e = reader("r1\nACGT\n+\nIIII\n")
            .read_record(&mut rec)
            .unwrap_err();
        assert!(e.to_string().contains("'@'"), "{e}");
        let e = reader("@r1\nACXT\n+\nIIII\n")
            .read_record(&mut rec)
            .unwrap_err();
        assert!(e.to_string().contains("'X'"), "{e}");
    }

    #[test]
    fn paired_reader_checks_names_and_lengths() {
        let (mut a, mut b) = (Record::default(), Record::default());
        let mut p = PairedReader::new(
            reader("@x/1\nA\n+\nI\n@y/1\nA\n+\nI\n"),
            reader("@x/2\nA\n+\nI\n@z/2\nA\n+\nI\n"),
        );
        assert!(p.read_pair(&mut a, &mut b).unwrap());
        assert!(matches!(
            p.read_pair(&mut a, &mut b),
            Err(FastqError::PairMismatch { index: 2, .. })
        ));

        let mut p = PairedReader::new(
            reader("@x\nA\n+\nI\n@y\nA\n+\nI\n"),
            reader("@x\nA\n+\nI\n"),
        );
        assert!(p.read_pair(&mut a, &mut b).unwrap());
        assert!(matches!(
            p.read_pair(&mut a, &mut b),
            Err(FastqError::PairLength { .. })
        ));
    }

    #[test]
    fn stats_count_q30_gc_and_n() {
        let mut s = Stats::default();
        s.add(&Record {
            seq: b"GCAN".to_vec(),
            qual: b"?#?#".to_vec(),
            ..Default::default()
        });
        assert_eq!(
            (s.reads, s.bases, s.q30_bases, s.gc_bases, s.n_bases),
            (1, 4, 2, 2, 1)
        );
    }

    #[test]
    fn opens_gzip_by_magic() {
        use flate2::{Compression, write::GzEncoder};
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("gcio-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("a.fq.gz");
        let mut enc = GzEncoder::new(File::create(&p).unwrap(), Compression::fast());
        enc.write_all(b"@r\nAC\n+\nII\n").unwrap();
        enc.finish().unwrap();
        let mut r = open(&p).unwrap();
        let mut rec = Record::default();
        assert!(r.read_record(&mut rec).unwrap());
        assert_eq!(rec.seq, b"AC");
        std::fs::remove_dir_all(dir).ok();
    }
}
