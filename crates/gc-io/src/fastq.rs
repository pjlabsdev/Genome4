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
    let (gz, bgzf) = {
        let head = file.fill_buf().map_err(io_err)?;
        let gz = head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b;
        // BGZF: gzip mit Zusatzfeld "BC" (SAM-Spezifikation, Abschnitt 4.1)
        let bgzf = gz
            && head.len() >= 16
            && head[3] & 4 != 0
            && head[12] == b'B'
            && head[13] == b'C'
            && head[14] == 2
            && head[15] == 0;
        (gz, bgzf)
    };
    let inner: Box<dyn BufRead + Send> = if gz && threaded_unpacking() {
        if bgzf {
            // Blöcke unabhängig: auf mehreren Kernen entpacken
            let workers = std::thread::available_parallelism()
                .map_or(4, |n| n.get())
                .clamp(1, 4);
            Box::new(noodles::bgzf::io::MultithreadedReader::with_worker_count(
                std::num::NonZero::new(workers).expect("mindestens 1"),
                file,
            ))
        } else {
            // Entpacken in eigenem Faden, das Zerlegen läuft daneben
            Box::new(ThreadedRead::new(MultiGzDecoder::new(file)))
        }
    } else if gz {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(file)
    };
    Ok(Reader::new(inner, path))
}

/// Entpacken in eigenen Fäden (Schalter zurück: `GC_LESEN_ENTPACKEN=0`).
fn threaded_unpacking() -> bool {
    std::env::var("GC_LESEN_ENTPACKEN").map_or(true, |v| v != "0")
}

/// Liest eine Quelle in einem eigenen Faden in Stücken von 1 MB voraus. Ein Lesefehler kommt an
/// derselben Stelle im Strom an wie ohne Faden.
pub struct ThreadedRead {
    rx: std::sync::mpsc::Receiver<io::Result<Vec<u8>>>,
    buf: Vec<u8>,
    pos: usize,
    done: bool,
}

impl ThreadedRead {
    pub fn new<R: Read + Send + 'static>(mut inner: R) -> Self {
        let (tx, rx) = std::sync::mpsc::sync_channel(4);
        std::thread::spawn(move || {
            loop {
                let mut chunk = vec![0u8; 1 << 20];
                let mut n = 0;
                while n < chunk.len() {
                    match inner.read(&mut chunk[n..]) {
                        Ok(0) => break,
                        Ok(k) => n += k,
                        Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                        Err(e) => {
                            if n > 0 {
                                chunk.truncate(n);
                                let _ = tx.send(Ok(chunk));
                            }
                            let _ = tx.send(Err(e));
                            return;
                        }
                    }
                }
                chunk.truncate(n);
                let last = n == 0;
                if tx.send(Ok(chunk)).is_err() || last {
                    return;
                }
            }
        });
        Self {
            rx,
            buf: Vec::new(),
            pos: 0,
            done: false,
        }
    }
}

impl Read for ThreadedRead {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let avail = self.fill_buf()?;
        let n = avail.len().min(out.len());
        out[..n].copy_from_slice(&avail[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for ThreadedRead {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.pos >= self.buf.len() && !self.done {
            match self.rx.recv() {
                Ok(Ok(chunk)) => {
                    self.done = chunk.is_empty();
                    self.buf = chunk;
                    self.pos = 0;
                }
                Ok(Err(e)) => {
                    self.done = true;
                    return Err(e);
                }
                Err(_) => {
                    self.done = true;
                    return Err(io::Error::other("Lesefaden vorzeitig beendet"));
                }
            }
        }
        Ok(&self.buf[self.pos.min(self.buf.len())..])
    }

    fn consume(&mut self, n: usize) {
        self.pos += n;
    }
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

    /// Alle Reads einer Datei, über [`open`].
    fn all_records(p: &Path) -> Vec<Record> {
        let mut r = open(p).unwrap();
        let mut out = Vec::new();
        loop {
            let mut rec = Record::default();
            if !r.read_record(&mut rec).unwrap() {
                return out;
            }
            out.push(rec);
        }
    }

    #[test]
    fn gzip_and_bgzf_give_the_same_reads_as_plain_text() {
        use flate2::{Compression, write::GzEncoder};
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("gcio-bgzf-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        // gut 3 MB: mehrere Stücke des Lesefadens und viele BGZF-Blöcke
        let mut text = Vec::new();
        let mut x = 7u64;
        for i in 0..20_000 {
            let seq: String = (0..120)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"ACGT"[(x & 3) as usize] as char
                })
                .collect();
            writeln!(text, "@r{i} c\n{seq}\n+\n{}", "I".repeat(120)).unwrap();
        }
        let plain = dir.join("a.fq");
        std::fs::write(&plain, &text).unwrap();
        let gz = dir.join("a.fq.gz");
        let mut enc = GzEncoder::new(File::create(&gz).unwrap(), Compression::fast());
        enc.write_all(&text).unwrap();
        enc.finish().unwrap();
        let bg = dir.join("b.fq.gz");
        let mut w = noodles::bgzf::io::Writer::new(File::create(&bg).unwrap());
        w.write_all(&text).unwrap();
        w.finish().unwrap();
        let want = all_records(&plain);
        assert_eq!(want.len(), 20_000);
        assert_eq!(all_records(&gz), want);
        assert_eq!(all_records(&bg), want);
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn threaded_read_passes_bytes_and_errors_in_order() {
        struct Failing(usize);
        impl Read for Failing {
            fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
                if self.0 == 0 {
                    return Err(io::Error::other("kaputt"));
                }
                let n = out.len().min(self.0).min(1000);
                out[..n].fill(b'x');
                self.0 -= n;
                Ok(n)
            }
        }
        let mut r = ThreadedRead::new(Failing(2_500_000));
        let mut got = Vec::new();
        let err = r.read_to_end(&mut got).unwrap_err();
        assert_eq!(err.to_string(), "kaputt");
        assert_eq!(got.len(), 2_500_000);
        let mut ok = ThreadedRead::new(&b"abc"[..]);
        let mut s = String::new();
        ok.read_to_string(&mut s).unwrap();
        assert_eq!(s, "abc");
    }
}
