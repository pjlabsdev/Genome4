//! FASTA lesen (auch gzip), Sequenz je Eintrag am Stück.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

use flate2::bufread::MultiGzDecoder;

/// Ein FASTA-Eintrag: Name bis zum ersten Leerzeichen, Sequenz ohne Zeilenenden.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Record {
    pub name: String,
    pub seq: Vec<u8>,
}

pub struct Reader<R> {
    inner: R,
    path: PathBuf,
    line: u64,
    pending: Option<String>,
    buf: Vec<u8>,
}

fn invalid(path: &Path, line: u64, message: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("{}: Zeile {line}: {message}", path.display()),
    )
}

/// Öffnet eine FASTA-Datei und erkennt gzip am Magic `1f 8b`.
pub fn open(path: impl AsRef<Path>) -> io::Result<Reader<Box<dyn BufRead + Send>>> {
    let path = path.as_ref().to_path_buf();
    let with_path = |e: io::Error| io::Error::new(e.kind(), format!("{}: {e}", path.display()));
    let mut file = BufReader::with_capacity(1 << 20, File::open(&path).map_err(with_path)?);
    let gz = {
        let head = file.fill_buf().map_err(with_path)?;
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
            pending: None,
            buf: Vec::new(),
        }
    }

    fn read_line(&mut self) -> io::Result<bool> {
        self.buf.clear();
        let n = self
            .inner
            .read_until(b'\n', &mut self.buf)
            .map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", self.path.display())))?;
        if n == 0 {
            return Ok(false);
        }
        self.line += 1;
        while matches!(self.buf.last(), Some(b'\n' | b'\r')) {
            self.buf.pop();
        }
        Ok(true)
    }

    fn header_name(&self) -> io::Result<String> {
        let name = self.buf[1..]
            .split(|b| b.is_ascii_whitespace())
            .next()
            .unwrap_or(&[]);
        if name.is_empty() {
            return Err(invalid(&self.path, self.line, "Kopfzeile ohne Namen"));
        }
        String::from_utf8(name.to_vec())
            .map_err(|_| invalid(&self.path, self.line, "Name ist kein gültiges UTF-8"))
    }

    /// Liest den nächsten Eintrag. `Ok(false)` am Dateiende.
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<bool> {
        record.seq.clear();
        record.name = match self.pending.take() {
            Some(name) => name,
            None => loop {
                if !self.read_line()? {
                    return Ok(false);
                }
                if self.buf.is_empty() {
                    continue;
                }
                if self.buf[0] != b'>' {
                    return Err(invalid(
                        &self.path,
                        self.line,
                        "erwartet '>' am Anfang des Eintrags",
                    ));
                }
                break self.header_name()?;
            },
        };
        while self.read_line()? {
            if self.buf.first() == Some(&b'>') {
                self.pending = Some(self.header_name()?);
                break;
            }
            if let Some(bad) = self
                .buf
                .iter()
                .find(|b| !b.is_ascii_alphabetic() && **b != b'*' && **b != b'-')
            {
                return Err(invalid(
                    &self.path,
                    self.line,
                    &format!("unerlaubtes Zeichen 0x{bad:02x} in der Sequenz"),
                ));
            }
            record.seq.extend_from_slice(&self.buf);
        }
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all(text: &[u8]) -> io::Result<Vec<Record>> {
        let mut r = Reader::new(text, "test.fa");
        let mut out = Vec::new();
        let mut rec = Record::default();
        while r.read_record(&mut rec)? {
            out.push(rec.clone());
        }
        Ok(out)
    }

    #[test]
    fn reads_multiline_and_crlf() {
        let recs = all(b">chr1 erster\r\nACGT\r\nNN\r\n\r\n>chr2\nacgt\n").unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "chr1");
        assert_eq!(recs[0].seq, b"ACGTNN");
        assert_eq!(recs[1].name, "chr2");
        assert_eq!(recs[1].seq, b"acgt");
    }

    #[test]
    fn empty_sequence_and_errors() {
        let recs = all(b">leer\n>voll\nA\n").unwrap();
        assert_eq!(recs[0].seq, b"");
        assert_eq!(recs[1].seq, b"A");
        let err = all(b"ACGT\n").unwrap_err().to_string();
        assert!(err.contains("Zeile 1"), "{err}");
        let err = all(b">a\nAC GT\n").unwrap_err().to_string();
        assert!(err.contains("Zeile 2") && err.contains("0x20"), "{err}");
    }
}
