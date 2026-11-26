//! BAM lesen und schreiben mit paralleler BGZF-(De-)Kompression über alle Kerne.

use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use std::num::NonZero;
use std::path::Path;

use noodles::bam;
use noodles::bgzf;
use noodles::sam;

pub type MtReader = bam::io::Reader<bgzf::io::MultithreadedReader<BufReader<File>>>;
pub type MtWriter = bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<File>>>;

fn workers(threads: usize) -> NonZero<usize> {
    NonZero::new(threads)
        .unwrap_or_else(|| std::thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN))
}

/// Öffnet eine BAM zum Lesen; `threads == 0` heißt alle Kerne.
pub fn open(path: impl AsRef<Path>, threads: usize) -> io::Result<(MtReader, sam::Header)> {
    let file = BufReader::with_capacity(1 << 20, File::open(path)?);
    let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
        workers(threads),
        file,
    ));
    let header = reader.read_header()?;
    Ok((reader, header))
}

/// Legt eine BAM zum Schreiben an und schreibt den Kopf.
pub fn create(
    path: impl AsRef<Path>,
    header: &sam::Header,
    threads: usize,
) -> io::Result<MtWriter> {
    let file = BufWriter::with_capacity(1 << 20, File::create(path)?);
    let mut writer = bam::io::Writer::from(bgzf::io::MultithreadedWriter::with_worker_count(
        workers(threads),
        file,
    ));
    writer.write_header(header)?;
    Ok(writer)
}

/// Schließt eine BAM sauber ab (letzte Blöcke und BGZF-Endmarke).
pub fn finish(mut writer: MtWriter) -> io::Result<()> {
    use std::io::Write;
    let mut inner = writer.get_mut().finish()?;
    inner.flush()
}

/// Öffnet eine BAM und liefert zusätzlich den Kopftext Byte für Byte wie in der Datei.
/// Nötig, wenn eine BAM unverändert weitergegeben wird: die Reihenfolge der Kopfzeilen
/// (z. B. die `@PG`-Kette) bleibt nur über den Rohtext erhalten.
pub fn open_raw(
    path: impl AsRef<Path>,
    threads: usize,
) -> io::Result<(MtReader, sam::Header, Vec<u8>)> {
    use std::io::Read;
    let file = BufReader::with_capacity(1 << 20, File::open(path)?);
    let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
        workers(threads),
        file,
    ));
    let mut hr = reader.header_reader();
    hr.read_magic_number()?;
    let mut raw = Vec::new();
    {
        let mut text = hr.raw_sam_header_reader()?;
        text.read_to_end(&mut raw)?;
        text.discard_to_end()?;
    }
    let references = hr.read_reference_sequences()?;
    // Der Text darf mit NUL aufgefüllt sein; für das Parsen zählt nur der Teil davor.
    let end = raw.iter().position(|&b| b == 0).unwrap_or(raw.len());
    let text = std::str::from_utf8(&raw[..end])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut header: sam::Header = text
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    if header.reference_sequences().is_empty() {
        *header.reference_sequences_mut() = references;
    } else if header.reference_sequences().len() != references.len()
        || header.reference_sequences().keys().ne(references.keys())
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Referenzliste im Kopftext passt nicht zur binären Referenzliste",
        ));
    }
    Ok((reader, header, raw))
}

/// Legt eine BAM an und schreibt den Kopftext unverändert (`raw_text`), dazu die
/// binäre Referenzliste aus `header`.
pub fn create_raw(
    path: impl AsRef<Path>,
    header: &sam::Header,
    raw_text: &[u8],
    threads: usize,
) -> io::Result<MtWriter> {
    create_raw_level(path, header, raw_text, threads, DEFAULT_LEVEL)
}

/// Stufe der Deflate-Kompression, wenn nichts anderes verlangt ist (wie samtools und noodles).
pub const DEFAULT_LEVEL: u8 = 6;

/// Wie [`create_raw`], mit Kompressionsstufe 0 (keine) bis 9 (kleinste Datei). Die Stufe ändert
/// nur die Bytes der Datei, nicht die Sätze.
pub fn create_raw_level(
    path: impl AsRef<Path>,
    header: &sam::Header,
    raw_text: &[u8],
    threads: usize,
    level: u8,
) -> io::Result<MtWriter> {
    use std::io::Write;
    let too_long = |_| io::Error::new(io::ErrorKind::InvalidInput, "BAM-Kopf zu groß");
    let level = bgzf::io::writer::CompressionLevel::new(level).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "Kompressionsstufe muss 0 bis 9 sein",
        )
    })?;
    let file = BufWriter::with_capacity(1 << 20, File::create(path)?);
    let mut writer = bam::io::Writer::from(
        bgzf::io::multithreaded_writer::Builder::default()
            .set_worker_count(workers(threads))
            .set_compression_level(level)
            .build_from_writer(file),
    );
    let w = writer.get_mut();
    w.write_all(b"BAM\x01")?;
    w.write_all(
        &u32::try_from(raw_text.len())
            .map_err(too_long)?
            .to_le_bytes(),
    )?;
    w.write_all(raw_text)?;
    let refs = header.reference_sequences();
    w.write_all(&u32::try_from(refs.len()).map_err(too_long)?.to_le_bytes())?;
    for (name, map) in refs {
        w.write_all(
            &u32::try_from(name.len() + 1)
                .map_err(too_long)?
                .to_le_bytes(),
        )?;
        w.write_all(name)?;
        w.write_all(&[0])?;
        let len = u32::try_from(usize::from(map.length())).map_err(too_long)?;
        w.write_all(&len.to_le_bytes())?;
    }
    Ok(writer)
}

/// Namen und Längen der Contigs aus dem Kopf, in der Reihenfolge der Contig-Nummern.
pub fn contigs(header: &sam::Header) -> Vec<(String, u64)> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, map)| (name.to_string(), usize::from(map.length()) as u64))
        .collect()
}

/// Legt eine BAM an; der Kopf kommt als SAM-Text (`@SQ`-Zeilen usw.) und wird unverändert
/// geschrieben.
pub fn create_from_text(
    path: impl AsRef<Path>,
    header_text: &str,
    threads: usize,
) -> io::Result<MtWriter> {
    create_from_text_level(path, header_text, threads, DEFAULT_LEVEL)
}

/// Wie [`create_from_text`], mit Kompressionsstufe 0 bis 9.
pub fn create_from_text_level(
    path: impl AsRef<Path>,
    header_text: &str,
    threads: usize,
    level: u8,
) -> io::Result<MtWriter> {
    let header: sam::Header = header_text
        .parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    create_raw_level(path, &header, header_text.as_bytes(), threads, level)
}

/// Schreibt fertig kodierte BAM-Sätze (jeder mit seinem Längenfeld vorn) unverändert.
pub fn write_encoded(writer: &mut MtWriter, records: &[u8]) -> io::Result<()> {
    use std::io::Write;
    writer.get_mut().write_all(records)
}

/// Liest den nächsten Datensatz roh (mit seinem Längenfeld) und hängt ihn an `out` an.
/// `false` am Ende der Datei. Der Kopf muss schon gelesen sein.
pub fn read_encoded(reader: &mut MtReader, out: &mut Vec<u8>) -> io::Result<bool> {
    use std::io::Read;
    let inner = reader.get_mut();
    let mut len = [0u8; 4];
    let mut got = 0;
    while got < 4 {
        let n = inner.read(&mut len[got..])?;
        if n == 0 {
            break;
        }
        got += n;
    }
    if got == 0 {
        return Ok(false);
    }
    let cut = || io::Error::new(io::ErrorKind::UnexpectedEof, "BAM-Datensatz abgeschnitten");
    if got < 4 {
        return Err(cut());
    }
    let n = u32::from_le_bytes(len) as usize;
    let start = out.len();
    out.extend_from_slice(&len);
    out.resize(start + 4 + n, 0);
    inner.read_exact(&mut out[start + 4..]).map_err(|_| cut())?;
    Ok(true)
}

/// Virtuelle Position der Lesestelle: Dateistelle des BGZF-Blocks mal 65.536 plus Stelle im
/// entpackten Block. So verweist ein Index auf Datensätze.
pub fn virtual_position(reader: &MtReader) -> u64 {
    u64::from(reader.get_ref().virtual_position())
}

/// Liest eine BAM und schreibt Kopf und alle Datensätze unverändert neu.
/// Grundlage des Rundreise-Tests gegen samtools (M1).
pub fn roundtrip(src: impl AsRef<Path>, dst: impl AsRef<Path>, threads: usize) -> io::Result<u64> {
    let (mut reader, header, raw) = open_raw(src, threads)?;
    let mut writer = create_raw(dst, &header, &raw, threads)?;
    let mut record = bam::Record::default();
    let mut n = 0u64;
    while reader.read_record(&mut record)? != 0 {
        writer.write_record(&header, &record)?;
        n += 1;
    }
    finish(writer)?;
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::io::Write as _;

    #[test]
    fn roundtrip_keeps_records() {
        let dir = std::env::temp_dir().join(format!("gcio-bam-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (a, b) = (dir.join("a.bam"), dir.join("b.bam"));

        let header = sam::Header::default();
        {
            let mut w = create(&a, &header, 2).unwrap();
            for i in 0..5000 {
                let rec = RecordBuf::builder()
                    .set_name(format!("r{i}"))
                    .set_sequence(b"ACGTACGT".to_vec().into())
                    .set_quality_scores(vec![30u8; 8].into())
                    .build();
                w.write_alignment_record(&header, &rec).unwrap();
            }
            finish(w).unwrap();
        }
        assert_eq!(roundtrip(&a, &b, 3).unwrap(), 5000);

        let read_all = |p: &Path| {
            let (mut r, h) = open(p, 1).unwrap();
            r.record_bufs(&h).map(|x| x.unwrap()).collect::<Vec<_>>()
        };
        assert_eq!(read_all(&a), read_all(&b));
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn roundtrip_keeps_header_line_order() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("gcio-bamhdr-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (a, b) = (dir.join("a.bam"), dir.join("b.bam"));

        // @PG vor @SQ: eine Reihenfolge, die beim Neuaufbau des Kopfes verloren ginge.
        let text =
            b"@HD\tVN:1.6\tSO:coordinate\n@PG\tID:x\tPN:x\n@SQ\tSN:chr1\tLN:1000\n@CO\thallo\n";
        let header: sam::Header = std::str::from_utf8(text).unwrap().parse().unwrap();
        let mut w = create_raw(&a, &header, text, 1).unwrap();
        w.get_mut().flush().unwrap();
        finish(w).unwrap();

        assert_eq!(roundtrip(&a, &b, 2).unwrap(), 0);
        let (_, parsed, raw) = open_raw(&b, 1).unwrap();
        assert_eq!(raw, text);
        assert_eq!(parsed.reference_sequences().len(), 1);
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn level_changes_bytes_not_records() {
        let dir = std::env::temp_dir().join(format!("gcio-bamstufe-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";
        let header: sam::Header = text.parse().unwrap();
        let write = |p: &Path, level: Option<u8>| {
            let mut w = match level {
                None => create_from_text(p, text, 2).unwrap(),
                Some(l) => create_from_text_level(p, text, 2, l).unwrap(),
            };
            for i in 0..20_000u32 {
                let rec = RecordBuf::builder()
                    .set_name(format!("r{i}"))
                    .set_sequence(b"ACGTACGTTTGACCA".to_vec().into())
                    .set_quality_scores(vec![(i % 40) as u8; 15].into())
                    .build();
                w.write_alignment_record(&header, &rec).unwrap();
            }
            finish(w).unwrap();
        };
        let (d, l6, l1) = (dir.join("d.bam"), dir.join("6.bam"), dir.join("1.bam"));
        write(&d, None);
        write(&l6, Some(DEFAULT_LEVEL));
        write(&l1, Some(1));
        // Vorgabe und Stufe 6 sind dieselbe Datei
        assert_eq!(std::fs::read(&d).unwrap(), std::fs::read(&l6).unwrap());
        assert_ne!(std::fs::read(&d).unwrap(), std::fs::read(&l1).unwrap());
        let read_all = |p: &Path| {
            let (mut r, h) = open(p, 1).unwrap();
            r.record_bufs(&h).map(|x| x.unwrap()).collect::<Vec<_>>()
        };
        assert_eq!(read_all(&d), read_all(&l1));
        assert!(create_from_text_level(dir.join("x.bam"), text, 1, 10).is_err());
        std::fs::remove_dir_all(dir).ok();
    }
}
