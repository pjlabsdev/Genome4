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
    use std::io::Write;
    let too_long = |_| io::Error::new(io::ErrorKind::InvalidInput, "BAM-Kopf zu groß");
    let file = BufWriter::with_capacity(1 << 20, File::create(path)?);
    let mut writer = bam::io::Writer::from(bgzf::io::MultithreadedWriter::with_worker_count(
        workers(threads),
        file,
    ));
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
}
