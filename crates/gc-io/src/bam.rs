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

/// Liest eine BAM und schreibt alle Datensätze unverändert neu.
/// Grundlage des Rundreise-Tests gegen samtools (M1).
pub fn roundtrip(src: impl AsRef<Path>, dst: impl AsRef<Path>, threads: usize) -> io::Result<u64> {
    let (mut reader, header) = open(src, threads)?;
    let mut writer = create(dst, &header, threads)?;
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
}
