//! gc-io: Lesen und Schreiben der Genomdaten für Genome 4.0.
//!
//! - [`fastq`]: FASTQ und FASTQ.gz lesen, gepaarte Dateien prüfen.
//! - [`readname`]: Herkunft einer Lesung (Gerät, Flowcell, Lane, Kachel) aus dem Namen.
//! - [`bam`]: BAM mit paralleler BGZF-Kompression lesen und schreiben.

pub mod bam;
pub mod fasta;
pub mod fastq;
pub mod readname;
