//! gc-io: Lesen und Schreiben der Genomdaten für Genome 4.0.
//!
//! - [`fastq`]: FASTQ und FASTQ.gz lesen, gepaarte Dateien prüfen.
//! - [`readname`]: Herkunft einer Lesung (Gerät, Flowcell, Lane, Kachel) aus dem Namen.
//! - [`bam`]: BAM mit paralleler BGZF-Kompression lesen und schreiben.
//! - [`trim`]: Read-Paare trimmen wie fastp 1.4.0.

pub mod bam;
pub mod fasta;
pub mod fastq;
mod known_adapters;
pub mod readname;
pub mod trim;
