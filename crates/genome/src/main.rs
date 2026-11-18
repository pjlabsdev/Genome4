//! `genome`: ein Programm für Oberfläche, Arbeitsprozesse und Werkzeuge (Genome 4.0).
//! Heute: erste Werkzeuge aus Meilenstein M1 des FASTQ→BAM-Schritts.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use gc_io::fastq::{self, PairedReader, Record, Stats};
use gc_io::readname;
use gc_schema::Envelope;
use serde_json::json;

#[derive(Parser)]
#[command(name = "genome", version, about = "Genome 4.0")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// FASTQ (einzeln oder gepaart) vollständig prüfen und Kennzahlen als JSON ausgeben.
    FastqCheck { r1: PathBuf, r2: Option<PathBuf> },
    /// BAM lesen und unverändert neu schreiben (Rundreise-Test gegen samtools).
    BamRoundtrip {
        src: PathBuf,
        dst: PathBuf,
        /// Kerne für BGZF (0 = alle)
        #[arg(short, long, default_value_t = 0)]
        threads: usize,
    },
}

fn stats_json(s: &Stats) -> serde_json::Value {
    let pct = |x: u64| {
        if s.bases == 0 {
            0.0
        } else {
            100.0 * x as f64 / s.bases as f64
        }
    };
    json!({
        "reads": s.reads, "bases": s.bases,
        "min_len": s.min_len, "max_len": s.max_len,
        "q30_pct": pct(s.q30_bases), "gc_pct": pct(s.gc_bases), "n_pct": pct(s.n_bases),
    })
}

fn origin_json(r: &Record) -> serde_json::Value {
    match readname::parse(&r.name) {
        Some(o) => json!({
            "platform": format!("{:?}", o.platform), "flowcell": o.flowcell, "lane": o.lane,
        }),
        None => serde_json::Value::Null,
    }
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::FastqCheck { r1, r2 } => {
            let started = std::time::Instant::now();
            let mut a = Record::default();
            let mut sa = Stats::default();
            let mut origin = serde_json::Value::Null;
            let result = match r2 {
                None => {
                    let mut r = fastq::open(&r1)?;
                    while r.read_record(&mut a)? {
                        if sa.reads == 0 {
                            origin = origin_json(&a);
                        }
                        sa.add(&a);
                    }
                    json!({ "paired": false, "r1": stats_json(&sa), "origin": origin })
                }
                Some(r2) => {
                    let mut b = Record::default();
                    let mut sb = Stats::default();
                    let mut p = PairedReader::new(fastq::open(&r1)?, fastq::open(&r2)?);
                    while p.read_pair(&mut a, &mut b)? {
                        if sa.reads == 0 {
                            origin = origin_json(&a);
                        }
                        sa.add(&a);
                        sb.add(&b);
                    }
                    json!({ "paired": true, "pairs": sa.reads, "r1": stats_json(&sa), "r2": stats_json(&sb), "origin": origin })
                }
            };
            let mut result = result;
            result["seconds"] = json!(started.elapsed().as_secs_f64());
            println!(
                "{}",
                Envelope::new("genome.fastq-check.v1", result).to_json()?
            );
        }
        Command::BamRoundtrip { src, dst, threads } => {
            let started = std::time::Instant::now();
            let n = gc_io::bam::roundtrip(&src, &dst, threads)
                .with_context(|| format!("{} → {}", src.display(), dst.display()))?;
            let result = json!({ "records": n, "seconds": started.elapsed().as_secs_f64() });
            println!(
                "{}",
                Envelope::new("genome.bam-roundtrip.v1", result).to_json()?
            );
        }
    }
    Ok(())
}
