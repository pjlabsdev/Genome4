//! genome-validate: misst Genome-4-Werkzeuge gegen eine bekannte Wahrheit.
//!
//! `reads-erzeugen` zieht Read-Paare aus einer Probe mit bekannten Varianten, der wahre Ort
//! steht im Read-Namen. `ausrichtung-pruefen` liest eine oder mehrere BAMs mit diesen Reads und
//! zählt je MAPQ, wie viele am wahren Ort liegen. Alle Zahlen in Messberichten und PRs zur
//! Genauigkeit der Ausrichtung kommen aus diesem Programm.

mod eval;
mod rng;
mod sim;

use std::io::Write;
use std::path::PathBuf;

use anyhow::{Result, bail};
use clap::{Parser, Subcommand, ValueEnum};
use serde_json::json;

#[derive(Parser)]
#[command(
    name = "genome-validate",
    version,
    about = "Genome 4: Prüfung gegen die Wahrheit"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Clone, Copy, ValueEnum)]
enum ProfileArg {
    Illumina,
    Mgi,
}

#[derive(Subcommand)]
enum Command {
    /// Read-Paare mit bekannter Herkunft erzeugen.
    ReadsErzeugen {
        /// Referenz als FASTA (auch gzip).
        #[arg(long)]
        referenz: PathBuf,
        /// Anfang der Ausgabedateien: <ausgabe>_1.fq.gz, _2.fq.gz, .contigs.tsv,
        /// .varianten.tsv, .json.
        #[arg(long)]
        ausgabe: PathBuf,
        #[arg(long, default_value_t = 1_000_000)]
        paare: u64,
        #[arg(long, default_value_t = 150)]
        laenge: usize,
        /// Mittlere Länge der Stücke.
        #[arg(long, default_value_t = 400.0)]
        insert: f64,
        #[arg(long, default_value_t = 60.0)]
        streuung: f64,
        #[arg(long, default_value_t = 1)]
        saat: u64,
        #[arg(long, value_enum, default_value = "illumina")]
        profil: ProfileArg,
        /// SNVs je Base der Probe (Mensch: etwa 1 auf 1.000).
        #[arg(long, default_value_t = 0.001)]
        snv_rate: f64,
        /// Kleine Indels je Base der Probe (Mensch: etwa 1 auf 8.000).
        #[arg(long, default_value_t = 0.000_125)]
        indel_rate: f64,
        #[arg(long, default_value_t = 10)]
        max_indel: u64,
    },
    /// BAMs gegen den wahren Ort im Read-Namen prüfen.
    AusrichtungPruefen {
        /// Contig-Liste der Simulation (<ausgabe>.contigs.tsv).
        #[arg(long)]
        contigs: PathBuf,
        /// name=pfad, mehrfach; die erste BAM wird mit jeder weiteren Read für Read verglichen.
        #[arg(long = "bam", required = true, value_parser = eval::parse_named)]
        bams: Vec<(String, PathBuf)>,
        /// name=bed, mehrfach: Ergebnis zusätzlich für Reads, deren wahrer Anfang dort liegt.
        #[arg(long = "bereich", value_parser = eval::parse_named)]
        regions: Vec<(String, PathBuf)>,
        /// Größter erlaubter Abstand zum wahren linken Ende.
        #[arg(long, default_value_t = 20)]
        toleranz: u64,
        #[arg(long)]
        json: Option<PathBuf>,
        /// Liste der sicher (MAPQ ab 20) falschen Reads als TSV, zum Ansehen jedes Falls.
        #[arg(long)]
        falsche: Option<PathBuf>,
        #[arg(long, default_value_t = 0)]
        kerne: usize,
    },
}

fn pct(a: u64, b: u64) -> String {
    if b == 0 {
        "-".into()
    } else {
        format!("{:.4} %", 100.0 * a as f64 / b as f64)
    }
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::ReadsErzeugen {
            referenz,
            ausgabe,
            paare,
            laenge,
            insert,
            streuung,
            saat,
            profil,
            snv_rate,
            indel_rate,
            max_indel,
        } => {
            let settings = sim::Settings {
                pairs: paare,
                length: laenge,
                insert_mean: insert,
                insert_sd: streuung,
                seed: saat,
                snv_rate,
                indel_rate,
                max_indel,
                profile: match profil {
                    ProfileArg::Illumina => sim::ILLUMINA,
                    ProfileArg::Mgi => sim::MGI,
                },
            };
            let s = sim::run(&referenz, &ausgabe, &settings)?;
            let result = json!({
                "referenz": referenz,
                "einstellungen": settings,
                "contigs": s.contigs,
                "varianten": s.variants,
                "snvs": s.snvs,
                "indels": s.indels,
                "paare": s.pairs,
                "verworfene_stuecke": s.rejected,
            });
            let env = gc_schema::Envelope::for_tool(
                "genome-validate",
                "genome.validate.reads.v1",
                result,
            );
            let mut p = ausgabe.into_os_string();
            p.push(".json");
            std::fs::write(&p, env.to_json()? + "\n")?;
            println!(
                "{} Paare, {} Varianten ({} SNVs, {} Indels), {} Stücke mit N oder am Contig-Ende verworfen",
                s.pairs, s.variants, s.snvs, s.indels, s.rejected
            );
        }
        Command::AusrichtungPruefen {
            contigs,
            bams,
            regions,
            toleranz,
            json: json_out,
            falsche,
            kerne,
        } => {
            let names = eval::read_contigs(&contigs)?;
            let regions = regions
                .iter()
                .map(|(n, p)| eval::Regions::load(n, p))
                .collect::<Result<Vec<_>>>()?;
            let mut list = match &falsche {
                Some(p) => {
                    let mut f = std::io::BufWriter::new(std::fs::File::create(p)?);
                    writeln!(f, "bam\tread\tnr\twahr\tgefunden\tmapq")?;
                    Some(f)
                }
                None => None,
            };
            let mut evals = Vec::new();
            for (label, path) in &bams {
                let w = list.as_mut().map(|f| f as &mut dyn std::io::Write);
                let ev = eval::evaluate(label, path, &names, &regions, toleranz, kerne, w)?;
                if ev.all.reads == 0 {
                    bail!("{}: keine Reads mit Wahrheit im Namen", path.display());
                }
                evals.push(ev);
            }
            let mut out = Vec::new();
            for ev in &evals {
                let a = &ev.all;
                println!("== {}", ev.label);
                println!(
                    "Reads {}  ohne Ort {}  falsch {} ({})  MAPQ ab 20: {} davon falsch {} ({})",
                    a.reads,
                    a.unmapped,
                    a.wrong,
                    pct(a.wrong, a.reads),
                    a.mapq20,
                    a.wrong_mapq20,
                    pct(a.wrong_mapq20, a.mapq20)
                );
                println!("MAPQ     Reads      falsch   erwartet");
                let bins = eval::bins(ev);
                for b in &bins {
                    println!(
                        "{:>2}-{:<2} {:>10} {:>10} {:>10.1}",
                        b.mapq_from, b.mapq_to, b.reads, b.wrong, b.expected_wrong
                    );
                }
                for (r, c) in regions.iter().zip(&ev.regions) {
                    println!(
                        "Bereich {}: Reads {}  falsch {} ({})  MAPQ ab 20 falsch {} von {} ({})",
                        r.name,
                        c.reads,
                        c.wrong,
                        pct(c.wrong, c.reads),
                        c.wrong_mapq20,
                        c.mapq20,
                        pct(c.wrong_mapq20, c.mapq20)
                    );
                }
                let per_mapq: Vec<_> = ev
                    .by_mapq
                    .iter()
                    .enumerate()
                    .filter(|(_, (n, _))| *n > 0)
                    .map(|(q, (n, w))| json!({"mapq": q, "reads": n, "falsch": w}))
                    .collect();
                out.push(json!({
                    "name": ev.label,
                    "gesamt": ev.all,
                    "stufen": bins,
                    "je_mapq": per_mapq,
                    "bereiche": regions.iter().zip(&ev.regions).map(|(r, c)| json!({
                        "name": r.name, "basen": r.bases, "zahlen": c
                    })).collect::<Vec<_>>(),
                    "uebersprungen_sekundaer_supplementaer": ev.records_skipped,
                    "namen_ohne_wahrheit": ev.names_unknown,
                }));
            }
            let mut pairs = Vec::new();
            for other in evals.iter().skip(1) {
                let p = eval::compare(&evals[0], other);
                println!(
                    "{} gegen {}: beide richtig {}, nur {} richtig {}, nur {} richtig {}, beide nicht {}; sicher falsch nur in {}: {}, nur in {}: {}",
                    p.a,
                    p.b,
                    p.both_right,
                    p.a,
                    p.only_a_right,
                    p.b,
                    p.only_b_right,
                    p.both_not_right,
                    p.a,
                    p.a_confident_wrong_b_right,
                    p.b,
                    p.b_confident_wrong_a_right
                );
                pairs.push(p);
            }
            if let Some(mut f) = list {
                f.flush()?;
            }
            if let Some(path) = json_out {
                let result = json!({
                    "toleranz": toleranz,
                    "bams": out,
                    "vergleiche": pairs,
                });
                let env = gc_schema::Envelope::for_tool(
                    "genome-validate",
                    "genome.validate.mapping.v1",
                    result,
                );
                std::fs::write(path, env.to_json()? + "\n")?;
            }
        }
    }
    Ok(())
}
