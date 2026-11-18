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
    /// FM-Index über ein Referenzgenom bauen, prüfen und abfragen (Meilenstein M2).
    Index {
        #[command(subcommand)]
        command: IndexCommand,
    },
}

#[derive(Subcommand)]
enum IndexCommand {
    /// Index aus FASTA bauen. Schreibt daneben das Lauf-Protokoll `<index>.log.jsonl`.
    Build {
        fasta: PathBuf,
        index: PathBuf,
        /// `kompakt` (wenig RAM) oder `schnell`
        #[arg(long, default_value = "kompakt")]
        variante: String,
        /// Nur diese Contigs aufnehmen (mehrfach möglich); ohne Angabe alle.
        #[arg(long = "contig")]
        contigs: Vec<String>,
        /// Selbstprüfung: Suffix-Array vollständig prüfen und so viele Zufallsabfragen
        /// gegen die binäre Suche im vollen Suffix-Array vergleichen (0 = aus).
        #[arg(long, default_value_t = 100_000)]
        pruefen: u64,
        /// Bauverfahren: `direkt` (volles Suffix-Array im RAM, für einzelne Chromosomen),
        /// `blockweise` (fester RAM-Bedarf, für das ganze Genom) oder `auto`.
        #[arg(long, default_value = "auto")]
        verfahren: String,
        /// Blockweise: Millionen Suffixe je Block (16 MB RAM je Million).
        #[arg(long, default_value_t = 96)]
        block_mio: u64,
        /// Kerne (0 = alle)
        #[arg(long, default_value_t = 0)]
        kerne: usize,
        /// Vollprüfung: die BWT einmal rückwärts durchgehen; sie muss Base für Base die Referenz
        /// ergeben und jede Suffix-Array-Stichprobe bestätigen. Dauert beim ganzen Genom lange.
        #[arg(long)]
        vollpruefung: bool,
        /// Protokollstufe: error, warn, info, debug, trace
        #[arg(long, default_value = "info")]
        log_stufe: String,
    },
    /// Exakte Treffer zählen: eine Sequenz je Zeile aus der Datei, Ausgabe `Sequenz<TAB>Anzahl`
    /// (beide Stränge, ohne Filter, wie die Breite des Suffix-Array-Bereichs).
    Count { index: PathBuf, sequences: PathBuf },
    /// Treffer einer Sequenz mit Contig, Position und Strang als JSON.
    Locate {
        index: PathBuf,
        sequence: String,
        #[arg(long, default_value_t = 1000)]
        limit: u64,
    },
}

fn load_index(path: &std::path::Path) -> Result<gc_index::FmIndex> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("Index {} nicht lesbar", path.display()))?;
    gc_index::FmIndex::read(std::io::BufReader::with_capacity(1 << 20, file))
        .with_context(|| format!("Index {}", path.display()))
}

fn encode_sequence(seq: &str) -> Result<Vec<u8>> {
    seq.bytes()
        .enumerate()
        .map(|(i, b)| {
            gc_index::reference::encode(b).with_context(|| {
                format!(
                    "Zeichen {} an Stelle {} ist keine der Basen A, C, G, T",
                    b as char,
                    i + 1
                )
            })
        })
        .collect()
}

struct IndexBuildOptions<'a> {
    variante: &'a str,
    wanted: &'a [String],
    pruefen: u64,
    verfahren: &'a str,
    block_mio: u64,
    kerne: usize,
    vollpruefung: bool,
}

fn index_build(
    fasta: &std::path::Path,
    out: &std::path::Path,
    o: &IndexBuildOptions,
) -> Result<serde_json::Value> {
    let IndexBuildOptions {
        variante,
        wanted,
        pruefen,
        verfahren,
        block_mio,
        kerne,
        vollpruefung,
    } = *o;
    let variant = match variante {
        "kompakt" => gc_index::Variant::Kompakt,
        "schnell" => gc_index::Variant::Schnell,
        other => anyhow::bail!("unbekannte Variante {other}: erlaubt sind kompakt und schnell"),
    };
    gc_log::input_file("referenz", fasta);

    let mut st = gc_log::stage(
        "index.fasta_lesen",
        json!({ "pfad": fasta.display().to_string() }),
    );
    let mut reader = gc_io::fasta::open(fasta)?;
    let mut record = gc_io::fasta::Record::default();
    let mut reference = gc_index::Reference::new();
    let mut seen = 0u64;
    while reader.read_record(&mut record)? {
        seen += 1;
        if wanted.is_empty() || wanted.contains(&record.name) {
            reference.push_contig(&record.name, &record.seq);
        }
        st.heartbeat(
            5,
            || json!({ "contigs_gelesen": seen, "basen_aufgenommen": reference.len() }),
        );
    }
    for w in wanted {
        if !reference.contigs.iter().any(|c| c.name == *w) {
            gc_log::failure(gc_log::Failure {
                code: "INDEX.CONTIG_FEHLT",
                what: &format!("Contig {w} steht nicht in der Referenz"),
                location: json!({ "datei": fasta.display().to_string(), "contigs_in_datei": seen }),
                cause: "Name falsch geschrieben oder andere Referenz (z. B. 20 statt chr20)",
                action: "Namen in der FASTA-Datei prüfen (Zeilen, die mit > beginnen)",
            });
            anyhow::bail!("Contig {w} nicht in {}", fasta.display());
        }
    }
    st.finish(json!({ "contigs_in_datei": seen, "contigs_aufgenommen": reference.contigs.len(), "basen": reference.len() }));

    // Der direkte Bau hält das volle Suffix-Array im RAM (rund 10 Bytes je Referenzbase).
    const DIREKT_BIS_BASEN: u64 = 300_000_000;
    let blockweise = match verfahren {
        "direkt" => false,
        "blockweise" => true,
        "auto" => {
            let b = reference.len() > DIREKT_BIS_BASEN;
            gc_log::decision(
                "INDEX.VERFAHREN",
                if b { "blockweise" } else { "direkt" },
                if b {
                    "Referenz zu groß für ein volles Suffix-Array im RAM"
                } else {
                    "Referenz klein genug für das volle Suffix-Array im RAM"
                },
                json!({ "referenz_basen": reference.len(), "grenze_basen": DIREKT_BIS_BASEN }),
            );
            b
        }
        other => {
            anyhow::bail!("unbekanntes Verfahren {other}: erlaubt sind auto, direkt und blockweise")
        }
    };
    let built = if blockweise {
        let options = gc_index::blockwise::Options {
            threads: kerne,
            block_suffixes: block_mio.max(1) * 1_000_000,
        };
        gc_index::blockwise::build(reference, variant, options).map(|index| (index, None))
    } else {
        gc_index::FmIndex::build(reference, variant, pruefen)
    };
    let (index, mut report) = match built {
        Ok(v) => v,
        Err(e) => {
            let (cause, action) = match &e {
                gc_index::BuildError::Empty => (
                    "Die FASTA-Datei enthält keine Basen oder kein gewähltes Contig",
                    "Datei und --contig prüfen",
                ),
                gc_index::BuildError::TooLong { .. } => (
                    "Der direkte Bau arbeitet mit 32-Bit-Positionen",
                    "--verfahren blockweise (oder auto) verwenden",
                ),
                _ => (
                    "Fehler im Programm, nicht in den Eingabedaten",
                    "Protokoll aufbewahren und melden",
                ),
            };
            gc_log::failure(gc_log::Failure {
                code: "INDEX.BAU_FEHLGESCHLAGEN",
                what: &e.to_string(),
                location: json!({ "datei": fasta.display().to_string() }),
                cause,
                action,
            });
            return Err(e.into());
        }
    };

    let check_failed = |what: &str, detail: &str| {
        gc_log::failure(gc_log::Failure {
            code: "INDEX.PRUEFUNG_FEHLGESCHLAGEN",
            what: &format!("{what}: {detail}"),
            location: json!({ "index": out.display().to_string() }),
            cause: "Fehler im Programm oder beschädigter Speicher, nicht in den Eingabedaten",
            action: "Index nicht verwenden; Protokoll aufbewahren und melden",
        });
    };
    if blockweise && pruefen > 0 {
        let st = gc_log::stage("index.stichprobenpruefung", json!({ "abfragen": pruefen }));
        match index.check_against_reference(pruefen) {
            Ok(r) => {
                st.finish(json!({
                    "abfragen": r.queries, "abfragen_mit_treffern": r.queries_with_hits,
                    "gepruefte_fundstellen": r.occurrences_compared, "abweichungen": 0,
                }));
                report = Some(r);
            }
            Err(detail) => {
                check_failed("Stichprobenprüfung", &detail);
                anyhow::bail!("Stichprobenprüfung fehlgeschlagen: {detail}");
            }
        }
    }
    let mut inversion = None;
    if vollpruefung {
        let mut st = gc_log::stage("index.vollpruefung", json!({ "zeichen": index.rows() - 1 }));
        let total = index.rows() - 1;
        let outcome = index.verify_by_inversion(|done| {
            st.heartbeat(
                30,
                || json!({ "zeichen_geprueft": done, "anteil": done as f64 / total as f64 }),
            );
        });
        match outcome {
            Ok(n) => {
                st.finish(json!({ "zeichen_geprueft": n, "abweichungen": 0 }));
                inversion = Some(n);
            }
            Err(detail) => {
                check_failed("Vollprüfung", &detail);
                anyhow::bail!("Vollprüfung fehlgeschlagen: {detail}");
            }
        }
    }

    let st = gc_log::stage(
        "index.schreiben",
        json!({ "pfad": out.display().to_string() }),
    );
    // Erst in eine Nachbardatei schreiben, dann umbenennen: ein Abbruch hinterlässt nie einen halben Index.
    let mut tmp = out.as_os_str().to_owned();
    tmp.push(".teil");
    let tmp = PathBuf::from(tmp);
    let file =
        std::fs::File::create(&tmp).with_context(|| format!("{} nicht anlegbar", tmp.display()))?;
    index.write(std::io::BufWriter::with_capacity(1 << 20, file))?;
    std::fs::rename(&tmp, out).with_context(|| format!("{} → {}", tmp.display(), out.display()))?;
    let bytes = std::fs::metadata(out)?.len();
    st.finish(json!({ "bytes": bytes }));

    Ok(json!({
        "index": out.display().to_string(),
        "variante": variant.name(),
        "referenz_basen": index.reference().len(),
        "contigs": index.reference().contigs.len(),
        "n_strecken": index.reference().holes.len(),
        "zeilen": index.rows(),
        "basen_je_zaehlblock": index.occ_interval(),
        "sa_abstand": index.sa_interval(),
        "datei_bytes": bytes,
        "bits_je_base": 8.0 * index.table_bytes() as f64 / index.reference().len() as f64,
        "verfahren": if blockweise { "blockweise" } else { "direkt" },
        "vollpruefung": inversion.map(|n| json!({
            "zeichen_geprueft": n, "abweichungen": 0,
            "aussage": "BWT ergibt rückwärts Base für Base die Referenz, jede Suffix-Array-Stichprobe bestätigt",
        })),
        "selbstpruefung": report.map(|r| json!({
            "art": if blockweise { "Stichproben gegen die Referenz" } else { "Suffix-Array vollständig geprüft, Abfragen gegen binäre Suche" },
            "abfragen": r.queries,
            "abfragen_mit_treffern": r.queries_with_hits,
            "verglichene_fundstellen": r.occurrences_compared,
            "abweichungen": 0,
        })),
    }))
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
        Command::Index { command } => match command {
            IndexCommand::Build {
                fasta,
                index,
                variante,
                contigs,
                pruefen,
                verfahren,
                block_mio,
                kerne,
                vollpruefung,
                log_stufe,
            } => {
                let started = std::time::Instant::now();
                let level = gc_log::Level::parse(&log_stufe)
                    .with_context(|| format!("unbekannte Protokollstufe {log_stufe}"))?;
                let mut log_path = index.as_os_str().to_owned();
                log_path.push(".log.jsonl");
                gc_log::init(gc_log::Config {
                    tool: "index.build".to_string(),
                    version: env!("CARGO_PKG_VERSION").to_string(),
                    level,
                    jsonl: Some(PathBuf::from(log_path)),
                    stderr: true,
                })?;
                gc_log::install_panic_hook();
                gc_log::run_header(json!({
                    "fasta": fasta.display().to_string(), "index": index.display().to_string(),
                    "variante": variante, "contigs": contigs, "pruefen": pruefen, "log_stufe": log_stufe,
                    "verfahren": verfahren, "block_mio": block_mio, "kerne": kerne, "vollpruefung": vollpruefung,
                }));
                let options = IndexBuildOptions {
                    variante: &variante,
                    wanted: &contigs,
                    pruefen,
                    verfahren: &verfahren,
                    block_mio,
                    kerne,
                    vollpruefung,
                };
                match index_build(&fasta, &index, &options) {
                    Ok(mut result) => {
                        result["seconds"] = json!(started.elapsed().as_secs_f64());
                        result["lauf"] = json!(gc_log::run_id());
                        gc_log::run_end(true, result.clone());
                        println!(
                            "{}",
                            Envelope::new("genome.index-build.v1", result).to_json()?
                        );
                    }
                    Err(e) => {
                        gc_log::run_end(false, json!({ "fehler": format!("{e:#}") }));
                        return Err(e);
                    }
                }
            }
            IndexCommand::Count { index, sequences } => {
                use std::io::{BufRead, Write};
                let fm = load_index(&index)?;
                let input = std::io::BufReader::new(
                    std::fs::File::open(&sequences)
                        .with_context(|| format!("{} nicht lesbar", sequences.display()))?,
                );
                let mut out = std::io::BufWriter::new(std::io::stdout().lock());
                for (n, line) in input.lines().enumerate() {
                    let line = line?;
                    let seq = line.trim();
                    if seq.is_empty() {
                        continue;
                    }
                    let pattern = encode_sequence(seq)
                        .with_context(|| format!("{}: Zeile {}", sequences.display(), n + 1))?;
                    writeln!(out, "{seq}\t{}", fm.count(&pattern))?;
                }
                out.flush()?;
            }
            IndexCommand::Locate {
                index,
                sequence,
                limit,
            } => {
                let fm = load_index(&index)?;
                let pattern = encode_sequence(&sequence)?;
                let hits: Vec<_> = fm
                    .hits(&pattern, limit)
                    .iter()
                    .map(|h| {
                        json!({
                            "contig": fm.reference().contigs[h.contig as usize].name,
                            "pos": h.pos + 1,
                            "strang": if h.reverse { "-" } else { "+" },
                        })
                    })
                    .collect();
                let result = json!({ "fundstellen_roh": fm.count(&pattern), "treffer": hits });
                println!(
                    "{}",
                    Envelope::new("genome.index-locate.v1", result).to_json()?
                );
            }
        },
    }
    Ok(())
}
