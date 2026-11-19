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
    /// BAM nach Ort sortieren (Contig, Position, Strang; ohne Ort zuletzt), mit festem Speicher.
    Sort {
        src: PathBuf,
        dst: PathBuf,
        /// Speicher für die Datensätze in MB; was nicht hineinpasst, geht sortiert in Teilen
        /// auf die Platte.
        #[arg(long, default_value_t = 2048)]
        speicher_mb: usize,
        /// Ordner für die Teile (Vorgabe: der Ordner der Ausgabe)
        #[arg(long)]
        temp: Option<PathBuf>,
        /// Kerne für BGZF (0 = alle)
        #[arg(long, default_value_t = 0)]
        kerne: usize,
        /// Lauf-Protokoll als JSON-Zeilen in diese Datei
        #[arg(long)]
        log: Option<PathBuf>,
    },
    /// Duplikate in einer nach Ort sortierten BAM markieren (Flag 1024), Verfahren wie
    /// `samtools markdup`. Die Eingabe braucht die Partner-Angaben `MC` und `ms`.
    Markdup {
        src: PathBuf,
        dst: PathBuf,
        /// Längster erwarteter Read einschließlich abgeschnittener Basen
        #[arg(long, default_value_t = 300)]
        max_laenge: i64,
        /// Kerne für BGZF (0 = alle)
        #[arg(long, default_value_t = 0)]
        kerne: usize,
        /// Lauf-Protokoll als JSON-Zeilen in diese Datei
        #[arg(long)]
        log: Option<PathBuf>,
    },
    /// Index `.bai` zu einer nach Ort sortierten BAM schreiben (wie `samtools index`).
    BamIndex {
        bam: PathBuf,
        /// Indexdatei (Vorgabe: `<bam>.bai`)
        #[arg(long)]
        aus: Option<PathBuf>,
        /// Kerne für BGZF (0 = alle)
        #[arg(long, default_value_t = 0)]
        kerne: usize,
        /// Lauf-Protokoll als JSON-Zeilen in diese Datei
        #[arg(long)]
        log: Option<PathBuf>,
    },
    /// Datensätze einer BAM nach ihren Flags zählen (wie `samtools flagstat`).
    Flagstat {
        bam: PathBuf,
        /// Kerne für BGZF (0 = alle)
        #[arg(long, default_value_t = 0)]
        kerne: usize,
    },
    /// Vom FASTQ zur sortierten, duplikatmarkierten, indizierten BAM in einem Lauf.
    Fq2bam {
        /// FM-Index des Referenzgenoms (`genome index build`)
        index: PathBuf,
        /// FASTQ mit den Reads (bei Paaren: Read 1)
        #[arg(long)]
        fastq: PathBuf,
        /// FASTQ mit den Partnern (Read 2)
        #[arg(long)]
        fastq2: Option<PathBuf>,
        /// Ausgabe-BAM; daneben entstehen `<aus>.bai` und, ohne `--log`, `<aus>.log.jsonl`
        #[arg(long)]
        aus: PathBuf,
        /// Lesegruppe wie bei `bwa mem -R`, zum Beispiel '@RG\tID:lauf1\tSM:probe\tPL:ILLUMINA'
        #[arg(long)]
        lesegruppe: Option<String>,
        /// Speicher für das Sortieren in MB
        #[arg(long, default_value_t = 2048)]
        speicher_mb: usize,
        /// Ordner für die Teile des Sortierens (Vorgabe: der Ordner der Ausgabe)
        #[arg(long)]
        temp: Option<PathBuf>,
        /// Kerne (0 = alle). Die Ausgabe hängt nicht davon ab.
        #[arg(long, default_value_t = 0)]
        kerne: usize,
        /// Paare: Die Insertgröße wird je Block von so vielen Basen geschätzt (wie bwa `-K`)
        #[arg(long, default_value_t = 10_000_000)]
        block_basen: u64,
        /// Längster erwarteter Read für das Markieren der Duplikate
        #[arg(long, default_value_t = 300)]
        max_laenge: i64,
        /// Lauf-Protokoll als JSON-Zeilen in diese Datei
        #[arg(long)]
        log: Option<PathBuf>,
        /// Paare vor dem Ausrichten trimmen wie fastp 1.4.0 (Adapter, PolyG, Qualität, Länge).
        /// Ohne die Option bleiben die Reads unverändert.
        #[arg(long)]
        trimmen: bool,
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
        /// `standard` (Vorgabe), `kompakt` (wenig RAM, langsamer) oder `schnell` (mehr RAM)
        #[arg(long, default_value = "standard")]
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
    /// Seeds eines Reads: längste exakte Treffer (SMEM) mit Lage im Read, Häufigkeit und Orten.
    /// Mit `--fastq` stattdessen eine Zusammenfassung über alle Reads der Datei.
    Smem {
        index: PathBuf,
        /// Einzelne Sequenz (A, C, G, T, N)
        #[arg(long)]
        sequenz: Option<String>,
        /// FASTQ-Datei, Ausgabe ist dann nur die Zusammenfassung
        #[arg(long)]
        fastq: Option<PathBuf>,
        /// Kürzester gemeldeter Treffer (BWA-MEM: 19)
        #[arg(long, default_value_t = 19)]
        min_laenge: u32,
        /// Höchstens so viele Orte je Treffer ausgeben
        #[arg(long, default_value_t = 20)]
        orte: u64,
        /// Mit `--fastq`: jeden Treffer mit der einfachen Suche nachprüfen (Bereich gleich,
        /// links und rechts nicht verlängerbar, kein Treffer in einem anderen enthalten)
        #[arg(long)]
        pruefen: bool,
    },
    /// Ketten aus Seeds für alle Reads einer FASTQ-Datei, als Zusammenfassung. Mit
    /// `--wahrheit-im-namen` wird geprüft, ob die schwerste Kette am wahren Ort liegt; der Name
    /// muss dann `contig:pos:strang:...` lauten (wie von `bench/m3/simulate_reads.py`).
    Ketten {
        index: PathBuf,
        #[arg(long)]
        fastq: PathBuf,
        #[arg(long)]
        wahrheit_im_namen: bool,
        /// Je Treffer höchstens so viele Fundstellen orten (BWA-MEM: 500)
        #[arg(long, default_value_t = 500)]
        max_orte: u64,
        /// Namen der Reads, bei denen kein Seed genau am wahren Ort liegt, in diese Datei schreiben
        #[arg(long)]
        fehlende: Option<PathBuf>,
        /// Je Read den primären Bereich als Tabelle schreiben: Name, Contig, Position (1-basiert,
        /// linkes Ende auf dem Vorwärtsstrang), Strang, MAPQ, Score, gemeldeter Score,
        /// Read-Anfang, Read-Ende, sub, sub_n, Anteil in sehr häufigen Seeds
        #[arg(long)]
        tabelle: Option<PathBuf>,
    },
    /// Reads ausrichten und als SAM auf die Standardausgabe schreiben (ein Kern): einzeln aus
    /// einer FASTQ-Datei oder als Paare aus zwei. Ein Prüfbefehl für die Entwicklung, Vorstufe
    /// von `fq2bam`.
    Sam {
        index: PathBuf,
        #[arg(long)]
        fastq: PathBuf,
        /// Zweite Datei mit den Partnern: Die Reads werden als Paare ausgerichtet.
        #[arg(long)]
        fastq2: Option<PathBuf>,
        /// Paare: Die Insertgröße wird je Block von so vielen Basen geschätzt (wie bwa `-K`).
        #[arg(long, default_value_t = 10_000_000)]
        block_basen: u64,
        /// Protokoll zusätzlich als JSONL in diese Datei (sonst nur lesbar auf der Fehlerausgabe).
        #[arg(long)]
        log: Option<PathBuf>,
        /// Zahl der Kerne; 0 heißt alle. Die Ausgabe hängt nicht davon ab.
        #[arg(long, default_value_t = 0)]
        kerne: usize,
        /// Als BAM in diese Datei schreiben statt als SAM-Text auf die Standardausgabe.
        #[arg(long)]
        bam: Option<PathBuf>,
        /// Lesegruppe als Kopfzeile, Tabulatoren als \t geschrieben (wie bwa `-R`), zum Beispiel
        /// '@RG\tID:lauf1\tSM:probe\tPL:ILLUMINA'. Jede Zeile bekommt die Kennung als `RG`.
        #[arg(long)]
        lesegruppe: Option<String>,
        /// Paare: Angaben zum Partner wie nach `samtools fixmate -m` ergänzen (`ms`, `MC` auch
        /// am Read ohne Ort). Nötig für `genome markdup`.
        #[arg(long)]
        partner_angaben: bool,
    },
}

/// Wohin die Zeilen gehen: SAM-Text auf die Standardausgabe oder BAM in eine Datei.
enum AlignOutput {
    Sam(std::io::BufWriter<std::io::Stdout>),
    Bam(Box<gc_io::bam::MtWriter>),
    /// Direkt in den Sortierer, ohne Zwischendatei (`fq2bam`).
    Sorter(Box<gc_bam::Sorter>),
}

/// Kodiert Zeilen als BAM-Sätze.
fn encode_lines<'a>(lines: impl Iterator<Item = &'a gc_align::Line>) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    for line in lines {
        // Das Längenfeld des Namens in BAM hat ein Byte (SAM-Spezifikation: 254 Zeichen).
        anyhow::ensure!(
            line.name.len() <= 254,
            "Read-Name mit {} Zeichen ist für BAM zu lang (höchstens 254): {}",
            line.name.len(),
            line.name
        );
        line.write_bam(&mut bytes);
    }
    Ok(bytes)
}

impl AlignOutput {
    fn open(
        fm: &gc_index::FmIndex,
        read_group: Option<&gc_align::ReadGroup>,
        bam: Option<&std::path::Path>,
        threads: usize,
    ) -> Result<Self> {
        use std::io::Write;
        let command_line = std::env::args().collect::<Vec<_>>().join(" ");
        let header = gc_align::sam_header(fm, read_group, env!("CARGO_PKG_VERSION"), &command_line);
        Ok(match bam {
            Some(path) => AlignOutput::Bam(Box::new(
                gc_io::bam::create_from_text(path, &header, threads)
                    .with_context(|| format!("{} nicht schreibbar", path.display()))?,
            )),
            None => {
                let mut out = std::io::BufWriter::with_capacity(1 << 20, std::io::stdout());
                out.write_all(header.as_bytes())?;
                AlignOutput::Sam(out)
            }
        })
    }

    /// Schreibt die Zeilen eines Blocks; liefert ihre Zahl.
    fn write(&mut self, fm: &gc_index::FmIndex, block: &[Vec<gc_align::Line>]) -> Result<u64> {
        use std::io::Write;
        let lines = block.iter().flatten();
        match self {
            AlignOutput::Sam(out) => {
                let mut text = String::new();
                for line in lines {
                    line.write_sam(&mut text, fm);
                }
                out.write_all(text.as_bytes())?;
            }
            AlignOutput::Bam(writer) => {
                gc_io::bam::write_encoded(writer, &encode_lines(lines)?)?;
            }
            AlignOutput::Sorter(sorter) => {
                let runs = sorter.stats().runs;
                sorter
                    .push(&encode_lines(lines)?)
                    .context("sortierten Teil auf die Platte schreiben")?;
                let now = sorter.stats();
                if now.runs != runs {
                    gc_log::info(
                        "SORT.TEIL",
                        "Speicher voll, sortierter Teil auf die Platte geschrieben",
                        json!({ "teile": now.runs, "bytes_auf_platte": now.run_bytes, "saetze_bisher": now.records }),
                    );
                }
            }
        }
        Ok(block.iter().map(|l| l.len() as u64).sum())
    }

    fn finish(self) -> Result<()> {
        use std::io::Write;
        match self {
            AlignOutput::Sam(mut out) => out.flush()?,
            AlignOutput::Bam(writer) => gc_io::bam::finish(*writer)?,
            AlignOutput::Sorter(_) => {}
        }
        Ok(())
    }
}

/// Einstellungen eines Laufs von `genome index sam`.
#[derive(Clone, Copy)]
struct SamRun<'a> {
    block_bases: u64,
    threads: usize,
    bam: Option<&'a std::path::Path>,
    read_group: Option<&'a gc_align::ReadGroup>,
    mate_tags: bool,
    /// Paare vor dem Ausrichten trimmen; verworfene Paare werden nicht ausgerichtet.
    trimmer: Option<&'a gc_io::trim::Trimmer>,
}

/// Richtet Reads aus, einzeln oder als Paare, Block für Block auf `threads` Kernen.
/// Bei Paaren wird die Insertgröße je Block geschätzt.
fn index_sam(
    fm: &gc_index::FmIndex,
    first: &std::path::Path,
    second: Option<&std::path::Path>,
    run: &SamRun<'_>,
) -> Result<()> {
    let SamRun {
        threads,
        bam,
        read_group,
        ..
    } = *run;
    let mut out = AlignOutput::open(fm, read_group, bam, threads)?;
    align_into(fm, first, second, run, &mut out)?;
    out.finish()
}

/// Richtet aus und gibt die Zeilen Block für Block an `out`. Liefert die Zahl der Reads und
/// der Zeilen.
fn align_into(
    fm: &gc_index::FmIndex,
    first: &std::path::Path,
    second: Option<&std::path::Path>,
    run: &SamRun<'_>,
    out: &mut AlignOutput,
) -> Result<(u64, u64, gc_io::trim::Counts)> {
    let SamRun {
        block_bases,
        threads,
        read_group,
        mate_tags,
        trimmer,
        ..
    } = *run;
    let mut trim_counts = gc_io::trim::Counts::default();
    let mut opt = gc_align::AlignOptions::default();
    opt.sam.read_group = read_group.map(|rg| rg.id.clone());
    opt.sam.mate_tags = mate_tags;
    let mut reader1 = fastq::open(first)?;
    let mut reader2 = match second {
        Some(p) => Some(fastq::open(p)?),
        None => None,
    };
    let paired = reader2.is_some();
    let started = std::time::Instant::now();
    let (mut n_reads, mut n_lines, mut block_nr) = (0u64, 0u64, 0u64);
    loop {
        // Einen Block lesen, bis die Zahl der Basen erreicht ist.
        let mut records: Vec<Record> = Vec::new();
        let mut bases = 0u64;
        while bases < block_bases {
            let mut a = Record::default();
            let more = reader1.read_record(&mut a)?;
            match reader2.as_mut() {
                Some(r2) => {
                    let mut b = Record::default();
                    let more2 = r2.read_record(&mut b)?;
                    if more != more2 {
                        anyhow::bail!(
                            "die beiden FASTQ-Dateien haben verschieden viele Reads ({} und {})",
                            first.display(),
                            r2.path().display()
                        );
                    }
                    if !more {
                        break;
                    }
                    if a.pair_key() != b.pair_key() {
                        anyhow::bail!(
                            "Paar {}: Namen passen nicht zusammen ({} und {})",
                            n_reads / 2 + records.len() as u64 / 2 + 1,
                            String::from_utf8_lossy(&a.name),
                            String::from_utf8_lossy(&b.name)
                        );
                    }
                    if let Some(t) = trimmer {
                        let report = t.trim_pair(&mut a, &mut b);
                        trim_counts.add(&report);
                        if !report.kept() {
                            continue;
                        }
                    }
                    bases += (a.seq.len() + b.seq.len()) as u64;
                    records.push(a);
                    records.push(b);
                }
                None => {
                    if !more {
                        break;
                    }
                    bases += a.seq.len() as u64;
                    records.push(a);
                }
            }
        }
        if records.is_empty() {
            break;
        }
        let seqs: Vec<Vec<u8>> = records.iter().map(|r| encode_read(&r.seq)).collect();
        let names: Vec<String> = records
            .iter()
            .map(|r| {
                let name = if paired { r.pair_key() } else { &r.name[..] };
                String::from_utf8_lossy(name).into_owned()
            })
            .collect();
        let inputs: Vec<gc_align::ReadInput<'_>> = (0..records.len())
            .map(|i| gc_align::ReadInput {
                name: &names[i],
                seq: &seqs[i],
                qual: &records[i].qual,
            })
            .collect();
        let lines = if paired {
            let (lines, pes) = gc_align::align_pairs(fm, &inputs, n_reads / 2, &opt, threads);
            let used: Vec<&str> = (0..4)
                .filter(|&d| !pes[d].failed)
                .map(|d| gc_align::DIRECTIONS[d])
                .collect();
            gc_log::decision(
                "ALIGN.INSERTGROESSE",
                &if used.is_empty() {
                    "keine Richtung".to_string()
                } else {
                    used.join(", ")
                },
                "aus den eindeutig liegenden Paaren dieses Blocks geschätzt",
                json!({
                    "block": block_nr,
                    "paare": records.len() / 2,
                    "richtungen": (0..4).map(|d| json!({
                        "richtung": gc_align::DIRECTIONS[d],
                        "paare": pes[d].n,
                        "benutzt": !pes[d].failed,
                        "perzentile_25_50_75": pes[d].percentiles,
                        "mittel": pes[d].avg,
                        "streuung": pes[d].std,
                        "richtiges_paar_von": pes[d].low,
                        "richtiges_paar_bis": pes[d].high,
                    })).collect::<Vec<_>>(),
                }),
            );
            lines
        } else {
            gc_align::align_reads(fm, &inputs, n_reads, &opt, threads)
        };
        n_lines += out.write(fm, &lines)?;
        n_reads += records.len() as u64;
        block_nr += 1;
    }
    let seconds = started.elapsed().as_secs_f64();
    gc_log::info(
        "ALIGN.FERTIG",
        "Ausrichten beendet",
        json!({
            "reads": n_reads,
            "zeilen": n_lines,
            "bloecke": block_nr,
            "kerne": threads,
            "sekunden": seconds,
            "reads_je_sekunde": if seconds > 0.0 { n_reads as f64 / seconds } else { 0.0 },
        }),
    );
    Ok((n_reads, n_lines, trim_counts))
}

/// Einstellungen eines Laufs von `genome fq2bam`.
struct Fq2bamRun<'a> {
    index: &'a std::path::Path,
    fastq: &'a std::path::Path,
    fastq2: Option<&'a std::path::Path>,
    dst: &'a std::path::Path,
    read_group: Option<&'a str>,
    memory_mb: usize,
    temp: &'a std::path::Path,
    threads: usize,
    block_bases: u64,
    max_length: i64,
    trim: bool,
}

/// Zahlen von `flagstat` als JSON, für eine der beiden Gruppen.
fn flag_counts_json(c: &gc_bam::Counts) -> serde_json::Value {
    let share = |part: u64, whole: u64| {
        if whole > 0 {
            json!(part as f64 / whole as f64)
        } else {
            json!(null)
        }
    };
    json!({
        "gesamt": c.total,
        "hauptzeilen": c.primary,
        "sekundaer": c.secondary,
        "ergaenzend": c.supplementary,
        "duplikate": c.duplicates,
        "duplikate_hauptzeilen": c.primary_duplicates,
        "mit_ort": c.mapped,
        "anteil_mit_ort": share(c.mapped, c.total),
        "hauptzeilen_mit_ort": c.primary_mapped,
        "anteil_hauptzeilen_mit_ort": share(c.primary_mapped, c.primary),
        "aus_paaren": c.paired,
        "read1": c.read1,
        "read2": c.read2,
        "richtiges_paar": c.properly_paired,
        "anteil_richtiges_paar": share(c.properly_paired, c.paired),
        "beide_mit_ort": c.both_mapped,
        "nur_selbst_mit_ort": c.singletons,
        "anteil_nur_selbst_mit_ort": share(c.singletons, c.paired),
        "partner_auf_anderem_contig": c.mate_on_other_contig,
        "partner_auf_anderem_contig_mapq5": c.mate_on_other_contig_mapq5,
    })
}

/// Zahlen des Markierens als JSON.
fn dup_stats_json(stats: &gc_bam::DupStats) -> serde_json::Value {
    let duplicates = stats.duplicate_pair + stats.duplicate_single;
    json!({
        "gelesen": stats.read,
        "geschrieben": stats.written,
        "nicht_bewertet": stats.excluded,
        "bewertet": stats.examined,
        "mit_partner": stats.paired,
        "ohne_partner": stats.single,
        "duplikate_mit_partner": stats.duplicate_pair,
        "duplikate_ohne_partner": stats.duplicate_single,
        "duplikate": duplicates,
        "anteil_duplikate": if stats.examined > 0 { duplicates as f64 / stats.examined as f64 } else { 0.0 },
        "geschaetzte_bibliotheksgroesse": stats.estimated_library_size(),
        "fenster_spitze_saetze": stats.window_peak,
    })
}

/// Entscheidungen vor dem Trimmen wie fastp 1.4.0: PolyG nach dem Namen des ersten Reads, je
/// Seite die Adapterfolge aus dem Anfang der Datei. Beide stehen mit Grund im Protokoll.
fn plan_trimming(
    first: &std::path::Path,
    second: &std::path::Path,
) -> Result<gc_io::trim::Trimmer> {
    use gc_io::trim;
    let st = gc_log::stage("fq2bam.trimmen_vorbereiten", json!({}));
    let mut adapters: [Option<Vec<u8>>; 2] = [None, None];
    // Beide Seiten werden gleichzeitig gelesen und untersucht; das Ergebnis hängt
    // nicht davon ab, nur die Wartezeit.
    let first_name_slot = std::sync::Mutex::new(Vec::new());
    let read_side =
        |side: usize, path: &std::path::Path| -> Result<(usize, Result<trim::Adapter, String>)> {
            let mut reader = fastq::open(path)?;
            let (mut reads, mut bases) = (Vec::new(), 0usize);
            let mut rec = Record::default();
            while reads.len() < trim::DETECT_READS && bases < trim::DETECT_BASES {
                if !reader.read_record(&mut rec)? {
                    break;
                }
                if side == 0 && reads.is_empty() {
                    *first_name_slot.lock().expect("Sperre") = rec.name.clone();
                }
                bases += rec.seq.len();
                reads.push(rec.seq.clone());
            }
            let found = trim::detect_adapter(&reads);
            Ok((reads.len(), found))
        };
    let results = std::thread::scope(|scope| {
        let h1 = scope.spawn(|| read_side(0, first));
        let h2 = scope.spawn(|| read_side(1, second));
        [h1.join(), h2.join()]
    });
    let first_name = first_name_slot.into_inner().expect("Sperre");
    for (side, result) in results.into_iter().enumerate() {
        let (n_reads, found) = result.map_err(|_| anyhow::anyhow!("Adaptersuche abgestürzt"))??;
        let label = if side == 0 { "read1" } else { "read2" };
        match found {
            Ok(a) => {
                let how = match &a.found {
                    trim::Found::Known { hits, checked } => json!({
                        "weg": "bekannter_adapter", "treffer": hits, "gepruefte_reads": checked,
                    }),
                    trim::Found::Seed { seed, known } => json!({
                        "weg": "haeufigste_teilfolge", "keim": seed, "bekannt": known,
                    }),
                };
                gc_log::decision(
                    "TRIM.ADAPTER",
                    &String::from_utf8_lossy(&a.seq),
                    "aus dem Anfang der Datei erkannt wie fastp 1.4.0",
                    json!({
                        "seite": label,
                        "name": a.name,
                        "reads_angesehen": n_reads,
                        "erkennung": how,
                    }),
                );
                adapters[side] = Some(a.seq);
            }
            Err(reason) => gc_log::decision(
                "TRIM.ADAPTER",
                "keiner",
                &reason,
                json!({ "seite": label, "reads_angesehen": n_reads }),
            ),
        }
    }
    let prefix = trim::two_color(&first_name);
    gc_log::decision(
        "TRIM.POLYG",
        if prefix.is_some() { "an" } else { "aus" },
        &match prefix {
            Some(p) => format!("erster Read-Name beginnt mit {p}: Gerät mit zweifarbiger Chemie"),
            None => "erster Read-Name passt zu keinem Gerät mit zweifarbiger Chemie".to_string(),
        },
        json!({ "erster_name": String::from_utf8_lossy(&first_name) }),
    );
    st.finish(json!({}));
    Ok(trim::Trimmer {
        settings: trim::Settings::default(),
        poly_g: prefix.is_some(),
        adapters,
    })
}

/// Zähler des Trimmens als JSON.
fn trim_counts_json(c: &gc_io::trim::Counts) -> serde_json::Value {
    json!({
        "paare": c.pairs,
        "behalten": c.kept,
        "verworfen_qualitaet": c.dropped_quality,
        "verworfen_zu_kurz": c.dropped_length,
        "verworfen_zu_viele_n": c.dropped_n,
        "verworfen_adapter_dimer": c.dropped_dimer,
        "polyg_reads": c.poly_g_reads,
        "polyg_basen": c.poly_g_bases,
        "ueberlappung_gekuerzt_paare": c.overlap_pairs,
        "adapterfolge_gekuerzt_read1": c.sequence_reads[0],
        "adapterfolge_gekuerzt_read2": c.sequence_reads[1],
        "adapter_basen": c.adapter_bases,
    })
}

/// Die ganze Kette in einem Lauf: ausrichten (mit Partner-Angaben) direkt in den Sortierer,
/// beim Zusammenführen Duplikate markieren und nach Flags zählen, dann den Index schreiben.
/// Zwischen den Schritten entsteht keine BAM; auf der Platte liegen nur die Teile des Sortierens.
fn fq2bam(run: &Fq2bamRun<'_>) -> Result<serde_json::Value> {
    let started = std::time::Instant::now();
    anyhow::ensure!(run.memory_mb > 0, "--speicher-mb muss größer als 0 sein");
    anyhow::ensure!(run.max_length > 0, "--max-laenge muss größer als 0 sein");
    let read_group = match run.read_group {
        Some(text) => Some(
            gc_align::ReadGroup::parse(text)
                .map_err(|cause| anyhow::anyhow!("--lesegruppe: {cause}"))?,
        ),
        None => None,
    };
    let trimmer = if run.trim {
        let Some(second) = run.fastq2 else {
            anyhow::bail!("--trimmen gibt es nur für Paare (--fastq2 fehlt)");
        };
        Some(plan_trimming(run.fastq, second)?)
    } else {
        None
    };
    let st = gc_log::stage("fq2bam.index_laden", json!({}));
    let fm = load_index(run.index)?;
    st.finish(json!({}));
    let seconds_load = started.elapsed().as_secs_f64();

    // 1. Ausrichten, die Sätze gehen sofort in den Sortierer.
    let t = std::time::Instant::now();
    let st = gc_log::stage("fq2bam.ausrichten", json!({ "kerne": run.threads }));
    let mut out = AlignOutput::Sorter(Box::new(gc_bam::Sorter::new(gc_bam::SortOptions {
        memory_bytes: run.memory_mb.saturating_mul(1 << 20),
        temp_dir: run.temp.to_path_buf(),
    })));
    let sam_run = SamRun {
        block_bases: run.block_bases,
        threads: run.threads,
        bam: None,
        read_group: read_group.as_ref(),
        mate_tags: true,
        trimmer: trimmer.as_ref(),
    };
    let (reads, lines, trim_counts) = align_into(&fm, run.fastq, run.fastq2, &sam_run, &mut out)?;
    let trim_json = trimmer.as_ref().map(|_| trim_counts_json(&trim_counts));
    if let Some(j) = &trim_json {
        anyhow::ensure!(
            trim_counts.balanced(),
            "Zähler des Trimmens gehen nicht auf: {j}"
        );
        gc_log::info("TRIM.FERTIG", "Trimmen beendet", j.clone());
    }
    let AlignOutput::Sorter(sorter) = out else {
        unreachable!("oben als Sortierer angelegt")
    };
    let sorted = sorter.stats();
    st.finish(json!({ "reads": reads, "zeilen": lines, "teile": sorted.runs }));
    let seconds_align = t.elapsed().as_secs_f64();

    // 2. Zusammenführen, dabei Duplikate markieren, zählen und schreiben.
    let t = std::time::Instant::now();
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let header_text = gc_bam::with_coordinate_order(&gc_align::sam_header(
        &fm,
        read_group.as_ref(),
        env!("CARGO_PKG_VERSION"),
        &command_line,
    ));
    let mut part = run.dst.as_os_str().to_owned();
    part.push(".unfertig");
    let part = PathBuf::from(part);
    let st = gc_log::stage("fq2bam.zusammenfuehren", json!({ "teile": sorted.runs }));
    let written = (|| -> Result<(gc_bam::SortStats, gc_bam::DupStats, gc_bam::FlagStats)> {
        let mut writer = gc_io::bam::create_from_text(&part, &header_text, run.threads)
            .with_context(|| format!("{} anlegen", part.display()))?;
        let mut marker = gc_bam::DupMarker::new(gc_bam::DupOptions {
            max_length: run.max_length,
        });
        let mut flags = gc_bam::FlagStats::default();
        let mut sink = |rec: &[u8]| {
            flags.add(&rec[4..]);
            gc_io::bam::write_encoded(&mut writer, rec)
        };
        let sort_stats = sorter
            .finish(|rec| marker.push(rec, &mut sink))
            .context("Teile zusammenführen und Duplikate markieren")?;
        let dup_stats = marker.finish(&mut sink)?;
        gc_io::bam::finish(writer).with_context(|| format!("{} abschließen", part.display()))?;
        std::fs::rename(&part, run.dst)
            .with_context(|| format!("{} umbenennen", part.display()))?;
        Ok((sort_stats, dup_stats, flags))
    })();
    let (sort_stats, dup_stats, flags) = match written {
        Ok(v) => v,
        Err(e) => {
            let _ = std::fs::remove_file(&part);
            return Err(e);
        }
    };
    st.finish(json!({ "saetze": sort_stats.records, "duplikate": dup_stats.duplicate_pair + dup_stats.duplicate_single }));
    let seconds_merge = t.elapsed().as_secs_f64();

    // 3. Index.
    let t = std::time::Instant::now();
    let mut bai = run.dst.as_os_str().to_owned();
    bai.push(".bai");
    let bai = PathBuf::from(bai);
    let index_result = index_bam(run.dst, &bai, run.threads)?;
    let seconds_index = t.elapsed().as_secs_f64();

    let seconds = started.elapsed().as_secs_f64();
    let result = json!({
        "bam": run.dst.display().to_string(),
        "bai": bai.display().to_string(),
        "reads": reads,
        "saetze": sort_stats.records,
        "paare": run.fastq2.is_some(),
        "trimmen": trim_json,
        "kerne": run.threads,
        "sortieren": {
            "bytes_saetze": sort_stats.bytes,
            "teile": sort_stats.runs,
            "bytes_teile": sort_stats.run_bytes,
            "speicher_mb": run.memory_mb,
        },
        "duplikate": dup_stats_json(&dup_stats),
        "flagstat": {
            "bestanden": flag_counts_json(&flags.passed),
            "nicht_bestanden": flag_counts_json(&flags.failed),
        },
        "index": {
            "faecher": index_result["faecher"],
            "bytes_index": index_result["bytes_index"],
            "ohne_contig": index_result["ohne_contig"],
        },
        "sekunden": seconds,
        "sekunden_schritte": {
            "index_laden": seconds_load,
            "ausrichten_und_teile": seconds_align,
            "zusammenfuehren_markieren_schreiben": seconds_merge,
            "bam_index": seconds_index,
        },
        "reads_je_sekunde": if seconds > 0.0 { reads as f64 / seconds } else { 0.0 },
    });
    gc_log::info("FQ2BAM.FERTIG", "fq2bam beendet", result.clone());
    Ok(result)
}

/// Sortiert eine BAM nach Ort. Die Ausgabe entsteht unter einem vorläufigen Namen und bekommt
/// ihren Namen erst, wenn alles geschrieben ist; nach einem Abbruch liegt also nie eine halbe
/// Datei unter dem Zielnamen.
fn sort_bam(
    src: &std::path::Path,
    dst: &std::path::Path,
    speicher_mb: usize,
    temp: &std::path::Path,
    kerne: usize,
) -> Result<serde_json::Value> {
    let started = std::time::Instant::now();
    anyhow::ensure!(speicher_mb > 0, "--speicher-mb muss größer als 0 sein");
    if let (Ok(a), Ok(b)) = (std::fs::canonicalize(src), std::fs::canonicalize(dst)) {
        anyhow::ensure!(a != b, "Eingabe und Ausgabe sind dieselbe Datei");
    }
    let (mut reader, header, raw) =
        gc_io::bam::open_raw(src, kerne).with_context(|| format!("{} öffnen", src.display()))?;
    let end = raw.iter().position(|&b| b == 0).unwrap_or(raw.len());
    let text = std::str::from_utf8(&raw[..end]).context("Kopftext ist kein UTF-8")?;
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let new_text = gc_bam::sorted_header(
        text,
        "genome-sort",
        env!("CARGO_PKG_VERSION"),
        &command_line,
    );

    let mut sorter = gc_bam::Sorter::new(gc_bam::SortOptions {
        memory_bytes: speicher_mb.saturating_mul(1 << 20),
        temp_dir: temp.to_path_buf(),
    });
    let mut st = gc_log::stage("sort.lesen", json!({}));
    let mut record = Vec::new();
    let mut runs = 0;
    loop {
        record.clear();
        if !gc_io::bam::read_encoded(&mut reader, &mut record)
            .with_context(|| format!("{} lesen", src.display()))?
        {
            break;
        }
        sorter
            .push(&record)
            .with_context(|| format!("Teil in {} schreiben", temp.display()))?;
        let now = sorter.stats();
        if now.runs != runs {
            runs = now.runs;
            gc_log::info(
                "SORT.TEIL",
                "Speicher voll, sortierter Teil auf die Platte geschrieben",
                json!({ "teile": now.runs, "bytes_auf_platte": now.run_bytes, "saetze_bisher": now.records }),
            );
        }
        st.heartbeat(30, || json!({ "saetze": now.records, "teile": now.runs }));
    }
    let read = sorter.stats();
    st.finish(json!({ "saetze": read.records, "bytes": read.bytes, "teile": read.runs }));
    let seconds_read = started.elapsed().as_secs_f64();

    let mut part = dst.as_os_str().to_owned();
    part.push(".unfertig");
    let part = PathBuf::from(part);
    let st = gc_log::stage("sort.schreiben", json!({ "teile": read.runs }));
    let written = (|| -> Result<gc_bam::SortStats> {
        let mut writer = gc_io::bam::create_raw(&part, &header, new_text.as_bytes(), kerne)
            .with_context(|| format!("{} anlegen", part.display()))?;
        let stats = sorter
            .finish(|rec| gc_io::bam::write_encoded(&mut writer, rec))
            .context("Teile zusammenführen")?;
        gc_io::bam::finish(writer).with_context(|| format!("{} abschließen", part.display()))?;
        std::fs::rename(&part, dst).with_context(|| format!("{} umbenennen", part.display()))?;
        Ok(stats)
    })();
    let stats = match written {
        Ok(s) => s,
        Err(e) => {
            let _ = std::fs::remove_file(&part);
            return Err(e);
        }
    };
    st.finish(json!({ "saetze": stats.records }));
    let seconds = started.elapsed().as_secs_f64();
    let result = json!({
        "saetze": stats.records,
        "bytes_saetze": stats.bytes,
        "teile": stats.runs,
        "bytes_teile": stats.run_bytes,
        "speicher_mb": speicher_mb,
        "sekunden": seconds,
        "sekunden_lesen": seconds_read,
        "saetze_je_sekunde": if seconds > 0.0 { stats.records as f64 / seconds } else { 0.0 },
    });
    gc_log::info("SORT.FERTIG", "Sortieren beendet", result.clone());
    Ok(result)
}

/// Markiert Duplikate in einer nach Ort sortierten BAM. Die Ausgabe entsteht wie beim Sortieren
/// unter einem vorläufigen Namen.
fn markdup_bam(
    src: &std::path::Path,
    dst: &std::path::Path,
    max_length: i64,
    kerne: usize,
) -> Result<serde_json::Value> {
    let started = std::time::Instant::now();
    anyhow::ensure!(max_length > 0, "--max-laenge muss größer als 0 sein");
    if let (Ok(a), Ok(b)) = (std::fs::canonicalize(src), std::fs::canonicalize(dst)) {
        anyhow::ensure!(a != b, "Eingabe und Ausgabe sind dieselbe Datei");
    }
    let (mut reader, header, raw) =
        gc_io::bam::open_raw(src, kerne).with_context(|| format!("{} öffnen", src.display()))?;
    let end = raw.iter().position(|&b| b == 0).unwrap_or(raw.len());
    let text = std::str::from_utf8(&raw[..end]).context("Kopftext ist kein UTF-8")?;
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let new_text = gc_bam::with_program_line(
        text,
        "genome-markdup",
        env!("CARGO_PKG_VERSION"),
        &command_line,
    );
    let mut part = dst.as_os_str().to_owned();
    part.push(".unfertig");
    let part = PathBuf::from(part);
    let mut st = gc_log::stage("markdup.durchlauf", json!({ "max_laenge": max_length }));
    let done = (|| -> Result<gc_bam::DupStats> {
        let mut writer = gc_io::bam::create_raw(&part, &header, new_text.as_bytes(), kerne)
            .with_context(|| format!("{} anlegen", part.display()))?;
        let mut marker = gc_bam::DupMarker::new(gc_bam::DupOptions { max_length });
        let mut record = Vec::new();
        loop {
            record.clear();
            if !gc_io::bam::read_encoded(&mut reader, &mut record)
                .with_context(|| format!("{} lesen", src.display()))?
            {
                break;
            }
            marker.push(&record, &mut |rec| {
                gc_io::bam::write_encoded(&mut writer, rec)
            })?;
            st.heartbeat(30, || {
                let now = marker.stats();
                json!({ "saetze": now.read, "duplikate": now.duplicate_pair + now.duplicate_single })
            });
        }
        let stats = marker.finish(&mut |rec| gc_io::bam::write_encoded(&mut writer, rec))?;
        gc_io::bam::finish(writer).with_context(|| format!("{} abschließen", part.display()))?;
        std::fs::rename(&part, dst).with_context(|| format!("{} umbenennen", part.display()))?;
        Ok(stats)
    })();
    let stats = match done {
        Ok(s) => s,
        Err(e) => {
            let _ = std::fs::remove_file(&part);
            return Err(e);
        }
    };
    st.finish(json!({ "saetze": stats.read }));
    let seconds = started.elapsed().as_secs_f64();
    let duplicates = stats.duplicate_pair + stats.duplicate_single;
    let result = json!({
        "gelesen": stats.read,
        "geschrieben": stats.written,
        "nicht_bewertet": stats.excluded,
        "bewertet": stats.examined,
        "mit_partner": stats.paired,
        "ohne_partner": stats.single,
        "duplikate_mit_partner": stats.duplicate_pair,
        "duplikate_ohne_partner": stats.duplicate_single,
        "duplikate": duplicates,
        "anteil_duplikate": if stats.examined > 0 { duplicates as f64 / stats.examined as f64 } else { 0.0 },
        "geschaetzte_bibliotheksgroesse": stats.estimated_library_size(),
        "fenster_spitze_saetze": stats.window_peak,
        "max_laenge": max_length,
        "sekunden": seconds,
        "saetze_je_sekunde": if seconds > 0.0 { stats.read as f64 / seconds } else { 0.0 },
    });
    gc_log::info("MARKDUP.FERTIG", "Duplikate markiert", result.clone());
    Ok(result)
}

/// Schreibt den Index `.bai` zu einer nach Ort sortierten BAM.
fn index_bam(
    bam: &std::path::Path,
    dst: &std::path::Path,
    kerne: usize,
) -> Result<serde_json::Value> {
    let started = std::time::Instant::now();
    let (mut reader, header, _) =
        gc_io::bam::open_raw(bam, kerne).with_context(|| format!("{} öffnen", bam.display()))?;
    let n_ref = header.reference_sequences().len();
    // Stelle hinter dem zuletzt gelesenen Datensatz; am Anfang die Stelle hinter dem Kopf.
    let mut position = gc_io::bam::virtual_position(&reader);
    let mut builder = gc_bam::BaiBuilder::new(n_ref, position);
    let mut st = gc_log::stage("bamindex.lesen", json!({ "contigs": n_ref }));
    let mut record = Vec::new();
    let mut records = 0u64;
    loop {
        record.clear();
        if !gc_io::bam::read_encoded(&mut reader, &mut record)
            .with_context(|| format!("{} lesen", bam.display()))?
        {
            break;
        }
        position = gc_io::bam::virtual_position(&reader);
        builder.push_record(&record[4..], position)?;
        records += 1;
        st.heartbeat(30, || json!({ "saetze": records }));
    }
    st.finish(json!({ "saetze": records }));
    // Nicht die Stelle nach dem Dateiende: Dort steht noch der leere Schlussblock der BAM.
    let index = builder.finish(position);
    let mut bytes = Vec::new();
    index.write(&mut bytes)?;
    let mut part = dst.as_os_str().to_owned();
    part.push(".unfertig");
    let part = PathBuf::from(part);
    std::fs::write(&part, &bytes).with_context(|| format!("{} schreiben", part.display()))?;
    std::fs::rename(&part, dst).with_context(|| format!("{} umbenennen", part.display()))?;
    let seconds = started.elapsed().as_secs_f64();
    let result = json!({
        "saetze": records,
        "ohne_contig": index.unplaced,
        "contigs": n_ref,
        "contigs_mit_daten": index.references.iter().filter(|r| !r.bins.is_empty()).count(),
        "faecher": index.references.iter().map(|r| r.bins.len()).sum::<usize>(),
        "bytes_index": bytes.len(),
        "sekunden": seconds,
    });
    gc_log::info("BAMINDEX.FERTIG", "Index geschrieben", result.clone());
    Ok(result)
}

/// Verkettet die Seeds aller Reads und fasst zusammen.
fn index_ketten(
    fm: &gc_index::FmIndex,
    path: &std::path::Path,
    truth_in_name: bool,
    max_occ: u64,
    missing_path: Option<&std::path::Path>,
    table_path: Option<&std::path::Path>,
) -> Result<serde_json::Value> {
    use std::io::Write;
    let mut table = match table_path {
        Some(p) => Some(std::io::BufWriter::new(
            std::fs::File::create(p)
                .with_context(|| format!("{} nicht schreibbar", p.display()))?,
        )),
        None => None,
    };
    let final_opt = gc_align::FinalizeOptions::default();
    let mut mapq_hist = [0u64; 61];
    let mut missing = match missing_path {
        Some(p) => Some(std::io::BufWriter::new(
            std::fs::File::create(p)
                .with_context(|| format!("{} nicht schreibbar", p.display()))?,
        )),
        None => None,
    };
    let seed_opt = gc_align::SeedOptions {
        max_occ,
        ..gc_align::SeedOptions::default()
    };
    let chain_opt = gc_align::ChainOptions::default();
    let filter_opt = gc_align::FilterOptions::default();
    let region_opt = gc_align::RegionOptions::default();
    let (mut regions_total, mut region_right, mut region_tied, mut region_exact) =
        (0u64, 0u64, 0u64, 0u64);
    let mut regions = Vec::new();
    let l = fm.reference().len();
    let mut reader = fastq::open(path)?;
    let mut record = Record::default();
    let started = std::time::Instant::now();
    let (mut reads, mut seeds_total, mut chains_total, mut without) = (0u64, 0u64, 0u64, 0u64);
    let (mut best_right, mut best_tied, mut any_right, mut with_truth) = (0u64, 0u64, 0u64, 0u64);
    let mut seed_right = 0u64;
    let mut chains_filtered = 0u64;
    let (mut region_wrong, mut region_missed) = (0u64, 0u64);
    let (mut region_equal, mut region_better) = (0u64, 0u64);
    // Wo ging der wahre Ort verloren? Zeile 0: gleich gut, Zeile 1: besser gewesen.
    // Spalten: kein Seed, keine Kette, Kette gefiltert, Bereich mit zu wenig Punkten.
    let mut lost_at = [[0u64; 4]; 2];
    while reader.read_record(&mut record)? {
        let read = encode_read(&record.seq);
        let intervals = gc_align::collect_intervals(fm, &read, &seed_opt);
        let frac_rep = gc_align::repeat_fraction(&intervals, seed_opt.max_occ, read.len());
        let seeds = gc_align::locate_seeds(fm, &intervals, &seed_opt);
        let chains = gc_align::chain_seeds(&seeds, l, &chain_opt);
        let filtered = gc_align::filter_chains(chains.clone(), &filter_opt);
        regions.clear();
        for c in &filtered {
            gc_align::chain_to_regions(fm, &read, c, &region_opt, &mut regions);
        }
        regions_total += regions.len() as u64;
        // Bereinigen, primär bestimmen, MAPQ. Die Zählungen weiter unten sehen die bereinigte Liste.
        gc_align::sort_dedup(
            &mut regions,
            fm,
            Some(&read),
            &region_opt.scoring,
            &final_opt,
        );
        gc_align::mark_primary(&mut regions, reads, &region_opt.scoring, &final_opt);
        if let Some(r) = regions.first() {
            let q = gc_align::mapq(r, frac_rep, &region_opt.scoring, &final_opt);
            mapq_hist[q as usize] += 1;
            if let Some(out) = table.as_mut() {
                let cg = &fm.reference().contigs[r.contig as usize];
                // Ort, CIGAR, NM und MD kommen aus dem fertigen Alignment, der Rest aus dem Bereich.
                let aln = gc_align::region_to_alignment(fm, &read, r, &region_opt.scoring)
                    .context("Bereich ohne gültiges Alignment")?;
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{q}\t{}\t{}\t{}\t{}\t{}\t{}\t{frac_rep:.4}\t{}\t{}\t{}",
                    String::from_utf8_lossy(&record.name),
                    cg.name,
                    aln.pos + 1,
                    if aln.is_reverse { '-' } else { '+' },
                    r.score,
                    r.truesc,
                    r.qb,
                    r.qe,
                    r.sub,
                    r.sub_n,
                    gc_align::cigar_string(&aln.cigar),
                    aln.nm,
                    aln.md
                )?;
            }
        } else if let Some(out) = table.as_mut() {
            writeln!(
                out,
                "{}\t*\t0\t*\t0\t0\t0\t0\t0",
                String::from_utf8_lossy(&record.name)
            )?;
        }
        chains_filtered += filtered.len() as u64;
        reads += 1;
        seeds_total += seeds.len() as u64;
        chains_total += chains.len() as u64;
        without += u64::from(chains.is_empty());
        if !truth_in_name {
            continue;
        }
        let name = String::from_utf8_lossy(&record.name).into_owned();
        let mut parts = name.split(':');
        let (Some(contig), Some(pos), Some(strand)) = (parts.next(), parts.next(), parts.next())
        else {
            anyhow::bail!("Read-Name {name} hat nicht die Form contig:pos:strang");
        };
        let pos: u64 = pos
            .parse()
            .with_context(|| format!("Position im Read-Namen {name}"))?;
        let reverse = strand == "-";
        with_truth += 1;
        // Linkes Ende des Reads auf dem Vorwärtsstrang, wie es die Kette nahelegt.
        let at_truth = |c: &gc_align::Chain| -> bool {
            let s = c.seeds[0];
            let start_text = s.rbeg as i64 - i64::from(s.qbeg);
            let forward_start = if c.is_reverse(l) {
                2 * l as i64 - start_text - read.len() as i64
            } else {
                start_text
            };
            let cg = &fm.reference().contigs[c.contig as usize];
            cg.name == contig
                && c.is_reverse(l) == reverse
                && (forward_start - (cg.offset + pos - 1) as i64).abs() <= chain_opt.band_width
        };
        let top = chains
            .iter()
            .map(gc_align::Chain::weight)
            .max()
            .unwrap_or(0);
        let top_chains: Vec<_> = chains.iter().filter(|c| c.weight() == top).collect();
        let right_on_top = top_chains.iter().any(|c| at_truth(c));
        best_right += u64::from(right_on_top && top_chains.len() == 1);
        best_tied += u64::from(right_on_top && top_chains.len() > 1);
        any_right += u64::from(chains.iter().any(at_truth));
        // Bereiche: liegt der mit dem höchsten Score am wahren Ort?
        let region_start = |r: &gc_align::Region| -> (bool, i64) {
            let start_text = r.rb as i64 - i64::from(r.qb);
            if r.rb >= l {
                (true, 2 * l as i64 - start_text - read.len() as i64)
            } else {
                (false, start_text)
            }
        };
        let region_off = |r: &gc_align::Region| -> Option<i64> {
            let cg = &fm.reference().contigs[r.contig as usize];
            let (rev, start) = region_start(r);
            (cg.name == contig && rev == reverse).then(|| start - (cg.offset + pos - 1) as i64)
        };
        let top_score = regions.iter().map(|r| r.score).max().unwrap_or(-1);
        let top_regions: Vec<_> = regions.iter().filter(|r| r.score == top_score).collect();
        let hit = top_regions
            .iter()
            .any(|r| region_off(r).is_some_and(|d| d.abs() <= chain_opt.band_width));
        // Gleichstand zählt nur, wenn die gleich guten Bereiche an verschiedenen Orten liegen.
        let places: std::collections::BTreeSet<(u32, u64)> =
            top_regions.iter().map(|r| (r.contig, r.rb)).collect();
        let mut lost_row = None;
        if !hit {
            region_wrong += 1;
            // Score am wahren Ort ohne Lücken: der beste zusammenhängende Abschnitt auf der wahren
            // Diagonale (lokal, wie der Score eines Bereichs). Die erzeugten Reads enthalten nur
            // Austausche, also ist das dort das richtige Alignment. Mit Lücken könnte ein Score am
            // wahren Ort höchstens noch größer sein, der Wert bleibt deshalb eine untere Schranke.
            let sc = &region_opt.scoring;
            let cg = fm
                .reference()
                .contigs
                .iter()
                .find(|c| c.name == contig)
                .with_context(|| format!("Contig {contig} aus dem Read-Namen fehlt im Index"))?;
            let fwd = cg.offset + pos - 1;
            let text_start = if reverse {
                2 * l - fwd - read.len() as u64
            } else {
                fwd
            };
            let (mut run, mut truth_score) = (0i32, 0i32);
            for (i, &b) in read.iter().enumerate() {
                let t = fm.text_base(text_start + i as u64);
                run += if b > 3 {
                    -sc.ambiguous
                } else if b == t {
                    sc.matching
                } else {
                    -sc.mismatch
                };
                run = run.max(0);
                truth_score = truth_score.max(run);
            }
            region_equal += u64::from(top_score == truth_score);
            region_better += u64::from(top_score > truth_score);
            if top_score <= truth_score {
                lost_row = Some(usize::from(top_score < truth_score));
            }
            region_missed += u64::from(top_score < truth_score);
        }
        region_right += u64::from(hit && places.len() == 1);
        region_tied += u64::from(hit && places.len() > 1);
        region_exact +=
            u64::from(places.len() == 1 && top_regions.iter().any(|r| region_off(r) == Some(0)));
        // Unabhängig vom Verketten: liegt wenigstens ein Seed genau auf der wahren Diagonale?
        let seed_at_truth = seeds.iter().any(|s| {
            let one = gc_align::Chain {
                contig: s.contig,
                seeds: vec![*s],
            };
            let cg = &fm.reference().contigs[s.contig as usize];
            let start_text = s.rbeg as i64 - i64::from(s.qbeg);
            let forward_start = if one.is_reverse(l) {
                2 * l as i64 - start_text - read.len() as i64
            } else {
                start_text
            };
            cg.name == contig
                && one.is_reverse(l) == reverse
                && forward_start == (cg.offset + pos - 1) as i64
        });
        seed_right += u64::from(seed_at_truth);
        if let Some(row) = lost_row {
            let col = if !seed_at_truth {
                0
            } else if !chains.iter().any(at_truth) {
                1
            } else if !filtered.iter().any(at_truth) {
                2
            } else {
                3
            };
            lost_at[row][col] += 1;
        }
        if let (false, Some(out)) = (seed_at_truth, missing.as_mut()) {
            writeln!(out, "{name}\t{}", String::from_utf8_lossy(&record.seq))?;
        }
    }
    if let Some(mut out) = missing {
        out.flush()?;
    }
    if let Some(mut out) = table {
        out.flush()?;
    }
    let seconds = started.elapsed().as_secs_f64();
    let mut result = json!({
        "reads": reads,
        "seeds": seeds_total,
        "ketten": chains_total,
        "ketten_nach_filter": chains_filtered,
        "bereiche": regions_total,
        "reads_ohne_kette": without,
        "sekunden": seconds,
        "reads_je_sekunde": if seconds > 0.0 { reads as f64 / seconds } else { 0.0 },
        "kerne": 1,
        "max_orte": max_occ,
        "mapq_0": mapq_hist[0],
        "mapq_1_bis_29": mapq_hist[1..30].iter().sum::<u64>(),
        "mapq_30_bis_59": mapq_hist[30..60].iter().sum::<u64>(),
        "mapq_60": mapq_hist[60],
    });
    let verloren = |a: &[u64; 4]| {
        json!({
            "kein_seed_genau_am_wahren_ort": a[0],
            "keine_kette_am_wahren_ort": a[1],
            "kette_am_wahren_ort_gefiltert": a[2],
            "bereich_am_wahren_ort_mit_weniger_punkten": a[3],
        })
    };
    if truth_in_name {
        result["wahrheit"] = json!({
            "reads": with_truth,
            "schwerste_kette_am_wahren_ort_eindeutig": best_right,
            "schwerste_kette_am_wahren_ort_gleichstand": best_tied,
            "irgendeine_kette_am_wahren_ort": any_right,
            "bester_bereich_am_wahren_ort_eindeutig": region_right,
            "bester_bereich_am_wahren_ort_gleichstand": region_tied,
            "bester_bereich_eindeutig_und_genau_am_wahren_ort": region_exact,
            "bester_bereich_woanders": region_wrong,
            "davon_score_unter_dem_des_wahren_orts": region_missed,
            "davon_score_gleich_dem_des_wahren_orts": region_equal,
            "davon_score_ueber_dem_des_wahren_orts": region_better,
            "wahrer_ort_verloren_bei_gleichem_score": verloren(&lost_at[0]),
            "wahrer_ort_verloren_bei_besserem_score": verloren(&lost_at[1]),
            "irgendein_seed_genau_am_wahren_ort": seed_right,
            "toleranz_basen": chain_opt.band_width,
        });
    }
    Ok(result)
}

/// Read-Basen als Codes 0..3, alles andere als 4 (N).
fn encode_read(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| gc_index::reference::encode(b).unwrap_or(4))
        .collect()
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
        "standard" => gc_index::Variant::Standard,
        "kompakt" => gc_index::Variant::Kompakt,
        "schnell" => gc_index::Variant::Schnell,
        other => {
            anyhow::bail!("unbekannte Variante {other}: erlaubt sind standard, kompakt und schnell")
        }
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
        Command::Sort {
            src,
            dst,
            speicher_mb,
            temp,
            kerne,
            log,
        } => {
            gc_log::init(gc_log::Config {
                tool: "sort".to_string(),
                version: env!("CARGO_PKG_VERSION").to_string(),
                level: gc_log::Level::Info,
                jsonl: log,
                stderr: true,
            })?;
            gc_log::install_panic_hook();
            let temp = temp.unwrap_or_else(|| match dst.parent() {
                Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
                _ => PathBuf::from("."),
            });
            gc_log::run_header(json!({
                "eingabe": src.display().to_string(),
                "ausgabe": dst.display().to_string(),
                "speicher_mb": speicher_mb,
                "temp": temp.display().to_string(),
                "kerne": kerne,
            }));
            gc_log::input_file("eingabe", &src);
            match sort_bam(&src, &dst, speicher_mb, &temp, kerne) {
                Ok(mut result) => {
                    result["lauf"] = json!(gc_log::run_id());
                    gc_log::run_end(true, result.clone());
                    println!("{}", Envelope::new("genome.sort.v1", result).to_json()?);
                }
                Err(e) => {
                    gc_log::failure(gc_log::Failure {
                        code: "SORT.FEHLGESCHLAGEN",
                        what: "Das Sortieren wurde abgebrochen, es gibt keine Ausgabedatei",
                        location: json!({
                            "eingabe": src.display().to_string(),
                            "ausgabe": dst.display().to_string(),
                            "temp": temp.display().to_string(),
                        }),
                        cause: &format!("{e:#}"),
                        action: "Eingabe mit samtools quickcheck prüfen, freien Platz im Temp-Ordner und am Ziel prüfen, dann neu starten",
                    });
                    gc_log::run_end(false, json!({ "fehler": format!("{e:#}") }));
                    return Err(e);
                }
            }
        }
        Command::Markdup {
            src,
            dst,
            max_laenge,
            kerne,
            log,
        } => {
            gc_log::init(gc_log::Config {
                tool: "markdup".to_string(),
                version: env!("CARGO_PKG_VERSION").to_string(),
                level: gc_log::Level::Info,
                jsonl: log,
                stderr: true,
            })?;
            gc_log::install_panic_hook();
            gc_log::run_header(json!({
                "eingabe": src.display().to_string(),
                "ausgabe": dst.display().to_string(),
                "max_laenge": max_laenge,
                "kerne": kerne,
            }));
            gc_log::input_file("eingabe", &src);
            match markdup_bam(&src, &dst, max_laenge, kerne) {
                Ok(mut result) => {
                    result["lauf"] = json!(gc_log::run_id());
                    gc_log::run_end(true, result.clone());
                    println!("{}", Envelope::new("genome.markdup.v1", result).to_json()?);
                }
                Err(e) => {
                    gc_log::failure(gc_log::Failure {
                        code: "MARKDUP.FEHLGESCHLAGEN",
                        what: "Das Markieren der Duplikate wurde abgebrochen, es gibt keine Ausgabedatei",
                        location: json!({
                            "eingabe": src.display().to_string(),
                            "ausgabe": dst.display().to_string(),
                        }),
                        cause: &format!("{e:#}"),
                        action: "Die Eingabe muss nach Ort sortiert sein (genome sort) und die Partner-Angaben MC und ms tragen (genome index sam --partner-angaben oder samtools fixmate -m)",
                    });
                    gc_log::run_end(false, json!({ "fehler": format!("{e:#}") }));
                    return Err(e);
                }
            }
        }
        Command::BamIndex {
            bam,
            aus,
            kerne,
            log,
        } => {
            gc_log::init(gc_log::Config {
                tool: "bam-index".to_string(),
                version: env!("CARGO_PKG_VERSION").to_string(),
                level: gc_log::Level::Info,
                jsonl: log,
                stderr: true,
            })?;
            gc_log::install_panic_hook();
            let dst = aus.unwrap_or_else(|| {
                let mut name = bam.as_os_str().to_owned();
                name.push(".bai");
                PathBuf::from(name)
            });
            gc_log::run_header(json!({
                "bam": bam.display().to_string(),
                "index": dst.display().to_string(),
                "kerne": kerne,
            }));
            gc_log::input_file("bam", &bam);
            match index_bam(&bam, &dst, kerne) {
                Ok(mut result) => {
                    result["lauf"] = json!(gc_log::run_id());
                    gc_log::run_end(true, result.clone());
                    println!(
                        "{}",
                        Envelope::new("genome.bam-index.v1", result).to_json()?
                    );
                }
                Err(e) => {
                    gc_log::failure(gc_log::Failure {
                        code: "BAMINDEX.FEHLGESCHLAGEN",
                        what: "Der Index wurde nicht geschrieben",
                        location: json!({
                            "bam": bam.display().to_string(),
                            "index": dst.display().to_string(),
                        }),
                        cause: &format!("{e:#}"),
                        action: "Die BAM muss vollständig und nach Ort sortiert sein (genome sort); Contigs über 512 Mbp passen nicht in einen BAI-Index",
                    });
                    gc_log::run_end(false, json!({ "fehler": format!("{e:#}") }));
                    return Err(e);
                }
            }
        }
        Command::Flagstat { bam, kerne } => {
            let started = std::time::Instant::now();
            let (mut reader, _, _) = gc_io::bam::open_raw(&bam, kerne)
                .with_context(|| format!("{} öffnen", bam.display()))?;
            let mut stats = gc_bam::FlagStats::default();
            let mut record = Vec::new();
            loop {
                record.clear();
                if !gc_io::bam::read_encoded(&mut reader, &mut record)
                    .with_context(|| format!("{} lesen", bam.display()))?
                {
                    break;
                }
                anyhow::ensure!(
                    stats.add(&record[4..]),
                    "{}: Datensatz ist kürzer als seine festen Felder",
                    bam.display()
                );
            }
            let group = flag_counts_json;
            let result = json!({
                "bestanden": group(&stats.passed),
                "nicht_bestanden": group(&stats.failed),
                "sekunden": started.elapsed().as_secs_f64(),
            });
            println!("{}", Envelope::new("genome.flagstat.v1", result).to_json()?);
        }
        Command::Fq2bam {
            index,
            fastq,
            fastq2,
            aus,
            lesegruppe,
            speicher_mb,
            temp,
            kerne,
            block_basen,
            max_laenge,
            log,
            trimmen,
        } => {
            let threads = if kerne == 0 {
                std::thread::available_parallelism().map_or(1, usize::from)
            } else {
                kerne
            };
            let log = log.unwrap_or_else(|| {
                let mut name = aus.as_os_str().to_owned();
                name.push(".log.jsonl");
                PathBuf::from(name)
            });
            gc_log::init(gc_log::Config {
                tool: "fq2bam".to_string(),
                version: env!("CARGO_PKG_VERSION").to_string(),
                level: gc_log::Level::Info,
                jsonl: Some(log),
                stderr: true,
            })?;
            gc_log::install_panic_hook();
            let temp = temp.unwrap_or_else(|| match aus.parent() {
                Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
                _ => PathBuf::from("."),
            });
            gc_log::run_header(json!({
                "index": index.display().to_string(),
                "fastq": fastq.display().to_string(),
                "fastq2": fastq2.as_ref().map(|p| p.display().to_string()),
                "ausgabe": aus.display().to_string(),
                "lesegruppe": lesegruppe,
                "speicher_mb": speicher_mb,
                "temp": temp.display().to_string(),
                "kerne": threads,
                "block_basen": block_basen,
                "max_laenge": max_laenge,
                "trimmen": trimmen,
            }));
            gc_log::input_file("fastq", &fastq);
            if let Some(p) = &fastq2 {
                gc_log::input_file("fastq2", p);
            }
            gc_log::input_file("index", &index);
            let run = Fq2bamRun {
                index: &index,
                fastq: &fastq,
                fastq2: fastq2.as_deref(),
                dst: &aus,
                read_group: lesegruppe.as_deref(),
                memory_mb: speicher_mb,
                temp: &temp,
                threads,
                block_bases: block_basen,
                max_length: max_laenge,
                trim: trimmen,
            };
            match fq2bam(&run) {
                Ok(mut result) => {
                    result["lauf"] = json!(gc_log::run_id());
                    gc_log::run_end(true, result.clone());
                    println!("{}", Envelope::new("genome.fq2bam.v1", result).to_json()?);
                }
                Err(e) => {
                    gc_log::failure(gc_log::Failure {
                        code: "FQ2BAM.FEHLGESCHLAGEN",
                        what: "fq2bam wurde abgebrochen, es gibt keine Ausgabe-BAM",
                        location: json!({
                            "fastq": fastq.display().to_string(),
                            "ausgabe": aus.display().to_string(),
                            "temp": temp.display().to_string(),
                        }),
                        cause: &format!("{e:#}"),
                        action: "Ursache im Feld ursache lesen; FASTQ mit genome fastq-check prüfen, freien Platz im Temp-Ordner und am Ziel prüfen, dann neu starten",
                    });
                    gc_log::run_end(false, json!({ "fehler": format!("{e:#}") }));
                    return Err(e);
                }
            }
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
            IndexCommand::Ketten {
                index,
                fastq,
                wahrheit_im_namen,
                max_orte,
                fehlende,
                tabelle,
            } => {
                let fm = load_index(&index)?;
                let result = index_ketten(
                    &fm,
                    &fastq,
                    wahrheit_im_namen,
                    max_orte,
                    fehlende.as_deref(),
                    tabelle.as_deref(),
                )?;
                println!(
                    "{}",
                    Envelope::new("genome.index-ketten.v1", result).to_json()?
                );
            }
            IndexCommand::Sam {
                index,
                fastq,
                fastq2,
                block_basen,
                log,
                kerne,
                bam,
                lesegruppe,
                partner_angaben,
            } => {
                let threads = if kerne == 0 {
                    std::thread::available_parallelism().map_or(1, usize::from)
                } else {
                    kerne
                };
                gc_log::init(gc_log::Config {
                    tool: "index.sam".to_string(),
                    version: env!("CARGO_PKG_VERSION").to_string(),
                    level: gc_log::Level::Info,
                    jsonl: log,
                    stderr: true,
                })?;
                gc_log::install_panic_hook();
                gc_log::run_header(json!({
                    "index": index.display().to_string(),
                    "fastq": fastq.display().to_string(),
                    "fastq2": fastq2.as_ref().map(|p| p.display().to_string()),
                    "block_basen": block_basen,
                    "kerne": threads,
                    "bam": bam.as_ref().map(|p| p.display().to_string()),
                    "lesegruppe": lesegruppe,
                    "partner_angaben": partner_angaben,
                }));
                let read_group = match &lesegruppe {
                    Some(text) => Some(gc_align::ReadGroup::parse(text).map_err(|cause| {
                        gc_log::failure(gc_log::Failure {
                            code: "ALIGN.LESEGRUPPE_UNGUELTIG",
                            what: "Die Angabe bei --lesegruppe ist keine gültige @RG-Zeile",
                            location: json!({ "angabe": text }),
                            cause: &cause,
                            action: "Form: '@RG\\tID:kennung\\tSM:probe', Tabulatoren als \\t geschrieben, ID ist Pflicht",
                        });
                        anyhow::anyhow!("--lesegruppe: {cause}")
                    })?),
                    None => None,
                };
                let fm = load_index(&index)?;
                index_sam(
                    &fm,
                    &fastq,
                    fastq2.as_deref(),
                    &SamRun {
                        block_bases: block_basen,
                        threads,
                        bam: bam.as_deref(),
                        read_group: read_group.as_ref(),
                        mate_tags: partner_angaben,
                        trimmer: None,
                    },
                )?;
            }
            IndexCommand::Smem {
                index,
                sequenz,
                fastq,
                min_laenge,
                orte,
                pruefen,
            } => {
                let fm = load_index(&index)?;
                let result = match (sequenz, fastq) {
                    (Some(seq), None) => {
                        let read = encode_read(seq.as_bytes());
                        let seeds: Vec<_> = fm
                            .smems(&read, min_laenge)
                            .iter()
                            .map(|m| {
                                let iv = m.interval;
                                let places: Vec<_> = (iv.k..iv.k + iv.size.min(orte))
                                    .filter_map(|row| {
                                        fm.hit_at(fm.locate(row), u64::from(m.len()))
                                    })
                                    .map(|h| {
                                        json!({
                                            "contig": fm.reference().contigs[h.contig as usize].name,
                                            "pos": h.pos + 1,
                                            "strang": if h.reverse { "-" } else { "+" },
                                        })
                                    })
                                    .collect();
                                json!({
                                    "read_von": m.start + 1,
                                    "read_bis": m.end,
                                    "laenge": m.len(),
                                    "fundstellen_roh": iv.size,
                                    "orte": places,
                                })
                            })
                            .collect();
                        json!({ "min_laenge": min_laenge, "read_laenge": read.len(), "seeds": seeds })
                    }
                    (None, Some(path)) => {
                        let mut reader = fastq::open(&path)?;
                        let mut record = Record::default();
                        let started = std::time::Instant::now();
                        let (mut reads, mut bases, mut seeds, mut covered, mut without) =
                            (0u64, 0u64, 0u64, 0u64, 0u64);
                        let mut rows = 0u64;
                        let (mut checked, mut deviations) = (0u64, 0u64);
                        while reader.read_record(&mut record)? {
                            let read = encode_read(&record.seq);
                            let found = fm.smems(&read, min_laenge);
                            reads += 1;
                            bases += read.len() as u64;
                            seeds += found.len() as u64;
                            without += u64::from(found.is_empty());
                            covered += found.iter().map(|m| u64::from(m.len())).sum::<u64>();
                            rows += found.iter().map(|m| m.interval.size).sum::<u64>();
                            if pruefen {
                                for (i, m) in found.iter().enumerate() {
                                    let (a, b) = (m.start as usize, m.end as usize);
                                    let iv = m.interval;
                                    let mut bad = fm.interval(&read[a..b])
                                        != (iv.k, iv.k + iv.size)
                                        || iv.size == 0;
                                    bad |=
                                        a > 0 && read[a - 1] <= 3 && fm.count(&read[a - 1..b]) > 0;
                                    bad |= b < read.len()
                                        && read[b] <= 3
                                        && fm.count(&read[a..=b]) > 0;
                                    bad |= i > 0
                                        && (found[i - 1].start >= m.start
                                            || found[i - 1].end >= m.end);
                                    deviations += u64::from(bad);
                                }
                                checked += found.len() as u64;
                            }
                        }
                        let seconds = started.elapsed().as_secs_f64();
                        json!({
                            "min_laenge": min_laenge,
                            "reads": reads,
                            "basen": bases,
                            "seeds": seeds,
                            "seed_basen": covered,
                            "fundstellen_roh": rows,
                            "reads_ohne_seed": without,
                            "sekunden": seconds,
                            "reads_je_sekunde": if seconds > 0.0 { reads as f64 / seconds } else { 0.0 },
                            "kerne": 1,
                            "nachgeprueft": checked,
                            "abweichungen": deviations,
                        })
                    }
                    _ => anyhow::bail!("genau eines von --sequenz und --fastq angeben"),
                };
                println!(
                    "{}",
                    Envelope::new("genome.index-smem.v1", result).to_json()?
                );
            }
        },
    }
    Ok(())
}
