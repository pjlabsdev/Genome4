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
    },
}

/// Wohin die Zeilen gehen: SAM-Text auf die Standardausgabe oder BAM in eine Datei.
enum AlignOutput {
    Sam(std::io::BufWriter<std::io::Stdout>),
    Bam(Box<gc_io::bam::MtWriter>),
}

impl AlignOutput {
    fn open(fm: &gc_index::FmIndex, bam: Option<&std::path::Path>, threads: usize) -> Result<Self> {
        use std::io::Write;
        let command_line = std::env::args().collect::<Vec<_>>().join(" ");
        let header = gc_align::sam_header(fm, env!("CARGO_PKG_VERSION"), &command_line);
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
                gc_io::bam::write_encoded(writer, &bytes)?;
            }
        }
        Ok(block.iter().map(|l| l.len() as u64).sum())
    }

    fn finish(self) -> Result<()> {
        use std::io::Write;
        match self {
            AlignOutput::Sam(mut out) => out.flush()?,
            AlignOutput::Bam(writer) => gc_io::bam::finish(*writer)?,
        }
        Ok(())
    }
}

/// Richtet Reads aus, einzeln oder als Paare, Block für Block auf `threads` Kernen.
/// Bei Paaren wird die Insertgröße je Block geschätzt.
fn index_sam(
    fm: &gc_index::FmIndex,
    first: &std::path::Path,
    second: Option<&std::path::Path>,
    block_bases: u64,
    threads: usize,
    bam: Option<&std::path::Path>,
) -> Result<()> {
    let opt = gc_align::AlignOptions::default();
    let mut out = AlignOutput::open(fm, bam, threads)?;
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
    out.finish()?;
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
    Ok(())
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
                }));
                let fm = load_index(&index)?;
                index_sam(
                    &fm,
                    &fastq,
                    fastq2.as_deref(),
                    block_basen,
                    threads,
                    bam.as_deref(),
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
