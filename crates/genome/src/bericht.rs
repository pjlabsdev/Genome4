//! Überblick über eine BAM als Faktendatei (`genome.bamqc.v1`) und als HTML-Seite.
//!
//! Die Zahlen kommen aus einem Durchlauf über die Datensätze ([`gc_bam::BamQc`]), aus dem Kopf
//! der BAM und, wenn vorhanden, aus dem Lauf-Protokoll von `fq2bam`. Die HTML-Seite wird nur aus
//! der Faktendatei gesetzt; sie rechnet nichts Neues. Sie ist eine einzige Datei ohne Netz,
//! umschaltbar zwischen Deutsch und Englisch.

use std::fmt::Write as _;
use std::path::Path;

use gc_bam::{BamQc, HistSummary};
use serde_json::{Value, json};

/// Schema der Faktendatei.
pub const SCHEMA: &str = "genome.bamqc.v1";

/// Gruppe eines Contigs nach seinem Namen (GRCh38 mit und ohne `chr`).
pub fn contig_group(name: &str) -> &'static str {
    let n = name.strip_prefix("chr").unwrap_or(name);
    if matches!(n.parse::<u32>(), Ok(1..=22)) && !n.starts_with('0') {
        "autosom"
    } else if n == "X" || n == "Y" {
        "geschlecht"
    } else if n == "M" || n == "MT" {
        "mito"
    } else if name.ends_with("_decoy") {
        "decoy"
    } else if name.ends_with("_alt") {
        "alt"
    } else if name.starts_with("HLA-") {
        "hla"
    } else if n == "EBV" {
        "ebv"
    } else if name.contains("_random") || name.starts_with("chrUn") || name.starts_with("Un_") {
        "unplatziert"
    } else {
        "andere"
    }
}

fn share(part: u64, whole: u64) -> Value {
    if whole > 0 {
        json!(part as f64 / whole as f64)
    } else {
        Value::Null
    }
}

/// Nur die besetzten Fächer als `[[wert, anzahl], …]`.
fn sparse(hist: &[u64]) -> Value {
    Value::Array(
        hist.iter()
            .enumerate()
            .filter(|(_, c)| **c > 0)
            .map(|(v, c)| json!([v, c]))
            .collect(),
    )
}

fn summary_json(s: &HistSummary) -> Value {
    json!({ "n": s.n, "mittel": s.mean, "streuung": s.sd, "median": s.median, "p05": s.p05, "p95": s.p95 })
}

/// Felder einer Kopfzeile (`@RG`, `@PG`, `@HD`) als Objekt.
fn header_fields(line: &str) -> Value {
    let mut m = serde_json::Map::new();
    for f in line.split('\t').skip(1) {
        if let Some((k, v)) = f.split_once(':') {
            m.insert(k.to_string(), json!(v));
        }
    }
    Value::Object(m)
}

/// Liest die Ereignisse des letzten Laufs aus einem Lauf-Protokoll `genome.log.v1`.
pub fn read_log(path: &Path) -> Vec<Value> {
    let Ok(text) = std::fs::read_to_string(path) else {
        return Vec::new();
    };
    let events: Vec<Value> = text
        .lines()
        .filter_map(|l| serde_json::from_str::<Value>(l).ok())
        .collect();
    // Nur der letzte Lauf zählt: ab dem letzten LAUF.START.
    let start = events
        .iter()
        .rposition(|e| e["code"] == "LAUF.START")
        .unwrap_or(0);
    events[start..].to_vec()
}

fn last_event<'a>(events: &'a [Value], code: &str) -> Option<&'a Value> {
    events.iter().rev().find(|e| e["code"] == code)
}

/// Faktendatei aus Zählung, Kopf und Protokoll.
pub fn facts(
    bam: &Path,
    header_text: &str,
    qc: &BamQc,
    events: &[Value],
    flag_counts: impl Fn(&gc_bam::Counts) -> Value,
) -> Value {
    let size = |p: &Path| std::fs::metadata(p).map(|m| m.len()).ok();
    let mut bai = bam.as_os_str().to_owned();
    bai.push(".bai");
    let bai = std::path::PathBuf::from(bai);

    // Kopf.
    let mut hd = Value::Null;
    let (mut read_groups, mut programs, mut comments) = (Vec::new(), Vec::new(), 0u64);
    for line in header_text.lines() {
        if line.starts_with("@HD") {
            hd = header_fields(line);
        } else if line.starts_with("@RG") {
            read_groups.push(header_fields(line));
        } else if line.starts_with("@PG") {
            programs.push(header_fields(line));
        } else if line.starts_with("@CO") {
            comments += 1;
        }
    }
    let mut groups: Vec<(&str, u64, u64)> = Vec::new();
    for c in &qc.contigs {
        let g = contig_group(&c.name);
        match groups.iter_mut().find(|x| x.0 == g) {
            Some(x) => {
                x.1 += 1;
                x.2 += c.length;
            }
            None => groups.push((g, 1, c.length)),
        }
    }

    // Abdeckung.
    let auto_len: u64 = qc
        .contigs
        .iter()
        .filter(|c| c.autosome)
        .map(|c| c.length)
        .sum();
    let auto_bases: u64 = qc
        .contigs
        .iter()
        .filter(|c| c.autosome)
        .map(|c| c.bases)
        .sum();
    let auto_mean = if auto_len > 0 {
        auto_bases as f64 / auto_len as f64
    } else {
        0.0
    };
    let depth_sum = BamQc::summary(&qc.depth_autosomes);
    let at_least = |t: usize| -> Value {
        share(
            qc.depth_autosomes[t.min(qc.depth_autosomes.len() - 1)..]
                .iter()
                .sum(),
            depth_sum.n,
        )
    };
    let mean_of = |name: &str| -> Option<f64> {
        qc.contigs
            .iter()
            .find(|c| c.name == name || c.name == format!("chr{name}"))
            .filter(|c| c.length > 0)
            .map(|c| c.bases as f64 / c.length as f64)
    };
    let x = mean_of("X").map(|v| v / auto_mean.max(1e-9));
    let y = mean_of("Y").map(|v| v / auto_mean.max(1e-9));
    let sex = match (x, y) {
        (Some(x), Some(y)) if x < 0.75 && y > 0.2 => "XY",
        (Some(x), Some(y)) if x > 0.8 && y < 0.05 => "XX",
        (Some(_), Some(_)) => "unklar",
        _ => "unbekannt",
    };
    let contigs: Vec<Value> = qc
        .contigs
        .iter()
        .map(|c| {
            let g = contig_group(&c.name);
            let main = matches!(g, "autosom" | "geschlecht" | "mito");
            let mut v = json!({
                "name": c.name,
                "gruppe": g,
                "laenge": c.length,
                "reads": c.reads,
                "duplikate": c.duplicates,
                "basen": c.bases,
                "mittel": if c.length > 0 { c.bases as f64 / c.length as f64 } else { 0.0 },
                "stellen_ab": c.covered,
                "mapq_mittel": if c.mapq_reads > 0 { json!(c.mapq_sum as f64 / c.mapq_reads as f64) } else { Value::Null },
            });
            if main {
                v["faecher"] = json!(c.bins);
                v["faecher_reads"] = json!(c.bin_reads);
                v["faecher_duplikate"] = json!(c.bin_dups);
                v["faecher_mapq0"] = json!(c.bin_mapq0);
                v["verteilung"] = sparse(&c.depth_hist);
            }
            v
        })
        .collect();

    // Tags.
    let mut tags: Vec<(String, u64)> = qc
        .tags
        .iter()
        .enumerate()
        .filter(|(_, n)| **n > 0)
        .map(|(i, n)| {
            let t = BamQc::tag_name(i);
            (String::from_utf8_lossy(&t).into_owned(), *n)
        })
        .collect();
    tags.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    // Lauf aus dem Protokoll.
    let lauf = match last_event(events, "FQ2BAM.FERTIG") {
        Some(done) => {
            let start = last_event(events, "LAUF.START");
            let adapters: Vec<Value> = events
                .iter()
                .filter(|e| e["code"] == "TRIM.ADAPTER")
                .map(|e| e["fields"].clone())
                .collect();
            json!({
                "werkzeug": done["tool"],
                "version": done["version"],
                "lauf": done["run"],
                "ende": done["time"],
                "beginn": start.map(|s| s["time"].clone()),
                "befehl": start.map(|s| s["fields"]["befehl"].clone()),
                "parameter": start.map(|s| s["fields"]["parameter"].clone()),
                "system": start.map(|s| json!({
                    "system": s["fields"]["system"],
                    "architektur": s["fields"]["architektur"],
                    "kerne": s["fields"]["kerne"],
                })),
                "ergebnis": done["fields"],
                "adapter": adapters,
                "polyg": last_event(events, "TRIM.POLYG").map(|e| e["fields"].clone()),
                "speicherplan": last_event(events, "FQ2BAM.SPEICHERPLAN").map(|e| e["fields"].clone()),
                "ausrichten": last_event(events, "ALIGN.FERTIG").map(|e| e["fields"].clone()),
            })
        }
        None => Value::Null,
    };

    let mapq_total: u64 = qc.mapq.iter().sum();
    json!({
        "schema": SCHEMA,
        "werkzeug": "genome",
        "version": env!("CARGO_PKG_VERSION"),
        "datei": {
            "bam": bam.display().to_string(),
            "bytes": size(bam),
            "bai": if bai.exists() { json!(bai.display().to_string()) } else { Value::Null },
            "bai_bytes": size(&bai),
        },
        "kopf": {
            "hd": hd,
            "contigs": qc.contigs.len(),
            "laenge": qc.contigs.iter().map(|c| c.length).sum::<u64>(),
            "gruppen": groups.iter().map(|(g, n, l)| json!({ "gruppe": g, "contigs": n, "laenge": l })).collect::<Vec<_>>(),
            "lesegruppen": read_groups,
            "programme": programs,
            "kommentare": comments,
        },
        "saetze": {
            "gesamt": qc.records,
            "bestanden": flag_counts(&qc.flags.passed),
            "nicht_bestanden": flag_counts(&qc.flags.failed),
            "hauptzeilen_ohne_ort": qc.unmapped,
        },
        "tags": tags.iter().map(|(t, n)| json!({ "tag": t, "saetze": n })).collect::<Vec<_>>(),
        "mapq": {
            "verteilung": sparse(&qc.mapq),
            "reads": mapq_total,
            "anteil_ab_20": share(qc.mapq[20..].iter().sum(), mapq_total),
            "anteil_0": share(qc.mapq[0], mapq_total),
            "anteil_60": share(qc.mapq[60..].iter().sum(), mapq_total),
        },
        "insertgroesse": {
            "kennzahlen": summary_json(&BamQc::summary(&qc.insert[..gc_bam::qc::MAX_INSERT])),
            "ueber_max": qc.insert[gc_bam::qc::MAX_INSERT],
            "max": gc_bam::qc::MAX_INSERT,
            "verteilung": sparse(&qc.insert),
        },
        "read_laenge": {
            "kennzahlen": summary_json(&BamQc::summary(&qc.read_length)),
            "verteilung": sparse(&qc.read_length),
        },
        "qualitaet": {
            "basen": qc.quality_bases,
            "mittel": if qc.quality_bases > 0 { json!(qc.quality_sum as f64 / qc.quality_bases as f64) } else { Value::Null },
            "anteil_q30": share(qc.q30_bases, qc.quality_bases),
        },
        "fehler": {
            "ausgerichtete_basen": qc.aligned_bases,
            "fehlpaarungen": qc.mismatches,
            "rate": share(qc.mismatches, qc.nm_aligned_bases),
            "nm_summe": qc.nm_sum,
            "rate_wie_samtools": share(qc.nm_sum, qc.nm_query_bases),
            "saetze_mit_nm": qc.nm_records,
            "insertierte_basen": qc.inserted_bases,
            "deletierte_basen": qc.deleted_bases,
            "soft_clip_basen": qc.soft_clipped_bases,
        },
        "abdeckung": {
            "autosomen": {
                "laenge": auto_len,
                "mittel": auto_mean,
                "median": depth_sum.median,
                "p05": depth_sum.p05,
                "p95": depth_sum.p95,
                "anteil_ab": { "1": at_least(1), "10": at_least(10), "20": at_least(20), "30": at_least(30) },
                "verteilung": sparse(&qc.depth_autosomes),
                "max": gc_bam::qc::MAX_DEPTH,
            },
            "fach_basen": gc_bam::qc::BIN,
            "schwellen": gc_bam::qc::DEPTH_STEPS,
            "geschlecht": { "x_zu_autosomen": x, "y_zu_autosomen": y, "deutung": sex },
            "gekuerzte_stuecke": qc.clipped_blocks,
            "contigs": contigs,
        },
        "lauf": lauf,
    })
}

// ---------------------------------------------------------------------------------------------
// HTML

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Zwei Sprachen nebeneinander; das Umschalten blendet eine aus.
fn l(de: &str, en: &str) -> String {
    if de == en {
        return de.to_string();
    }
    format!("<span class=\"de\">{de}</span><span class=\"en\">{en}</span>")
}

fn group_digits(int: &str, sep: char) -> String {
    let (sign, digits) = match int.strip_prefix('-') {
        Some(d) => ("-", d),
        None => ("", int),
    };
    let mut out = String::new();
    for (i, ch) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(sep);
        }
        out.push(ch);
    }
    format!("{sign}{out}")
}

/// Zahl mit `d` Nachkommastellen in beiden Schreibweisen.
fn num_pair(x: f64, d: usize) -> (String, String) {
    if !x.is_finite() {
        return ("–".into(), "–".into());
    }
    let s = format!("{x:.d$}");
    let (int, frac) = match s.split_once('.') {
        Some((i, f)) => (i.to_string(), Some(f.to_string())),
        None => (s.clone(), None),
    };
    let de = group_digits(&int, '.') + &frac.as_ref().map(|f| format!(",{f}")).unwrap_or_default();
    let en = group_digits(&int, ',') + &frac.map(|f| format!(".{f}")).unwrap_or_default();
    (de, en)
}

fn num(x: f64, d: usize) -> String {
    let (de, en) = num_pair(x, d);
    l(&de, &en)
}

fn int(v: &Value) -> String {
    match v.as_f64() {
        Some(x) => num(x, 0),
        None => "–".into(),
    }
}

fn pct(v: &Value, d: usize) -> String {
    match v.as_f64() {
        Some(x) => {
            let (de, en) = num_pair(100.0 * x, d);
            l(&format!("{de} %"), &format!("{en} %"))
        }
        None => "–".into(),
    }
}

fn bytes(v: &Value) -> String {
    let Some(b) = v.as_f64() else {
        return "–".into();
    };
    let (unit, div) = if b >= 1e9 {
        ("GB", 1e9)
    } else if b >= 1e6 {
        ("MB", 1e6)
    } else if b >= 1e3 {
        ("kB", 1e3)
    } else {
        ("B", 1.0)
    };
    let (de, en) = num_pair(b / div, if div > 1.0 { 1 } else { 0 });
    l(&format!("{de} {unit}"), &format!("{en} {unit}"))
}

fn duration(secs: f64) -> String {
    let s = secs.round() as u64;
    let (h, m, sec) = (s / 3600, (s % 3600) / 60, s % 60);
    if h > 0 {
        format!("{h} h {m} min")
    } else if m > 0 {
        format!("{m} min {sec} s")
    } else {
        let (de, en) = num_pair(secs, 1);
        l(&format!("{de} s"), &format!("{en} s"))
    }
}

/// `[[wert, anzahl], …]` zurück in ein Feld bis `max`.
fn dense(v: &Value, max: usize) -> Vec<u64> {
    let mut out = vec![0u64; max + 1];
    for p in v.as_array().into_iter().flatten() {
        let (Some(k), Some(n)) = (p[0].as_u64(), p[1].as_u64()) else {
            continue;
        };
        if (k as usize) <= max {
            out[k as usize] += n;
        }
    }
    out
}

fn tile(value: &str, de: &str, en: &str, note: &str) -> String {
    format!(
        "<div class=\"tile\"><div class=\"v\">{value}</div><div class=\"k\">{}</div><div class=\"n\">{note}</div></div>",
        l(de, en)
    )
}

fn row(cells: &[String]) -> String {
    let mut s = String::from("<tr>");
    for (i, c) in cells.iter().enumerate() {
        if i == 0 {
            let _ = write!(s, "<td>{c}</td>");
        } else {
            let _ = write!(s, "<td class=\"r\">{c}</td>");
        }
    }
    s.push_str("</tr>");
    s
}

fn tag_meaning(tag: &str) -> (&'static str, &'static str) {
    match tag {
        "NM" => (
            "Editierabstand zur Referenz (Fehlpaarungen plus eingefügte und fehlende Basen)",
            "Edit distance to the reference (mismatches plus inserted and deleted bases)",
        ),
        "MD" => (
            "Fehlpaarungen und Deletionen als Text, damit die Referenz ohne FASTA rekonstruierbar ist",
            "Mismatches and deletions as text, so the reference can be rebuilt without the FASTA",
        ),
        "AS" => ("Punktzahl des Alignments", "Alignment score"),
        "XS" => (
            "Punktzahl des besten anderen Orts (Hinweis auf Wiederholungen)",
            "Score of the best other location (hint of repeats)",
        ),
        "MC" => ("CIGAR des Partners", "CIGAR of the mate"),
        "MQ" => ("MAPQ des Partners", "MAPQ of the mate"),
        "ms" => (
            "Summe der Basenqualitäten des Partners (für das Markieren der Duplikate)",
            "Sum of the mate's base qualities (used for duplicate marking)",
        ),
        "RG" => ("Lesegruppe", "Read group"),
        "SA" => (
            "weitere Teile eines geteilten Alignments (ergänzende Alignments)",
            "Other parts of a split alignment (supplementary alignments)",
        ),
        "XA" => (
            "andere gleich gute oder fast gleich gute Orte",
            "Alternative equally or nearly equally good locations",
        ),
        "pa" => ("Anteil der Punktzahl im Paar", "Share of the pair score"),
        "MI" => ("Kennung des Moleküls (UMI)", "Molecule identifier (UMI)"),
        "BC" => ("Barcode", "Barcode"),
        "OQ" => ("ursprüngliche Basenqualitäten", "Original base qualities"),
        _ => ("", ""),
    }
}

fn group_name(g: &str) -> (&'static str, &'static str) {
    match g {
        "autosom" => ("Autosomen (chr1–22)", "Autosomes (chr1–22)"),
        "geschlecht" => ("Geschlechtschromosomen (X, Y)", "Sex chromosomes (X, Y)"),
        "mito" => ("Mitochondrium (chrM)", "Mitochondrion (chrM)"),
        "decoy" => ("Decoy-Contigs", "Decoy contigs"),
        "alt" => ("ALT-Contigs", "ALT contigs"),
        "hla" => ("HLA-Contigs", "HLA contigs"),
        "ebv" => ("Epstein-Barr-Virus (chrEBV)", "Epstein–Barr virus (chrEBV)"),
        "unplatziert" => (
            "nicht zugeordnete Stücke (chrUn, _random)",
            "Unplaced pieces (chrUn, _random)",
        ),
        _ => ("andere", "Other"),
    }
}

/// Summen einer Art von Contigs: Länge, Reads, Duplikate, Basen, Stellen ab Schwelle, Contigs.
type GroupSum = (u64, u64, u64, u64, [u64; 4], u64);

const CSS: &str = include_str!("bericht/seite.css");
const JS: &str = include_str!("bericht/seite.js");

fn round_to(x: f64, d: i32) -> f64 {
    let f = 10f64.powi(d);
    (x * f).round() / f
}

/// Daten der Grafiken, aus der Faktendatei verdichtet: je Chromosom die Fächer als mittlere
/// Abdeckung und Anteile, Verteilungen als dichte Felder bis zu einem sinnvollen Rand.
fn viz_data(f: &Value) -> Value {
    let cov = &f["abdeckung"]["autosomen"];
    let auto = cov["mittel"].as_f64().unwrap_or(0.0);
    let bin = f["abdeckung"]["fach_basen"]
        .as_f64()
        .unwrap_or(gc_bam::qc::BIN as f64);
    let depth_cut = ((auto * 2.5).ceil() as usize).clamp(10, gc_bam::qc::MAX_DEPTH);
    let mut chromosomes = Vec::new();
    for c in f["abdeckung"]["contigs"].as_array().into_iter().flatten() {
        if !matches!(
            c["gruppe"].as_str(),
            Some("autosom" | "geschlecht" | "mito")
        ) {
            continue;
        }
        let len = c["laenge"].as_f64().unwrap_or(0.0);
        let arr = |k: &str| -> Vec<f64> {
            c[k].as_array()
                .into_iter()
                .flatten()
                .map(|v| v.as_f64().unwrap_or(0.0))
                .collect()
        };
        let (sums, reads, dups, mq0) = (
            arr("faecher"),
            arr("faecher_reads"),
            arr("faecher_duplikate"),
            arr("faecher_mapq0"),
        );
        let depth: Vec<f64> = sums
            .iter()
            .enumerate()
            .map(|(i, s)| {
                let width = (len - i as f64 * bin).min(bin).max(1.0);
                round_to(s / width, 2)
            })
            .collect();
        let ratio = |a: &[f64], b: &[f64], d: i32| -> Vec<f64> {
            a.iter()
                .zip(b)
                .map(|(x, n)| if *n > 0.0 { round_to(x / n, d) } else { 0.0 })
                .collect()
        };
        let non_dup: Vec<f64> = reads.iter().zip(&dups).map(|(r, d)| r - d).collect();
        let mean = c["mittel"].as_f64().unwrap_or(0.0);
        let cut = ((auto.max(mean) * 2.5).ceil() as usize).clamp(10, gc_bam::qc::MAX_DEPTH);
        let hist = dense(&c["verteilung"], gc_bam::qc::MAX_DEPTH);
        let len_u = c["laenge"].as_u64().unwrap_or(0);
        let ab: Vec<Value> = (0..4)
            .map(|i| share(c["stellen_ab"][i].as_u64().unwrap_or(0), len_u))
            .collect();
        chromosomes.push(json!({
            "name": c["name"],
            "laenge": len_u,
            "mittel": round_to(mean, 3),
            "reads": c["reads"],
            "duplikate": c["duplikate"],
            "mapq": c["mapq_mittel"],
            "ab": ab,
            "abdeckung": depth,
            "fach_reads": reads,
            "mapq0": ratio(&mq0, &non_dup, 4),
            "dup": ratio(&dups, &reads, 4),
            "verteilung": if hist.iter().any(|v| *v > 0) { json!(hist[..=cut]) } else { Value::Null },
        }));
    }
    let mapq = dense(&f["mapq"]["verteilung"], 255);
    let range = |a: usize, b: usize| mapq[a..=b].iter().sum::<u64>();
    let mapq_bins: Vec<Value> = [
        ("0", 0, 0),
        ("1–9", 1, 9),
        ("10–19", 10, 19),
        ("20–29", 20, 29),
        ("30–39", 30, 39),
        ("40–49", 40, 49),
        ("50–59", 50, 59),
        ("≥ 60", 60, 255),
    ]
    .iter()
    .map(|(label, a, b)| json!({ "label": label, "n": range(*a, *b) }))
    .collect();
    let ins = &f["insertgroesse"]["kennzahlen"];
    let insert = if ins["n"].as_u64().unwrap_or(0) > 0 {
        let hist = dense(&f["insertgroesse"]["verteilung"], gc_bam::qc::MAX_INSERT);
        let p95 = ins["p95"].as_f64().unwrap_or(1000.0);
        let cut = ((p95 * 1.4) as usize).clamp(50, gc_bam::qc::MAX_INSERT - 1);
        json!({ "verteilung": hist[..=cut], "median": ins["median"] })
    } else {
        Value::Null
    };
    let lengths = dense(&f["read_laenge"]["verteilung"], gc_bam::qc::MAX_READ_LENGTH);
    let read_length = if lengths.iter().filter(|c| **c > 0).count() > 1 {
        let maxlen = lengths.iter().rposition(|c| *c > 0).unwrap_or(1);
        json!(lengths[..=maxlen])
    } else {
        Value::Null
    };
    let depth = dense(&cov["verteilung"], gc_bam::qc::MAX_DEPTH);
    let total: u64 = depth.iter().sum();
    let beyond: u64 = depth[depth_cut + 1..].iter().sum();
    json!({
        "autosomen_mittel": auto,
        "fach": bin,
        "zusammenfassen": ((1_000_000.0 / bin).round() as usize).max(1),
        "chromosomen": chromosomes,
        "mapq": mapq_bins,
        "insert": insert,
        "read_laenge": read_length,
        "tiefe": {
            "verteilung": depth[..=depth_cut],
            "median": cov["median"],
            "ab_ende": beyond,
            "stellen": total,
        },
    })
}

/// Setzt die HTML-Seite aus der Faktendatei.
pub fn html(f: &Value) -> String {
    let mut s = String::new();
    let bam_path = f["datei"]["bam"].as_str().unwrap_or("");
    let file_name = Path::new(bam_path)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let rg = f["kopf"]["lesegruppen"]
        .as_array()
        .and_then(|a| a.first())
        .cloned()
        .unwrap_or(Value::Null);
    let sample = rg["SM"].as_str().unwrap_or("");
    let passed = &f["saetze"]["bestanden"];
    let cov = &f["abdeckung"]["autosomen"];
    let lauf = &f["lauf"];
    let viz = viz_data(f);

    let _ = write!(
        s,
        "<!DOCTYPE html><html lang=\"de\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>{} – Genome</title><style>{CSS}</style></head><body><main>",
        esc(&file_name)
    );
    let _ = write!(
        s,
        "<header><div><h1>{}{}</h1><div class=\"sub\">{} · <code>{}</code></div></div><div class=\"lang\"><button data-l=\"de\" onclick=\"setLang('de')\">Deutsch</button><button data-l=\"en\" onclick=\"setLang('en')\">English</button></div></header>",
        l("Überblick über die BAM", "BAM overview"),
        if sample.is_empty() {
            String::new()
        } else {
            format!(" · {}", esc(sample))
        },
        l(
            &format!(
                "erstellt von Genome {}",
                esc(f["version"].as_str().unwrap_or(""))
            ),
            &format!(
                "made by Genome {}",
                esc(f["version"].as_str().unwrap_or(""))
            )
        ),
        esc(bam_path)
    );

    // Kernaussagen.
    s.push_str("<div class=\"tiles\">");
    let primary = passed["hauptzeilen"].as_f64().unwrap_or(0.0);
    s.push_str(&tile(
        &int(&passed["hauptzeilen"]),
        "Reads (Hauptzeilen)",
        "Reads (primary)",
        &if primary > 0.0 {
            let (de, en) = num_pair(primary / 2e6, 1);
            l(&format!("{de} Mio. Paare"), &format!("{en} million pairs"))
        } else {
            String::new()
        },
    ));
    s.push_str(&tile(
        &pct(&passed["anteil_hauptzeilen_mit_ort"], 2),
        "mit Ort im Genom",
        "mapped",
        &l("Hauptzeilen", "primary reads"),
    ));
    s.push_str(&tile(
        &pct(&passed["anteil_richtiges_paar"], 2),
        "richtige Paare",
        "properly paired",
        "",
    ));
    let dup_share = if primary > 0.0 {
        json!(passed["duplikate_hauptzeilen"].as_f64().unwrap_or(0.0) / primary)
    } else {
        Value::Null
    };
    s.push_str(&tile(
        &pct(&dup_share, 2),
        "Duplikate",
        "duplicates",
        &l(
            "von allen Hauptzeilen (Flag 1024)",
            "of all primary reads (flag 1024)",
        ),
    ));
    s.push_str(&tile(
        &format!("{}x", num(cov["mittel"].as_f64().unwrap_or(0.0), 1)),
        "mittlere Abdeckung",
        "mean coverage",
        &l("Autosomen, ohne Duplikate", "autosomes, without duplicates"),
    ));
    s.push_str(&tile(
        &pct(&cov["anteil_ab"]["20"], 1),
        "Stellen mit ≥ 20x",
        "positions with ≥ 20x",
        &l(
            "Autosomen, N-Strecken eingerechnet",
            "autosomes, N gaps included",
        ),
    ));
    s.push_str(&tile(
        &int(&f["insertgroesse"]["kennzahlen"]["median"]),
        "Insertgröße (Median)",
        "insert size (median)",
        &l("Basen", "bases"),
    ));
    s.push_str(&tile(
        &pct(&f["fehler"]["rate"], 3),
        "Fehlpaarungen",
        "mismatch rate",
        &l("je ausgerichtete Base", "per aligned base"),
    ));
    let sex = f["abdeckung"]["geschlecht"]["deutung"]
        .as_str()
        .unwrap_or("");
    s.push_str(&tile(
        &match sex {
            "XY" | "XX" => sex.to_string(),
            _ => l("unklar", "unclear"),
        },
        "Geschlechtschromosomen",
        "sex chromosomes",
        &l("aus der Abdeckung von X und Y", "from X and Y coverage"),
    ));
    s.push_str("</div>");

    // Genom auf einen Blick und Chromosomen im Detail.
    let _ = write!(
        s,
        "<h2>{}</h2><p class=\"note\">{}</p><div class=\"viz\" id=\"v-karyo\"></div><div class=\"legend\"><span>{}</span><span class=\"ramp\"></span><span>{}</span><span><i class=\"sw\" style=\"background:var(--nodata)\"></i>{}</span></div>",
        l("Das Genom auf einen Blick", "The genome at a glance"),
        l(
            "Jede Leiste ist ein Chromosom, jedes Feld 1 Mio. Basen. Die Farbe zeigt die Abdeckung im Verhältnis zum Mittel der Autosomen. Rechts die mittlere Abdeckung. Ein Klick öffnet das Chromosom im Detail.",
            "Each bar is a chromosome, each cell 1 million bases. Colour shows coverage relative to the autosome mean. Mean coverage on the right. Click to open the chromosome in detail."
        ),
        l("weniger (0x)", "less (0x)"),
        l("mehr (2× Mittel)", "more (2× mean)"),
        l("keine Reads (N-Strecken)", "no reads (N gaps)")
    );
    let _ = write!(
        s,
        "<h3>{}</h3><div class=\"viz\" id=\"v-rel\"></div><p class=\"note\">{}</p>",
        l(
            "Abdeckung je Chromosom im Verhältnis zu den Autosomen",
            "Coverage per chromosome relative to the autosomes"
        ),
        l(
            "1,0 entspricht zwei Kopien. Bei XY liegen X und Y nahe 0,5 (eine Kopie); Y liegt oft darunter, weil große Teile von Y nicht eindeutig kartierbar sind. chrM ist abgeschnitten (viele Mitochondrien je Zelle).",
            "1.0 means two copies. In XY, X and Y are near 0.5 (one copy); Y is often lower because much of Y cannot be mapped uniquely. chrM is clipped (many mitochondria per cell)."
        )
    );
    let _ = write!(
        s,
        "<h2 id=\"chromosomen\">{}</h2><div class=\"chips\" id=\"v-chips\"></div><div class=\"detail\" id=\"v-detail\"></div>",
        l("Chromosomen im Detail", "Chromosomes in detail")
    );

    // Datei und Kopf.
    let _ = write!(
        s,
        "<h2>{}</h2>",
        l("Was in der Datei steckt", "What the file contains")
    );
    let mut t = String::from("<table>");
    t.push_str(&row(&[l("BAM", "BAM"), bytes(&f["datei"]["bytes"])]));
    t.push_str(&row(&[
        l("Index (.bai)", "Index (.bai)"),
        bytes(&f["datei"]["bai_bytes"]),
    ]));
    let hd = &f["kopf"]["hd"];
    t.push_str(&row(&[
        l("Sortierung (@HD SO)", "Sort order (@HD SO)"),
        esc(hd["SO"].as_str().unwrap_or("–")),
    ]));
    t.push_str(&row(&[
        l("SAM-Version (@HD VN)", "SAM version (@HD VN)"),
        esc(hd["VN"].as_str().unwrap_or("–")),
    ]));
    t.push_str(&row(&[
        l("Contigs der Referenz (@SQ)", "Reference contigs (@SQ)"),
        format!("{} · {}", int(&f["kopf"]["contigs"]), {
            let (de, en) = num_pair(f["kopf"]["laenge"].as_f64().unwrap_or(0.0) / 1e9, 3);
            l(&format!("{de} Mrd. Basen"), &format!("{en} billion bases"))
        }),
    ]));
    t.push_str(&row(&[
        l("Datensätze insgesamt", "Records in total"),
        int(&f["saetze"]["gesamt"]),
    ]));
    t.push_str("</table>");
    s.push_str(&t);

    let _ = write!(
        s,
        "<h3>{}</h3><table><tr><th>{}</th><th class=\"r\">{}</th><th class=\"r\">{}</th></tr>",
        l("Referenz nach Art der Contigs", "Reference by contig type"),
        l("Art", "Type"),
        l("Contigs", "Contigs"),
        l("Basen", "Bases")
    );
    for g in f["kopf"]["gruppen"].as_array().into_iter().flatten() {
        let (de, en) = group_name(g["gruppe"].as_str().unwrap_or(""));
        s.push_str(&row(&[l(de, en), int(&g["contigs"]), int(&g["laenge"])]));
    }
    s.push_str("</table>");

    if !rg.is_null() {
        let _ = write!(
            s,
            "<h3>{}</h3><table>",
            l("Lesegruppe (@RG)", "Read group (@RG)")
        );
        for rg in f["kopf"]["lesegruppen"].as_array().into_iter().flatten() {
            for (k, v) in rg.as_object().into_iter().flatten() {
                let (de, en) = match k.as_str() {
                    "ID" => ("Kennung", "Identifier"),
                    "SM" => ("Probe", "Sample"),
                    "PL" => ("Plattform", "Platform"),
                    "LB" => ("Bibliothek", "Library"),
                    "PU" => ("Lauf/Spur", "Platform unit"),
                    "CN" => ("Labor", "Sequencing centre"),
                    "DT" => ("Datum", "Date"),
                    _ => ("", ""),
                };
                let label = if de.is_empty() {
                    k.clone()
                } else {
                    format!("{k} · {}", l(de, en))
                };
                s.push_str(&row(&[
                    label,
                    format!("<code>{}</code>", esc(v.as_str().unwrap_or(""))),
                ]));
            }
        }
        s.push_str("</table>");
    }
    let programs = f["kopf"]["programme"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    if !programs.is_empty() {
        let _ = write!(
            s,
            "<h3>{}</h3><table><tr><th>{}</th><th>Version</th><th>{}</th></tr>",
            l(
                "Programme, die die Datei erzeugt haben (@PG)",
                "Programs that made the file (@PG)"
            ),
            l("Programm", "Program"),
            l("Aufruf", "Command line")
        );
        for p in &programs {
            let _ = write!(
                s,
                "<tr><td>{}</td><td>{}</td><td><code>{}</code></td></tr>",
                esc(p["PN"].as_str().or(p["ID"].as_str()).unwrap_or("")),
                esc(p["VN"].as_str().unwrap_or("")),
                esc(p["CL"].as_str().unwrap_or(""))
            );
        }
        s.push_str("</table>");
    }

    // Satzarten.
    let _ = write!(
        s,
        "<h2>{}</h2>",
        l("Datensätze nach Art", "Records by kind")
    );
    let _ = write!(
        s,
        "<table><tr><th></th><th class=\"r\">{}</th><th class=\"r\">{}</th></tr>",
        l("Anzahl", "Count"),
        l("Anteil", "Share")
    );
    let total = passed["gesamt"].as_u64().unwrap_or(0);
    let of_total = |k: &str| share(passed[k].as_u64().unwrap_or(0), total);
    let of = |k: &str, w: &str| {
        share(
            passed[k].as_u64().unwrap_or(0),
            passed[w].as_u64().unwrap_or(0),
        )
    };
    for (k, de, en, part) in [
        (
            "gesamt",
            "alle Datensätze (Prüfung bestanden)",
            "all records (QC passed)",
            of_total("gesamt"),
        ),
        (
            "hauptzeilen",
            "Hauptzeilen (ein Satz je Read)",
            "primary (one record per read)",
            of_total("hauptzeilen"),
        ),
        (
            "ergaenzend",
            "ergänzende Alignments (geteilte Reads)",
            "supplementary (split reads)",
            of_total("ergaenzend"),
        ),
        (
            "sekundaer",
            "Nebenzeilen (sekundär)",
            "secondary",
            of_total("sekundaer"),
        ),
        (
            "hauptzeilen_mit_ort",
            "Hauptzeilen mit Ort",
            "primary mapped",
            of("hauptzeilen_mit_ort", "hauptzeilen"),
        ),
        (
            "richtiges_paar",
            "richtige Paare",
            "properly paired",
            of("richtiges_paar", "aus_paaren"),
        ),
        (
            "beide_mit_ort",
            "beide Partner mit Ort",
            "both mates mapped",
            of("beide_mit_ort", "aus_paaren"),
        ),
        (
            "nur_selbst_mit_ort",
            "nur ein Partner mit Ort",
            "singletons",
            of("nur_selbst_mit_ort", "aus_paaren"),
        ),
        (
            "partner_auf_anderem_contig",
            "Partner auf anderem Contig",
            "mate on other contig",
            of("partner_auf_anderem_contig", "aus_paaren"),
        ),
        (
            "partner_auf_anderem_contig_mapq5",
            "… davon MAPQ ≥ 5",
            "… of which MAPQ ≥ 5",
            of("partner_auf_anderem_contig_mapq5", "aus_paaren"),
        ),
        (
            "duplikate_hauptzeilen",
            "Duplikate unter den Hauptzeilen",
            "duplicates among primary",
            of("duplikate_hauptzeilen", "hauptzeilen"),
        ),
    ] {
        s.push_str(&row(&[l(de, en), int(&passed[k]), pct(&part, 3)]));
    }
    let failed = f["saetze"]["nicht_bestanden"]["gesamt"]
        .as_u64()
        .unwrap_or(0);
    s.push_str(&row(&[
        l(
            "Prüfung des Geräts nicht bestanden (Flag 512)",
            "Failed vendor QC (flag 512)",
        ),
        num(failed as f64, 0),
        String::new(),
    ]));
    s.push_str("</table>");
    let _ = write!(
        s,
        "<p class=\"note\">{}</p>",
        l(
            "Zählung wie <code>samtools flagstat</code>. Anteile der Paar-Zeilen beziehen sich auf die Hauptzeilen aus Paaren.",
            "Counted like <code>samtools flagstat</code>. Pair shares refer to primary reads from pairs."
        )
    );

    // Tags.
    let _ = write!(
        s,
        "<h3>{}</h3><table><tr><th>Tag</th><th>{}</th><th class=\"r\">{}</th><th class=\"r\">{}</th></tr>",
        l(
            "Zusatzfelder (Tags) in den Datensätzen",
            "Optional fields (tags) in the records"
        ),
        l("Bedeutung", "Meaning"),
        l("Datensätze", "Records"),
        l("Anteil", "Share")
    );
    let records = f["saetze"]["gesamt"].as_u64().unwrap_or(0);
    for tag in f["tags"].as_array().into_iter().flatten() {
        let name = tag["tag"].as_str().unwrap_or("");
        let (de, en) = tag_meaning(name);
        let _ = write!(
            s,
            "<tr><td><code>{}</code></td><td>{}</td><td class=\"r\">{}</td><td class=\"r\">{}</td></tr>",
            esc(name),
            l(de, en),
            int(&tag["saetze"]),
            pct(&share(tag["saetze"].as_u64().unwrap_or(0), records), 1)
        );
    }
    s.push_str("</table>");

    // Reads.
    let _ = write!(s, "<h2>{}</h2>", l("Reads und Basen", "Reads and bases"));
    let rl = &f["read_laenge"]["kennzahlen"];
    let q = &f["qualitaet"];
    let e = &f["fehler"];
    s.push_str("<table>");
    s.push_str(&row(&[
        l(
            "Read-Länge (Median, 5–95 %)",
            "Read length (median, 5–95 %)",
        ),
        format!(
            "{} ({}–{})",
            int(&rl["median"]),
            int(&rl["p05"]),
            int(&rl["p95"])
        ),
    ]));
    s.push_str(&row(&[
        l("mittlere Basenqualität", "Mean base quality"),
        match q["mittel"].as_f64() {
            Some(v) => format!("Q{}", num(v, 1)),
            None => "–".into(),
        },
    ]));
    s.push_str(&row(&[
        l("Basen mit Qualität ≥ Q30", "Bases with quality ≥ Q30"),
        pct(&q["anteil_q30"], 2),
    ]));
    s.push_str(&row(&[
        l("ausgerichtete Basen", "Aligned bases"),
        int(&e["ausgerichtete_basen"]),
    ]));
    s.push_str(&row(&[
        l("Fehlpaarungen", "Mismatches"),
        format!("{} ({})", int(&e["fehlpaarungen"]), pct(&e["rate"], 3)),
    ]));
    s.push_str(&row(&[
        l(
            "Fehlerrate mit Indels (wie <code>samtools stats</code>)",
            "Error rate with indels (like <code>samtools stats</code>)",
        ),
        pct(&e["rate_wie_samtools"], 3),
    ]));
    s.push_str(&row(&[
        l("eingefügte Basen (I)", "Inserted bases (I)"),
        int(&e["insertierte_basen"]),
    ]));
    s.push_str(&row(&[
        l("fehlende Basen (D)", "Deleted bases (D)"),
        int(&e["deletierte_basen"]),
    ]));
    s.push_str(&row(&[
        l("weich abgeschnittene Basen (S)", "Soft-clipped bases (S)"),
        int(&e["soft_clip_basen"]),
    ]));
    s.push_str("</table>");
    let _ = write!(
        s,
        "<p class=\"note\">{}</p>",
        l(
            "Länge und Qualität über alle Hauptzeilen. Fehlpaarungen = NM minus eingefügte und fehlende Basen, über alle Sätze mit Ort ohne Duplikate und Nebenzeilen.",
            "Length and quality over all primary reads. Mismatches = NM minus inserted and deleted bases, over all mapped records without duplicates and secondary lines."
        )
    );
    if viz["read_laenge"].is_array() {
        let _ = write!(
            s,
            "<h3>{}</h3><div class=\"viz\" id=\"v-readlen\"></div>",
            l(
                "Verteilung der Read-Längen (nach dem Trimmen)",
                "Read length distribution (after trimming)"
            )
        );
    }

    // MAPQ.
    let _ = write!(
        s,
        "<h2>{}</h2>",
        l(
            "Wie sicher der Ort ist (MAPQ)",
            "How certain the location is (MAPQ)"
        )
    );
    s.push_str("<div class=\"viz\" id=\"v-mapq\"></div>");
    let _ = write!(
        s,
        "<p class=\"note\">{} {} · MAPQ 0: {} · MAPQ ≥ 20: {} · MAPQ 60: {}</p>",
        l(
            "Hauptzeilen mit Ort, ohne Duplikate:",
            "Primary mapped reads without duplicates:"
        ),
        int(&f["mapq"]["reads"]),
        pct(&f["mapq"]["anteil_0"], 2),
        pct(&f["mapq"]["anteil_ab_20"], 2),
        pct(&f["mapq"]["anteil_60"], 2)
    );
    let _ = write!(
        s,
        "<p class=\"note\">{}</p>",
        l(
            "MAPQ 0 heißt: Der Read passt gleich gut an mehrere Orte (Wiederholungen, Segmentduplikationen). 60 ist der höchste Wert, den bwa und Genome vergeben.",
            "MAPQ 0 means the read fits equally well at several places (repeats, segmental duplications). 60 is the highest value bwa and Genome assign."
        )
    );

    // Insertgröße.
    let ins = &f["insertgroesse"]["kennzahlen"];
    if ins["n"].as_u64().unwrap_or(0) > 0 {
        let _ = write!(s, "<h2>{}</h2>", l("Insertgröße", "Insert size"));
        s.push_str("<div class=\"viz\" id=\"v-insert\"></div>");
        let _ = write!(
            s,
            "<p class=\"note\">{}: {} · {}: {} · {}: {} ± {} · 5–95 %: {}–{} · {}: {}</p>",
            l("Paare", "Pairs"),
            int(&ins["n"]),
            l("Median", "Median"),
            int(&ins["median"]),
            l(
                "Mittel ± Streuung (Paare bis 2.000 Basen; samtools schneidet Ausreißer anders ab)",
                "Mean ± SD (pairs up to 2,000 bases; samtools trims outliers differently)"
            ),
            num(ins["mittel"].as_f64().unwrap_or(0.0), 1),
            num(ins["streuung"].as_f64().unwrap_or(0.0), 1),
            int(&ins["p05"]),
            int(&ins["p95"]),
            l(
                &format!("über {}", gc_bam::qc::MAX_INSERT),
                &format!("above {}", gc_bam::qc::MAX_INSERT)
            ),
            int(&f["insertgroesse"]["ueber_max"])
        );
    }

    // Abdeckung.
    let _ = write!(s, "<h2>{}</h2>", l("Abdeckung", "Coverage"));
    let mean = cov["mittel"].as_f64().unwrap_or(0.0);
    let _ = write!(
        s,
        "<h3>{}</h3><div class=\"viz\" id=\"v-depth\"></div><h3>{}</h3><div class=\"viz\" id=\"v-depth-cum\"></div>",
        l(
            "Verteilung der Abdeckung je Stelle (Autosomen)",
            "Per-position coverage distribution (autosomes)"
        ),
        l(
            "Anteil der Stellen mit mindestens so viel Abdeckung (Autosomen)",
            "Share of positions with at least this coverage (autosomes)"
        )
    );
    let _ = write!(
        s,
        "<p class=\"note\">{} {}x · {} {}x · ≥ 1x: {} · ≥ 10x: {} · ≥ 20x: {} · ≥ 30x: {}</p>",
        l("Mittel", "Mean"),
        num(mean, 2),
        l("Median", "median"),
        int(&cov["median"]),
        pct(&cov["anteil_ab"]["1"], 2),
        pct(&cov["anteil_ab"]["10"], 2),
        pct(&cov["anteil_ab"]["20"], 2),
        pct(&cov["anteil_ab"]["30"], 2)
    );

    let _ = write!(
        s,
        "<h3>{}</h3><table><tr><th>Contig</th><th class=\"r\">{}</th><th class=\"r\">Reads</th><th class=\"r\">{}</th><th class=\"r\">{}</th><th class=\"r\">≥ 1x</th><th class=\"r\">≥ 10x</th><th class=\"r\">≥ 20x</th><th class=\"r\">≥ 30x</th><th class=\"r\">{}</th></tr>",
        l("Je Chromosom", "Per chromosome"),
        l("Länge", "Length"),
        l("Duplikate", "Duplicates"),
        l("Abdeckung", "Coverage"),
        l("MAPQ (Mittel)", "MAPQ (mean)")
    );
    let contig_row = |name: String,
                      len: u64,
                      reads: u64,
                      dups: u64,
                      bases: u64,
                      covered: [u64; 4],
                      mapq: Option<f64>|
     -> String {
        let mut cells = vec![
            name,
            num(len as f64, 0),
            num(reads as f64, 0),
            pct(&share(dups, reads), 2),
            if len > 0 {
                format!("{}x", num(bases as f64 / len as f64, 2))
            } else {
                "–".into()
            },
        ];
        for c in covered {
            cells.push(pct(&share(c, len), 1));
        }
        cells.push(mapq.map(|m| num(m, 1)).unwrap_or_else(|| "–".into()));
        row(&cells)
    };
    let mut others: Vec<(String, GroupSum)> = Vec::new();
    for c in f["abdeckung"]["contigs"].as_array().into_iter().flatten() {
        let g = c["gruppe"].as_str().unwrap_or("");
        let covered: [u64; 4] = std::array::from_fn(|i| c["stellen_ab"][i].as_u64().unwrap_or(0));
        let (len, reads, dups, bases) = (
            c["laenge"].as_u64().unwrap_or(0),
            c["reads"].as_u64().unwrap_or(0),
            c["duplikate"].as_u64().unwrap_or(0),
            c["basen"].as_u64().unwrap_or(0),
        );
        if matches!(g, "autosom" | "geschlecht" | "mito") {
            let name = esc(c["name"].as_str().unwrap_or(""));
            let line = contig_row(
                name.clone(),
                len,
                reads,
                dups,
                bases,
                covered,
                c["mapq_mittel"].as_f64(),
            );
            s.push_str(&line.replacen(
                "<tr>",
                &format!("<tr class=\"pick\" data-c=\"{name}\">"),
                1,
            ));
        } else {
            let entry = match others.iter_mut().find(|o| o.0 == g) {
                Some(e) => e,
                None => {
                    others.push((g.to_string(), (0, 0, 0, 0, [0; 4], 0)));
                    others.last_mut().expect("eben angelegt")
                }
            };
            let v = &mut entry.1;
            v.0 += len;
            v.1 += reads;
            v.2 += dups;
            v.3 += bases;
            for i in 0..4 {
                v.4[i] += covered[i];
            }
            v.5 += 1;
        }
    }
    for (g, (len, reads, dups, bases, covered, n)) in others {
        let (de, en) = group_name(&g);
        let name = format!(
            "{} <span class=\"sub\">({})</span>",
            l(de, en),
            l(&format!("{n} Contigs"), &format!("{n} contigs"))
        );
        s.push_str(&contig_row(name, len, reads, dups, bases, covered, None));
    }
    let unmapped = f["saetze"]["hauptzeilen_ohne_ort"].as_u64().unwrap_or(0);
    s.push_str(&row(&[
        l("ohne Ort", "Unmapped"),
        String::new(),
        num(unmapped as f64, 0),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
    ]));
    s.push_str("</table>");
    let gx = &f["abdeckung"]["geschlecht"];
    let _ = write!(
        s,
        "<p class=\"note\">{} X/A = {} · Y/A = {} → {}. {}</p>",
        l(
            "Abdeckung der Geschlechtschromosomen im Verhältnis zu den Autosomen:",
            "Sex chromosome coverage relative to the autosomes:"
        ),
        num(gx["x_zu_autosomen"].as_f64().unwrap_or(f64::NAN), 2),
        num(gx["y_zu_autosomen"].as_f64().unwrap_or(f64::NAN), 2),
        match sex {
            "XY" => l("passt zu XY", "consistent with XY"),
            "XX" => l("passt zu XX", "consistent with XX"),
            _ => l("nicht eindeutig", "not clear"),
        },
        l(
            "Abdeckung = ausgerichtete Basen ohne Duplikate, Nebenzeilen und nicht bestandene Sätze, geteilt durch die Länge (N-Strecken eingerechnet, wie <code>samtools coverage</code>). Reads je Contig: Hauptzeilen mit Ort, Duplikate eingeschlossen.",
            "Coverage = aligned bases without duplicates, secondary and QC-failed records, divided by the length (N gaps included, like <code>samtools coverage</code>). Reads per contig: primary mapped reads, duplicates included."
        )
    );

    // Lauf.
    if !lauf.is_null() {
        let r = &lauf["ergebnis"];
        let _ = write!(
            s,
            "<h2>{}</h2>",
            l("Wie die BAM entstanden ist", "How the BAM was made")
        );
        s.push_str("<table>");
        let p = &lauf["parameter"];
        s.push_str(&row(&[
            l("Werkzeug", "Tool"),
            format!(
                "Genome {} · <code>fq2bam</code>",
                esc(lauf["version"].as_str().unwrap_or(""))
            ),
        ]));
        for (k, de, en) in [
            ("fastq", "FASTQ Read 1", "FASTQ read 1"),
            ("fastq2", "FASTQ Read 2", "FASTQ read 2"),
            ("index", "Referenz-Index", "Reference index"),
        ] {
            if let Some(v) = p[k].as_str() {
                s.push_str(&row(&[l(de, en), format!("<code>{}</code>", esc(v))]));
            }
        }
        if let (Some(b), Some(e)) = (lauf["beginn"].as_str(), lauf["ende"].as_str()) {
            s.push_str(&row(&[
                l("Beginn – Ende (UTC)", "Start – end (UTC)"),
                format!("{} – {}", esc(b), esc(e)),
            ]));
        }
        if let Some(sys) = lauf["system"].as_object() {
            s.push_str(&row(&[
                l("Rechner", "Machine"),
                format!(
                    "{} {} · {} {}",
                    esc(sys.get("system").and_then(|v| v.as_str()).unwrap_or("")),
                    esc(sys
                        .get("architektur")
                        .and_then(|v| v.as_str())
                        .unwrap_or("")),
                    int(&p["kerne"]),
                    l("Kerne benutzt", "cores used")
                ),
            ]));
        }
        s.push_str(&row(&[
            l("Trimmen", "Trimming"),
            if p["trimmen"] == true {
                l("an", "on")
            } else {
                l("aus", "off")
            },
        ]));
        s.push_str("</table>");

        if let Some(steps) = r["sekunden_schritte"].as_object() {
            let _ = write!(s, "<h3>{}</h3><table>", l("Laufzeit", "Runtime"));
            for (k, de, en) in [
                ("index_laden", "Index laden", "Load index"),
                (
                    "ausrichten_und_teile",
                    "Trimmen, Ausrichten, Sortierteile",
                    "Trim, align, sort parts",
                ),
                (
                    "zusammenfuehren_markieren_schreiben",
                    "Zusammenführen, Duplikate markieren, BAM schreiben",
                    "Merge, mark duplicates, write BAM",
                ),
                if p["bericht"] == true {
                    (
                        "bam_index",
                        "Index und Zahlen für den Überblick",
                        "Index and overview numbers",
                    )
                } else {
                    ("bam_index", "Index", "Index")
                },
            ] {
                if let Some(v) = steps.get(k).and_then(|v| v.as_f64()) {
                    s.push_str(&row(&[l(de, en), duration(v)]));
                }
            }
            if let Some(v) = r["sekunden"].as_f64() {
                s.push_str(&row(&[
                    format!("<b>{}</b>", l("gesamt", "total")),
                    format!("<b>{}</b>", duration(v)),
                ]));
            }
            if let Some(v) = r["reads_je_sekunde"].as_f64() {
                s.push_str(&row(&[
                    l("Reads je Sekunde", "Reads per second"),
                    num(v, 0),
                ]));
            }
            s.push_str("</table>");
        }

        let tr = &r["trimmen"];
        if tr.is_object() {
            let _ = write!(s, "<h3>{}</h3><table>", l("Trimmen", "Trimming"));
            let pairs = tr["paare"].as_u64().unwrap_or(0);
            for (k, de, en) in [
                ("paare", "Paare gelesen", "Pairs read"),
                ("behalten", "behalten", "kept"),
                (
                    "verworfen_qualitaet",
                    "verworfen: Qualität",
                    "dropped: quality",
                ),
                (
                    "verworfen_zu_kurz",
                    "verworfen: zu kurz",
                    "dropped: too short",
                ),
                (
                    "verworfen_zu_viele_n",
                    "verworfen: zu viele N",
                    "dropped: too many N",
                ),
                (
                    "verworfen_adapter_dimer",
                    "verworfen: Adapter-Dimer",
                    "dropped: adapter dimer",
                ),
                (
                    "ueberlappung_gekuerzt_paare",
                    "über die Überlappung gekürzt (Paare)",
                    "trimmed by overlap (pairs)",
                ),
                (
                    "adapterfolge_gekuerzt_read1",
                    "Adapterfolge gekürzt, Read 1",
                    "adapter sequence trimmed, read 1",
                ),
                (
                    "adapterfolge_gekuerzt_read2",
                    "Adapterfolge gekürzt, Read 2",
                    "adapter sequence trimmed, read 2",
                ),
                (
                    "adapter_basen",
                    "Adapterbasen entfernt",
                    "adapter bases removed",
                ),
                (
                    "polyg_reads",
                    "PolyG gekürzt (Reads)",
                    "polyG trimmed (reads)",
                ),
            ] {
                let v = &tr[k];
                let extra = if k != "adapter_basen" && k != "paare" && pairs > 0 {
                    format!(
                        " <span class=\"sub\">({})</span>",
                        pct(&share(v.as_u64().unwrap_or(0), pairs), 3)
                    )
                } else {
                    String::new()
                };
                s.push_str(&row(&[l(de, en), format!("{}{extra}", int(v))]));
            }
            for a in lauf["adapter"].as_array().into_iter().flatten() {
                s.push_str(&row(&[
                    format!(
                        "{} {}",
                        l("Adapter", "Adapter"),
                        esc(a["seite"].as_str().unwrap_or(""))
                    ),
                    format!(
                        "<code>{}</code> {}",
                        esc(a["wahl"].as_str().unwrap_or("")),
                        esc(a["name"].as_str().unwrap_or(""))
                    ),
                ]));
            }
            if let Some(pg) = lauf["polyg"].as_object() {
                s.push_str(&row(&[
                    l("PolyG kürzen", "polyG trimming"),
                    format!(
                        "{} <span class=\"sub\">({})</span>",
                        esc(pg.get("wahl").and_then(|v| v.as_str()).unwrap_or("")),
                        esc(pg.get("grund").and_then(|v| v.as_str()).unwrap_or(""))
                    ),
                ]));
            }
            s.push_str("</table>");
        }

        let d = &r["duplikate"];
        if d.is_object() {
            let _ = write!(
                s,
                "<h3>{}</h3><table>",
                l("Duplikate markieren", "Duplicate marking")
            );
            s.push_str(&row(&[
                l("bewertete Reads", "Reads examined"),
                int(&d["bewertet"]),
            ]));
            s.push_str(&row(&[
                l("Duplikate aus Paaren", "Duplicates from pairs"),
                int(&d["duplikate_mit_partner"]),
            ]));
            s.push_str(&row(&[
                l("Duplikate ohne Partner", "Duplicates without mate"),
                int(&d["duplikate_ohne_partner"]),
            ]));
            s.push_str(&row(&[
                l("Anteil", "Share"),
                pct(&d["anteil_duplikate"], 3),
            ]));
            s.push_str(&row(&[
                l(
                    "geschätzte Bibliotheksgröße (Moleküle)",
                    "Estimated library size (molecules)",
                ),
                int(&d["geschaetzte_bibliotheksgroesse"]),
            ]));
            s.push_str("</table>");
        }
        if let Some(cmd) = lauf["befehl"].as_array() {
            let text = cmd
                .iter()
                .filter_map(|v| v.as_str())
                .collect::<Vec<_>>()
                .join(" ");
            let _ = write!(
                s,
                "<details><summary>{}</summary><p><code>{}</code></p></details>",
                l("Vollständiger Aufruf", "Full command line"),
                esc(&text)
            );
        }
    }

    let _ = write!(
        s,
        "<footer>{} <code>{}</code>.</footer>",
        l(
            &format!(
                "Erstellt von Genome {} in einem Durchlauf über alle Datensätze. Alle Zahlen stehen auch maschinenlesbar in",
                esc(f["version"].as_str().unwrap_or(""))
            ),
            &format!(
                "Made by Genome {} in one pass over all records. All numbers are also available machine-readable in",
                esc(f["version"].as_str().unwrap_or(""))
            )
        ),
        SCHEMA
    );
    // Daten der Grafiken; `</` darf im Skript-Block nicht vorkommen.
    let data = serde_json::to_string(&viz)
        .unwrap_or_else(|_| "{}".into())
        .replace("</", "<\\/");
    let _ = writeln!(
        s,
        "</main><script type=\"application/json\" id=\"g-daten\">{data}</script><script>{JS}</script></body></html>"
    );
    s
}

/// Pfade der Ausgaben zu einer BAM: `<name ohne .bam>.html` und `<name ohne .bam>.bamqc.json`.
pub fn output_paths(bam: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let base = if bam.extension().is_some_and(|e| e == "bam") {
        bam.with_extension("")
    } else {
        bam.to_path_buf()
    };
    let mut html = base.as_os_str().to_owned();
    html.push(".html");
    let mut json = base.as_os_str().to_owned();
    json.push(".bamqc.json");
    (html.into(), json.into())
}

/// Schreibt beide Dateien unter vorläufigem Namen und benennt sie dann um.
pub fn write(facts: &Value, html_path: &Path, json_path: &Path) -> anyhow::Result<()> {
    use anyhow::Context;
    for (path, text) in [
        (json_path, serde_json::to_string_pretty(facts)? + "\n"),
        (html_path, html(facts)),
    ] {
        let mut part = path.as_os_str().to_owned();
        part.push(".unfertig");
        let part = std::path::PathBuf::from(part);
        std::fs::write(&part, text).with_context(|| format!("{} schreiben", part.display()))?;
        std::fs::rename(&part, path).with_context(|| format!("{} umbenennen", part.display()))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zahlen_in_beiden_sprachen() {
        assert_eq!(
            num_pair(1234567.891, 2),
            ("1.234.567,89".into(), "1,234,567.89".into())
        );
        assert_eq!(num_pair(-1234.0, 0), ("-1.234".into(), "-1,234".into()));
        assert_eq!(num_pair(12.0, 0), ("12".into(), "12".into()));
    }

    #[test]
    fn contig_gruppen() {
        assert_eq!(contig_group("chr1"), "autosom");
        assert_eq!(contig_group("chrX"), "geschlecht");
        assert_eq!(contig_group("chrM"), "mito");
        assert_eq!(contig_group("chrUn_JTFH01000001v1_decoy"), "decoy");
        assert_eq!(contig_group("chr6_GL000250v2_alt"), "alt");
        assert_eq!(contig_group("HLA-A*01:01:01:01"), "hla");
        assert_eq!(contig_group("chrEBV"), "ebv");
        assert_eq!(contig_group("chr1_KI270706v1_random"), "unplatziert");
        assert_eq!(contig_group("chrUn_GL000220v1"), "unplatziert");
    }

    #[test]
    fn ausgabepfade() {
        let (h, j) = output_paths(Path::new("/x/PJ.g4.bam"));
        assert_eq!(h, Path::new("/x/PJ.g4.html"));
        assert_eq!(j, Path::new("/x/PJ.g4.bamqc.json"));
    }
}
