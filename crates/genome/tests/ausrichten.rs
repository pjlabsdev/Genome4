//! `genome index sam` als Ganzes: Die Ausgabe hängt nicht von der Zahl der Kerne ab, und die
//! BAM-Ausgabe enthält dieselben Zeilen wie die SAM-Ausgabe.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::sam;
use noodles::sam::alignment::io::Write as _;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// Liest ein Stück mit einzelnen falschen Basen.
fn read_from(rng: &mut Rng, piece: &[u8]) -> String {
    let mut read = piece.to_vec();
    for base in &mut read {
        if rng.below(50) == 0 {
            *base = b"ACGT"[rng.below(4) as usize];
        }
    }
    String::from_utf8(read).unwrap()
}

/// Zwei Contigs, 400 Paare mit 100 Basen je Read, Stücke von 250 bis 350 Basen. Jedes 40. Paar
/// hat einen Partner aus Zufallsbasen (ohne Ort).
fn make_data(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let contigs: Vec<Vec<u8>> = (0..2)
        .map(|_| {
            (0..20_000)
                .map(|_| b"ACGT"[rng.below(4) as usize])
                .collect()
        })
        .collect();
    let mut fasta = String::new();
    for (i, contig) in contigs.iter().enumerate() {
        writeln!(fasta, ">c{}", i + 1).unwrap();
        for line in contig.chunks(60) {
            fasta.push_str(std::str::from_utf8(line).unwrap());
            fasta.push('\n');
        }
    }
    let (mut r1, mut r2) = (String::new(), String::new());
    let qual = "I".repeat(100);
    for i in 0..400 {
        let contig = &contigs[rng.below(2) as usize];
        let len = 250 + rng.below(101) as usize;
        let start = rng.below((contig.len() - len) as u64) as usize;
        let piece = &contig[start..start + len];
        let mut first = read_from(&mut rng, &piece[..100]);
        let mut second = read_from(&mut rng, &revcomp(&piece[len - 100..]));
        if i % 40 == 7 {
            second = (0..100)
                .map(|_| b"ACGT"[rng.below(4) as usize] as char)
                .collect();
        }
        if rng.below(2) == 0 {
            std::mem::swap(&mut first, &mut second);
        }
        writeln!(r1, "@p{i}\n{first}\n+\n{qual}").unwrap();
        writeln!(r2, "@p{i}\n{second}\n+\n{qual}").unwrap();
    }
    let paths = (dir.join("ref.fa"), dir.join("r1.fq"), dir.join("r2.fq"));
    std::fs::write(&paths.0, fasta).unwrap();
    std::fs::write(&paths.1, r1).unwrap();
    std::fs::write(&paths.2, r2).unwrap();
    paths
}

fn genome(args: &[&str]) -> Vec<u8> {
    let out = Command::new(env!("CARGO_BIN_EXE_genome"))
        .args(args)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "genome {args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

/// Die Zeilen ohne den Kopf (dort steht der Aufruf, der sich unterscheidet).
fn records(sam_text: &[u8]) -> Vec<String> {
    String::from_utf8(sam_text.to_vec())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(str::to_string)
        .collect()
}

/// Liest eine BAM und gibt ihre Sätze als SAM-Zeilen zurück.
fn bam_records(path: &Path) -> Vec<String> {
    let (mut reader, header) = gc_io::bam::open(path, 1).unwrap();
    let mut writer = sam::io::Writer::new(Vec::new());
    for record in reader.record_bufs(&header) {
        writer
            .write_alignment_record(&header, &record.unwrap())
            .unwrap();
    }
    records(writer.get_ref())
}

#[test]
fn output_is_the_same_for_any_core_count_and_as_bam() {
    let dir = std::env::temp_dir().join(format!("genome-ausrichten-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (fasta, r1, r2) = make_data(&dir);
    let index = dir.join("ref.gcfm");
    let s = |p: &Path| p.to_str().unwrap().to_string();
    genome(&["index", "build", &s(&fasta), &s(&index)]);

    for paired in [false, true] {
        let run = |extra: &[&str]| {
            let (index, r1, r2) = (s(&index), s(&r1), s(&r2));
            // Kleine Blöcke, damit die Insertgröße mehrmals geschätzt wird.
            let mut args = vec!["index", "sam", &index, "--fastq", &r1];
            args.extend(["--block-basen", "20000"]);
            if paired {
                args.extend(["--fastq2", &r2]);
            }
            args.extend(extra);
            genome(&args)
        };
        let one = records(&run(&["--kerne", "1"]));
        let three = records(&run(&["--kerne", "3"]));
        assert!(one.len() >= if paired { 800 } else { 400 });
        assert_eq!(one, three, "paired={paired}");

        let mapped = one
            .iter()
            .filter(|l| l.split('\t').nth(1).unwrap().parse::<u32>().unwrap() & 4 == 0)
            .count();
        assert!(mapped * 10 >= one.len() * 9, "paired={paired}: {mapped}");

        let bam = dir.join("aus.bam");
        let text = run(&["--kerne", "2", "--bam", &s(&bam)]);
        assert!(text.is_empty(), "mit --bam bleibt die Standardausgabe leer");
        assert_eq!(bam_records(&bam), one, "paired={paired}");
    }
    std::fs::remove_dir_all(dir).ok();
}

/// Paralleles Lesen und Trimmen ergibt dieselben Blöcke und damit dieselbe BAM wie nacheinander,
/// auch mit Stückgrenzen mitten in einem Block.
#[test]
fn parallel_reading_gives_the_same_bam() {
    let dir = std::env::temp_dir().join(format!("genome-lesen-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (fasta, r1, r2) = make_data(&dir);
    let index = dir.join("ref.gcfm");
    let s = |p: &Path| p.to_str().unwrap().to_string();
    genome(&["index", "build", &s(&fasta), &s(&index)]);
    let run = |parallel: &str, chunk: &str, trim: bool| {
        let bam = dir.join(format!("p{parallel}-{chunk}-{trim}.bam"));
        let (i, a, b, o) = (s(&index), s(&r1), s(&r2), s(&bam));
        let mut args = vec!["fq2bam", "--fastq", &a, "--fastq2", &b, "--aus", &o];
        args.extend(["--block-basen", "20000", "--kerne", "3", "--ohne-bericht"]);
        if trim {
            args.push("--trimmen");
        }
        args.push(&i);
        let out = Command::new(env!("CARGO_BIN_EXE_genome"))
            .args(&args)
            .env("GC_LESEN_PARALLEL", parallel)
            .env("GC_LESEN_STUECK", chunk)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        bam_records(&bam)
    };
    for trim in [false, true] {
        let serial = run("0", "8192", trim);
        assert!(serial.len() >= 790, "{}", serial.len());
        for chunk in ["7", "64", "8192"] {
            assert_eq!(
                run("1", chunk, trim),
                serial,
                "Stück {chunk}, trimmen {trim}"
            );
        }
    }
    std::fs::remove_dir_all(dir).ok();
}

/// N in der Referenz zählt als unbekannt: Abzug 1, `N` in `MD`, gezählt in `NM`. Das Ergebnis
/// hängt nicht von der Füllbase ab, die der Index dort für die Suche führt.
#[test]
fn n_in_the_reference_is_unknown_not_a_filler_base() {
    let dir = std::env::temp_dir().join(format!("genome-n-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut rng = Rng(0x1234_5678_9ABC_DEF1);
    let mut contig: Vec<u8> = (0..3000).map(|_| b"ACGT"[rng.below(4) as usize]).collect();
    let read = String::from_utf8(contig[1000..1100].to_vec()).unwrap();
    contig[1050] = b'N';
    contig[1051] = b'N';
    let fasta = dir.join("ref.fa");
    std::fs::write(
        &fasta,
        format!(">c\n{}\n", std::str::from_utf8(&contig).unwrap()),
    )
    .unwrap();
    let fastq = dir.join("r.fq");
    let rc = String::from_utf8(revcomp(read.as_bytes())).unwrap();
    let qual = "I".repeat(100);
    std::fs::write(
        &fastq,
        format!("@vor\n{read}\n+\n{qual}\n@rueck\n{rc}\n+\n{qual}\n"),
    )
    .unwrap();
    let index = dir.join("ref.gcfm");
    let s = |p: &Path| p.to_str().unwrap().to_string();
    genome(&["index", "build", &s(&fasta), &s(&index)]);
    let lines = records(&genome(&[
        "index",
        "sam",
        &s(&index),
        "--fastq",
        &s(&fastq),
    ]));
    assert_eq!(lines.len(), 2);
    for (line, flag) in lines.iter().zip(["0", "16"]) {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!((f[1], f[3], f[5]), (flag, "1001", "100M"), "{line}");
        for tag in ["NM:i:2", "MD:Z:50N0N48", "AS:i:96"] {
            assert!(f[11..].contains(&tag), "{tag} fehlt: {line}");
        }
    }
    std::fs::remove_dir_all(dir).ok();
}

fn tiny_case(dir: &Path) -> (String, String) {
    let mut rng = Rng(0x0BAD_CAFE_1234_0001);
    let contig: Vec<u8> = (0..3000).map(|_| b"ACGT"[rng.below(4) as usize]).collect();
    let fasta = dir.join("ref.fa");
    std::fs::write(
        &fasta,
        format!(">c\n{}\n", std::str::from_utf8(&contig).unwrap()),
    )
    .unwrap();
    let qual = "I".repeat(100);
    let mut fq = String::new();
    // ein Read mit Ort, einer aus Zufallsbasen ohne Ort
    let read = std::str::from_utf8(&contig[500..600]).unwrap().to_string();
    let other: String = (0..100)
        .map(|_| b"ACGT"[rng.below(4) as usize] as char)
        .collect();
    writeln!(fq, "@a\n{read}\n+\n{qual}\n@b\n{other}\n+\n{qual}").unwrap();
    let fastq = dir.join("r.fq");
    std::fs::write(&fastq, fq).unwrap();
    let index = dir.join("ref.gcfm");
    let s = |p: &Path| p.to_str().unwrap().to_string();
    genome(&["index", "build", &s(&fasta), &s(&index)]);
    (s(&index), s(&fastq))
}

/// Lesegruppe: Kopfzeile `@RG` zwischen `@SQ` und `@PG`, `RG` in jeder Zeile (auch ohne Ort),
/// in SAM und BAM gleich. Eine Angabe ohne ID wird mit festem Fehlercode abgelehnt.
#[test]
fn read_group_goes_into_header_and_every_line() {
    let dir = std::env::temp_dir().join(format!("genome-rg-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (index, fastq) = tiny_case(&dir);
    let rg = "@RG\\tID:lauf1\\tSM:probe\\tPL:ILLUMINA";
    let base = [
        "index",
        "sam",
        &index,
        "--fastq",
        &fastq,
        "--lesegruppe",
        rg,
    ];
    let text = String::from_utf8(genome(&base)).unwrap();
    let head: Vec<&str> = text.lines().filter(|l| l.starts_with('@')).collect();
    assert_eq!(head[0], "@SQ\tSN:c\tLN:3000");
    assert_eq!(head[1], "@RG\tID:lauf1\tSM:probe\tPL:ILLUMINA");
    assert!(head[2].starts_with("@PG\t"));
    let lines = records(text.as_bytes());
    assert_eq!(lines.len(), 2);
    for line in &lines {
        assert!(line.split('\t').any(|f| f == "RG:Z:lauf1"), "{line}");
    }
    // bei bwa steht RG hinter XS
    assert!(lines[0].contains("\tXS:i:0\tRG:Z:lauf1"), "{}", lines[0]);

    let bam = dir.join("aus.bam");
    let mut args = base.to_vec();
    let bam_path = bam.to_str().unwrap().to_string();
    args.extend(["--bam", &bam_path]);
    genome(&args);
    assert_eq!(bam_records(&bam), lines);
    let (_, header) = gc_io::bam::open(&bam, 1).unwrap();
    assert_eq!(header.read_groups().len(), 1);

    let bad = Command::new(env!("CARGO_BIN_EXE_genome"))
        .args(["index", "sam", &index, "--fastq", &fastq])
        .args(["--lesegruppe", "@RG\\tSM:probe"])
        .output()
        .unwrap();
    assert!(!bad.status.success());
    assert!(bad.stdout.is_empty(), "bei einem Fehler keine Ausgabe");
    let err = String::from_utf8_lossy(&bad.stderr);
    assert!(err.contains("ALIGN.LESEGRUPPE_UNGUELTIG"), "{err}");
    std::fs::remove_dir_all(dir).ok();
}

/// Die ganze Kette bis zur markierten BAM: ausrichten mit Partner-Angaben, sortieren, Duplikate
/// markieren. Jedes Paar kommt zweimal vor, die Kopie mit schlechterer Qualität.
#[test]
fn chain_to_marked_duplicates() {
    let dir = std::env::temp_dir().join(format!("genome-duplikate-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (fasta, r1, r2) = make_data(&dir);
    let index = dir.join("ref.gcfm");
    let s = |p: &Path| p.to_str().unwrap().to_string();
    genome(&["index", "build", &s(&fasta), &s(&index)]);
    // Kopien anhängen: gleicher Inhalt, Name mit "k" davor, Qualität 5 Stufen tiefer.
    for path in [&r1, &r2] {
        let text = std::fs::read_to_string(path).unwrap();
        let mut copy = String::new();
        for (i, line) in text.lines().enumerate() {
            match i % 4 {
                0 => writeln!(copy, "@k{}", &line[1..]).unwrap(),
                3 => writeln!(copy, "{}", "D".repeat(line.len())).unwrap(),
                _ => writeln!(copy, "{line}").unwrap(),
            }
        }
        std::fs::write(path, text + &copy).unwrap();
    }
    let (raw, sorted, marked) = (
        dir.join("roh.bam"),
        dir.join("sortiert.bam"),
        dir.join("markiert.bam"),
    );
    let (i, a, b) = (s(&index), s(&r1), s(&r2));
    genome(&[
        "index",
        "sam",
        &i,
        "--fastq",
        &a,
        "--fastq2",
        &b,
        "--bam",
        &s(&raw),
        "--partner-angaben",
    ]);
    genome(&["sort", &s(&raw), &s(&sorted)]);
    let out = String::from_utf8(genome(&["markdup", &s(&sorted), &s(&marked)])).unwrap();
    assert!(out.contains("\"genome.markdup.v1\""), "{out}");

    let before = bam_records(&sorted);
    let after = bam_records(&marked);
    assert_eq!(before.len(), after.len());
    let mut dup_originals = 0;
    let mut dup_copies = 0;
    for (x, y) in before.iter().zip(&after) {
        let (fx, fy): (Vec<&str>, Vec<&str>) = (x.split('\t').collect(), y.split('\t').collect());
        // Nur das Flag darf sich ändern, und nur um 1024.
        let (flag_x, flag_y): (u32, u32) = (fx[1].parse().unwrap(), fy[1].parse().unwrap());
        assert_eq!(flag_x, flag_y & !1024, "{x}\n{y}");
        assert_eq!(fx[0], fy[0]);
        assert_eq!(fx[2..], fy[2..]);
        assert!(flag_y & 4 != 0 || y.contains("\tMC:Z:"), "{y}");
        if flag_y & 1024 != 0 {
            if fy[0].starts_with('k') {
                dup_copies += 1;
            } else {
                dup_originals += 1;
            }
        }
    }
    // 400 Paare, 10 davon mit einem Read ohne Ort: Die Kopie verliert immer, das Original nie.
    assert_eq!(dup_originals, 0);
    assert!(dup_copies >= 780, "{dup_copies} Kopien markiert");

    // Der eine Befehl gibt dieselben Zeilen wie die einzelnen Schritte, auch mit wenig Speicher,
    // und schreibt Index und Protokoll daneben.
    for (name, memory) in [("ganz.bam", "2048"), ("ganz-klein.bam", "1")] {
        let whole = dir.join(name);
        let out = String::from_utf8(genome(&[
            "fq2bam",
            &i,
            "--fastq",
            &a,
            "--fastq2",
            &b,
            "--aus",
            &s(&whole),
            "--speicher-mb",
            memory,
            // Vergleich mit den einzelnen Schritten: bwa-gleicher Weg.
            "--strobeweg",
            "aus",
        ]))
        .unwrap();
        assert!(out.contains("\"genome.fq2bam.v1\""), "{out}");
        assert!(
            bam_records(&whole) == after,
            "{name} weicht von den Schritten ab"
        );
        let bai = std::fs::read(dir.join(format!("{name}.bai"))).unwrap();
        assert_eq!(&bai[..4], b"BAI\x01");
        let log = std::fs::read_to_string(dir.join(format!("{name}.log.jsonl"))).unwrap();
        assert!(log.contains("FQ2BAM.FERTIG"), "{log}");
        // Schreibt das Sortieren Teile auf die Platte, steht das im Protokoll.
        assert_eq!(
            log.contains("SORT.TEIL"),
            !out.contains("\"teile\": 0"),
            "{out}"
        );
        assert!(!dir.join(format!("{name}.unfertig")).exists());
        // Der Überblick liegt daneben; `bam-bericht` setzt aus BAM und Protokoll dieselbe Datei.
        let stem = name.trim_end_matches(".bam");
        let html = std::fs::read_to_string(dir.join(format!("{stem}.html"))).unwrap();
        assert!(
            html.contains("<span class=\"en\">BAM overview</span>"),
            "{html}"
        );
        let facts = std::fs::read_to_string(dir.join(format!("{stem}.bamqc.json"))).unwrap();
        assert!(facts.contains("\"genome.bamqc.v1\""));
        assert!(facts.contains("\"FQ2BAM.FERTIG\"") || facts.contains("\"ergebnis\""));
        let again = dir.join(format!("{stem}-neu.html"));
        genome(&["bam-bericht", &s(&whole), "--aus", &s(&again)]);
        let facts_again =
            std::fs::read_to_string(dir.join(format!("{stem}-neu.bamqc.json"))).unwrap();
        assert_eq!(
            facts, facts_again,
            "bam-bericht weicht vom Überblick aus fq2bam ab"
        );
    }
    // Mit --ohne-bericht entsteht kein Überblick, die BAM bleibt gleich.
    let plain_run = dir.join("ohne-bericht.bam");
    genome(&[
        "fq2bam",
        &i,
        "--fastq",
        &a,
        "--fastq2",
        &b,
        "--aus",
        &s(&plain_run),
        "--ohne-bericht",
        "--strobeweg",
        "aus",
    ]);
    assert!(bam_records(&plain_run) == after);
    assert!(!dir.join("ohne-bericht.html").exists());
    assert!(!dir.join("ohne-bericht.bamqc.json").exists());
    // Ungleich lange FASTQ-Dateien: Abbruch mit Fehlercode, keine BAM.
    let short = dir.join("kurz.fq");
    let text = std::fs::read_to_string(&r2).unwrap();
    let cut: Vec<&str> = text.lines().take(40).collect();
    std::fs::write(&short, cut.join("\n") + "\n").unwrap();
    let broken = dir.join("kaputt.bam");
    let failed = Command::new(env!("CARGO_BIN_EXE_genome"))
        .args([
            "fq2bam",
            &i,
            "--fastq",
            &a,
            "--fastq2",
            &s(&short),
            "--aus",
            &s(&broken),
        ])
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(!broken.exists());
    let log = std::fs::read_to_string(dir.join("kaputt.bam.log.jsonl")).unwrap();
    assert!(log.contains("FQ2BAM.FEHLGESCHLAGEN"), "{log}");

    // Ohne Partner-Angaben: klarer Fehler, keine Ausgabe.
    let plain = dir.join("ohne.bam");
    let plain_sorted = dir.join("ohne-sortiert.bam");
    genome(&["index", "sam", &i, "--fastq", &a, "--bam", &s(&plain)]);
    genome(&["sort", &s(&plain), &s(&plain_sorted)]);
    // einzelne Reads brauchen keine Partner-Angaben
    genome(&[
        "markdup",
        &s(&plain_sorted),
        &s(&dir.join("ohne-markiert.bam")),
    ]);
    let stripped = dir.join("paare-ohne-angaben.bam");
    let stripped_sorted = dir.join("paare-ohne-angaben-sortiert.bam");
    genome(&[
        "index",
        "sam",
        &i,
        "--fastq",
        &a,
        "--fastq2",
        &b,
        "--bam",
        &s(&stripped),
    ]);
    genome(&["sort", &s(&stripped), &s(&stripped_sorted)]);
    let log = dir.join("fehler.jsonl");
    let target = dir.join("darf-nicht-entstehen.bam");
    let failed = Command::new(env!("CARGO_BIN_EXE_genome"))
        .args([
            "markdup",
            &s(&stripped_sorted),
            &s(&target),
            "--log",
            &s(&log),
        ])
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(!target.exists());
    let log_text = std::fs::read_to_string(&log).unwrap();
    assert!(log_text.contains("MARKDUP.FEHLGESCHLAGEN"), "{log_text}");
    assert!(log_text.contains("ms"), "{log_text}");
    let _ = std::fs::remove_dir_all(&dir);
}
