//! `genome sort` als Ganzes: richtige Reihenfolge, kein Satz verloren, Kopf angepasst, keine
//! Reste im Temp-Ordner, und nach einem Fehler keine Ausgabedatei.

use std::path::Path;
use std::process::Command;

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

/// Ein BAM-Satz ohne Sequenz, mit Längenfeld. Der Name ist die laufende Nummer der Eingabe.
fn record(nr: u32, ref_id: i32, pos: i32, flag: u16) -> Vec<u8> {
    let name = format!("r{nr:07}");
    let mut body = Vec::new();
    body.extend_from_slice(&ref_id.to_le_bytes());
    body.extend_from_slice(&pos.to_le_bytes());
    body.push(name.len() as u8 + 1);
    body.push(0); // MAPQ
    body.extend_from_slice(&4680u16.to_le_bytes()); // bin
    body.extend_from_slice(&0u16.to_le_bytes()); // CIGAR-Schritte
    body.extend_from_slice(&flag.to_le_bytes());
    body.extend_from_slice(&0u32.to_le_bytes()); // Länge der Sequenz
    body.extend_from_slice(&(-1i32).to_le_bytes());
    body.extend_from_slice(&(-1i32).to_le_bytes());
    body.extend_from_slice(&0i32.to_le_bytes());
    body.extend_from_slice(name.as_bytes());
    body.push(0);
    let mut out = (body.len() as u32).to_le_bytes().to_vec();
    out.extend_from_slice(&body);
    out
}

const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted\tGO:query\n@SQ\tSN:a\tLN:100000\n@SQ\tSN:b\tLN:100000\n@SQ\tSN:c\tLN:100000\n@PG\tID:genome\tPN:genome\tVN:0\tCL:x\n";

/// 60.000 Sätze auf drei Contigs mit vielen gleichen Orten; jeder zwanzigste ohne Contig.
fn make_input(path: &Path) -> Vec<Vec<u8>> {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let mut writer = gc_io::bam::create_from_text(path, HEADER, 1).unwrap();
    let mut all = Vec::new();
    for nr in 0..60_000u32 {
        let flag = if rng.below(2) == 0 { 16 } else { 0 };
        let rec = if rng.below(20) == 0 {
            record(nr, -1, -1, 4)
        } else {
            record(nr, rng.below(3) as i32, rng.below(500) as i32, flag)
        };
        gc_io::bam::write_encoded(&mut writer, &rec).unwrap();
        all.push(rec);
    }
    gc_io::bam::finish(writer).unwrap();
    all
}

fn read_all(path: &Path) -> (String, Vec<Vec<u8>>) {
    let (mut reader, _, raw) = gc_io::bam::open_raw(path, 1).unwrap();
    let mut out = Vec::new();
    loop {
        let mut rec = Vec::new();
        if !gc_io::bam::read_encoded(&mut reader, &mut rec).unwrap() {
            break;
        }
        out.push(rec);
    }
    (String::from_utf8(raw).unwrap(), out)
}

fn name(rec: &[u8]) -> &[u8] {
    &rec[36..36 + rec[12] as usize]
}

#[test]
fn sort_orders_like_samtools_and_loses_nothing() {
    let dir = std::env::temp_dir().join(format!("genome-sortieren-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let src = dir.join("ein.bam");
    let input = make_input(&src);

    // Erwartung: stabil nach dem Schlüssel sortiert.
    let mut expected = input.clone();
    expected.sort_by_key(|r| gc_bam::sort_key(&r[4..]));

    let mut outputs = Vec::new();
    // 1 MB erzwingt Teile auf der Platte (die Sätze brauchen rund 2,8 MB), 64 MB nicht.
    for mb in ["1", "64"] {
        let dst = dir.join(format!("aus-{mb}.bam"));
        let temp = dir.join(format!("temp-{mb}"));
        std::fs::create_dir_all(&temp).unwrap();
        let out = Command::new(env!("CARGO_BIN_EXE_genome"))
            .args(["sort", src.to_str().unwrap(), dst.to_str().unwrap()])
            .args(["--speicher-mb", mb, "--kerne", "2"])
            .args(["--temp", temp.to_str().unwrap()])
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8(out.stdout).unwrap();
        assert!(stdout.contains("\"genome.sort.v1\""), "{stdout}");
        assert!(stdout.contains("\"saetze\": 60000"), "{stdout}");
        assert_eq!(
            stdout.contains("\"teile\": 0"),
            mb == "64",
            "Teile nur bei kleinem Speicher: {stdout}"
        );
        let (header, records) = read_all(&dst);
        assert_eq!(records.len(), expected.len());
        if let Some(i) = (0..records.len()).find(|&i| records[i] != expected[i]) {
            panic!(
                "Satz {i}: {:?} statt {:?}",
                String::from_utf8_lossy(name(&records[i])),
                String::from_utf8_lossy(name(&expected[i]))
            );
        }
        let lines: Vec<&str> = header.lines().collect();
        assert_eq!(lines[0], "@HD\tVN:1.6\tSO:coordinate\tGO:query");
        assert_eq!(lines[1..5], HEADER.lines().collect::<Vec<_>>()[1..5]);
        assert!(
            lines[5].starts_with("@PG\tID:genome-sort\tPN:genome\tPP:genome\t"),
            "{}",
            lines[5]
        );
        assert_eq!(lines.len(), 6);
        assert_eq!(
            std::fs::read_dir(&temp).unwrap().count(),
            0,
            "Temp-Ordner ist nicht leer"
        );
        assert!(!dir.join(format!("aus-{mb}.bam.unfertig")).exists());
        outputs.push(records);
    }
    assert!(outputs[0] == outputs[1]);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_broken_input_leaves_no_output_and_a_clear_error() {
    let dir = std::env::temp_dir().join(format!("genome-sortieren-kaputt-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let whole = dir.join("ganz.bam");
    make_input(&whole);
    // Mitten in den Daten abgeschnitten.
    let bytes = std::fs::read(&whole).unwrap();
    let src = dir.join("ein.bam");
    std::fs::write(&src, &bytes[..bytes.len() / 2]).unwrap();
    let dst = dir.join("aus.bam");
    let log = dir.join("lauf.jsonl");
    let out = Command::new(env!("CARGO_BIN_EXE_genome"))
        .args(["sort", src.to_str().unwrap(), dst.to_str().unwrap()])
        .args(["--speicher-mb", "1", "--log", log.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(!dst.exists(), "trotz Fehler gibt es eine Ausgabedatei");
    assert!(!dir.join("aus.bam.unfertig").exists());
    let log_text = std::fs::read_to_string(&log).unwrap();
    assert!(log_text.contains("SORT.FEHLGESCHLAGEN"), "{log_text}");
    let left: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name())
        .filter(|n| n.to_string_lossy().starts_with("genome-sort-"))
        .collect();
    assert!(left.is_empty(), "Teile blieben liegen: {left:?}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// `genome bam-index`: Der Index nennt für jedes Contig die richtigen Zähler, und jeder
/// vermerkte Dateibereich lässt sich lesen.
#[test]
fn index_counts_match_the_file() {
    let dir = std::env::temp_dir().join(format!("genome-index-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let src = dir.join("ein.bam");
    let input = make_input(&src);
    let sorted = dir.join("sortiert.bam");
    let run = |args: &[&str]| {
        let out = Command::new(env!("CARGO_BIN_EXE_genome"))
            .args(args)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        String::from_utf8(out.stdout).unwrap()
    };
    run(&["sort", src.to_str().unwrap(), sorted.to_str().unwrap()]);
    let out = run(&["bam-index", sorted.to_str().unwrap()]);
    assert!(out.contains("\"genome.bam-index.v1\""), "{out}");
    let bai = gc_bam::Bai::read(&std::fs::read(dir.join("sortiert.bam.bai")).unwrap()).unwrap();
    assert_eq!(bai.references.len(), 3);
    // Zähler im Sonderfach 37450: (mit Ort, ohne Ort) je Contig. Alle Testsätze mit Contig
    // tragen kein Flag 4.
    let mut total = 0;
    for (tid, reference) in bai.references.iter().enumerate() {
        let expected = input
            .iter()
            .filter(|r| i32::from_le_bytes([r[4], r[5], r[6], r[7]]) == tid as i32)
            .count() as u64;
        assert_eq!(reference.bins[&37450][1], (expected, 0), "Contig {tid}");
        // Positionen 0 bis 499 liegen alle im ersten Fenster von 16.384 Basen.
        assert_eq!(reference.linear.len(), 1);
        total += expected;
    }
    assert_eq!(bai.unplaced, input.len() as u64 - total);
    // Eine unsortierte Datei wird abgelehnt, ohne dass ein Index entsteht.
    let failed = Command::new(env!("CARGO_BIN_EXE_genome"))
        .args(["bam-index", src.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("BAMINDEX.FEHLGESCHLAGEN"));
    assert!(!dir.join("ein.bam.bai").exists());
    let _ = std::fs::remove_dir_all(&dir);
}

/// `genome flagstat` zählt die Datei aus `make_input` richtig.
#[test]
fn flagstat_counts_the_file() {
    let dir = std::env::temp_dir().join(format!("genome-flagstat-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let src = dir.join("ein.bam");
    let input = make_input(&src);
    let out = Command::new(env!("CARGO_BIN_EXE_genome"))
        .args(["flagstat", src.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(out.status.success());
    let text = String::from_utf8(out.stdout).unwrap();
    let unplaced = input.iter().filter(|r| r[18] & 4 != 0).count();
    assert!(text.contains("\"genome.flagstat.v1\""), "{text}");
    assert!(text.contains("\"gesamt\": 60000"), "{text}");
    assert!(
        text.contains(&format!("\"mit_ort\": {}", 60000 - unplaced)),
        "{text}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
