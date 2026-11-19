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
