//! Misst die Seed-Suche allein: Read für Read gegen verschränkt. Nur zum Messen.
//! Aufruf: seedbench INDEX FASTQ [GRUPPE]
use gc_align::seed::{SeedOptions, collect_intervals, collect_intervals_many};
use std::io::BufRead;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut index = gc_index::FmIndex::read(std::io::BufReader::new(
        std::fs::File::open(&args[1]).unwrap(),
    ))
    .unwrap();
    // OCC_PLAIN=1: Zählen mit Schleifen wie früher.
    index.set_occ_plain(std::env::var("OCC_PLAIN").is_ok_and(|v| v == "1"));
    let group: usize = args.get(3).map_or(32, |g| g.parse().unwrap());
    let reads: Vec<Vec<u8>> = std::io::BufReader::new(std::fs::File::open(&args[2]).unwrap())
        .lines()
        .skip(1)
        .step_by(4)
        .take(std::env::var("N").map_or(200_000, |n| n.parse().unwrap()))
        .map(|l| {
            l.unwrap()
                .bytes()
                .map(|b| gc_index::reference::encode(b).unwrap_or(4))
                .collect()
        })
        .collect();
    let opt = SeedOptions::default();
    let t = std::time::Instant::now();
    let mut n = 0;
    for r in &reads {
        n += collect_intervals(&index, r, &opt).len();
    }
    let single = t.elapsed();
    let t = std::time::Instant::now();
    let mut m = 0;
    for chunk in reads.chunks(group) {
        let v: Vec<&[u8]> = chunk.iter().map(|r| &r[..]).collect();
        m += collect_intervals_many(&index, &v, &opt)
            .iter()
            .map(Vec::len)
            .sum::<usize>();
    }
    let many = t.elapsed();
    println!(
        "reads {} treffer {n} {m} einzeln {single:?} verschraenkt({group}) {many:?}",
        reads.len()
    );
}
