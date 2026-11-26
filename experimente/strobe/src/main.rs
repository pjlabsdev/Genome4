// Messwerkzeug: Randstrobe-Seeds (Syncmere k=20,s=16) auf einer Referenz, Kandidatenort je Read.
use rayon::prelude::*;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const K: usize = 20;
const S: usize = 16;
const T: usize = 3; // Syncmer: kleinstes s-mer an Stelle T-1 (offen)

fn code(b: u8) -> u8 {
    match b { b'A' | b'a' => 0, b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' => 3, _ => 4 }
}
#[inline]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 33; x = x.wrapping_mul(0xff51afd7ed558ccd); x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53); x ^= x >> 33; x
}

/// Syncmere: (Position, Hash des k-mers). Strecken mit N werden übersprungen.
fn syncmers(seq: &[u8], out: &mut Vec<(u32, u64)>) {
    out.clear();
    if seq.len() < K { return; }
    let smask = (1u64 << (2 * S)) - 1;
    let kmask = (1u64 << (2 * K)) - 1;
    let (mut sv, mut kv) = (0u64, 0u64);
    let mut valid = 0usize;
    let mut sh: Vec<u64> = vec![0; seq.len()];
    for (i, &b) in seq.iter().enumerate() {
        let c = code(b);
        if c > 3 { valid = 0; sv = 0; kv = 0; continue; }
        sv = ((sv << 2) | c as u64) & smask;
        kv = ((kv << 2) | c as u64) & kmask;
        valid += 1;
        if valid >= S { sh[i + 1 - S] = mix(sv); }
        if valid >= K {
            let start = i + 1 - K;
            let n = K - S + 1;
            let w = &sh[start..start + n];
            let mut m = 0;
            for j in 1..n { if w[j] < w[m] { m = j; } }
            if m == T - 1 { out.push((start as u32, mix(kv ^ 0x9e3779b97f4a7c15))); }
        }
    }
}

const CANDS: usize = 16;
const WMIN: usize = 2;
const WMAX: usize = 9;
/// Randstrobes: (Position des ersten, Hash).
fn randstrobes(sy: &[(u32, u64)], out: &mut Vec<(u64, u32)>) {
    out.clear();
    for i in 0..sy.len() {
        let lo = i + WMIN;
        if lo >= sy.len() { break; }
        let hi = (i + WMAX).min(sy.len() - 1);
        let h1 = sy[i].1;
        let mut best = lo; let mut bv = u32::MAX;
        for j in lo..=hi {
            let v = (h1 ^ sy[j].1).count_ones();
            if v < bv { bv = v; best = j; }
        }
        out.push(((h1 >> 1).wrapping_add(sy[best].1 / 3), sy[i].0));
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let fa = &args[1];
    let maxocc: usize = std::env::var("MAXOCC").ok().map_or(1000, |v| v.parse().unwrap());
    let t = Instant::now();
    let mut genome = Vec::with_capacity(3_300_000_000usize.min(std::env::var("GENOM_BP").ok().map_or(3_300_000_000, |v| v.parse().unwrap())));
    let name = String::from("alle");
    let mut starts: std::collections::HashMap<String, i64> = Default::default();
    let input: Box<dyn std::io::Read> = if fa == "-" { Box::new(std::io::stdin()) } else { Box::new(std::fs::File::open(fa).unwrap()) };
    for line in BufReader::with_capacity(1 << 20, input).lines() {
        let l = line.unwrap();
        if let Some(n) = l.strip_prefix('>') {
            genome.extend(std::iter::repeat(b'N').take(200));
            starts.insert(n.split_whitespace().next().unwrap().to_string(), genome.len() as i64);
            continue;
        }
        genome.extend(l.bytes().map(|b| b.to_ascii_uppercase()));
    }
    // Contig-Stücke (zwischen den N-Trennern), je Stück Strobes rechnen; zweimal: zählen, füllen.
    let mut cuts: Vec<i64> = starts.values().copied().collect();
    cuts.sort();
    let pieces: Vec<(usize, usize)> = cuts.iter().enumerate().map(|(i, &s0)| {
        let e = if i + 1 < cuts.len() { cuts[i + 1] as usize - 200 } else { genome.len() };
        (s0 as usize, e)
    }).collect();
    const B: u32 = 24;
    let strobes_of = |&(s0, e): &(usize, usize)| -> Vec<(u64, u32)> {
        let mut sy = Vec::new(); let mut st = Vec::new();
        syncmers(&genome[s0..e], &mut sy);
        randstrobes(&sy, &mut st);
        for x in st.iter_mut() { x.1 += s0 as u32; }
        st
    };
    let counts: Vec<u32> = pieces.par_iter().fold(|| vec![0u32; 1 << B], |mut c, p| {
        for (h, _) in strobes_of(p) { c[(h >> (64 - B)) as usize] += 1; }
        c
    }).reduce(|| vec![0u32; 1 << B], |mut a, b| { for (x, y) in a.iter_mut().zip(b) { *x += y; } a });
    let mut bucket = vec![0u64; (1usize << B) + 1];
    for (i, &v) in counts.iter().enumerate() { bucket[i + 1] = v as u64; }
    drop(counts);
    for i in 1..bucket.len() { bucket[i] += bucket[i - 1]; }
    let total = bucket[1 << B] as usize;
    let mut entries = vec![0u64; total];
    {
        let mut fill: Vec<u64> = bucket[..1 << B].to_vec();
        for p in &pieces {
            for (h, pos) in strobes_of(p) {
                let bi = (h >> (64 - B)) as usize;
                entries[fill[bi] as usize] = (((h >> 8) & 0xffff_ffff) << 32) | pos as u64;
                fill[bi] += 1;
            }
        }
    }
    // je Fach sortieren
    {
        let mut rest: &mut [u64] = &mut entries;
        let mut slices = Vec::with_capacity(1 << B);
        for i in 0..(1usize << B) {
            let n = (bucket[i + 1] - bucket[i]) as usize;
            let (a, b) = rest.split_at_mut(n);
            slices.push(a); rest = b;
        }
        slices.par_iter_mut().for_each(|s| s.sort_unstable());
    }
    eprintln!("Index {name}: {} bp, {} Strobes, {:.1} s, {:.0} MB",
        genome.len(), total, t.elapsed().as_secs_f64(), (total * 8 + bucket.len() * 8) as f64 / 1e6);
    for (fi, fq) in args[2..].iter().enumerate() {
    // Reads
    let mut reads = Vec::new();
    let r = BufReader::new(std::fs::File::open(fq).unwrap());
    let mut it = r.lines();
    while let Some(Ok(h)) = it.next() {
        let s = it.next().unwrap().unwrap(); it.next(); it.next();
        reads.push((h, s.into_bytes()));
    }
    let t = Instant::now();
    let res: Vec<(i64, u32, u32, u32, u32)> = reads.par_iter().map_init(
        || (Vec::new(), Vec::new(), Vec::new(), Vec::new()),
        |(sy, st, hits, rc), (_, seq)| {
            rc.clear();
            rc.extend(seq.iter().rev().map(|&b| match b { b'A'=>b'T', b'C'=>b'G', b'G'=>b'C', b'T'=>b'A', x=>x }));
            hits.clear();
            let mut skipped = 0u32;
            for (strand, s) in [&seq[..], &rc[..]].into_iter().enumerate() {
                syncmers(s, sy);
                randstrobes(sy, st);
                for &(h, qp) in st.iter() {
                    let bi = (h >> (64 - B)) as usize;
                    let (lo, hi) = (bucket[bi] as usize, bucket[bi + 1] as usize);
                    let low = (h >> 8) & 0xffff_ffff;
                    let e = &entries[lo..hi];
                    let a = e.partition_point(|&k| (k >> 32) < low);
                    let b = a + e[a..].iter().take(maxocc + 1).take_while(|&&k| (k >> 32) == low).count();
                    if b - a > maxocc { skipped += 1; continue; }
                    for &k in &e[a..b] {
                        let p = (k & 0xffff_ffff) as i64;
                        hits.push(((strand as i64) << 40 | (p - qp as i64 + (1 << 32)), 1u32));
                    }
                }
            }
            hits.sort_unstable();
            // Gruppen gleicher Diagonale (±10), dann die besten Gruppen ohne Lücken prüfen
            let mut groups: Vec<(u32, i64)> = Vec::new();
            let mut i = 0;
            while i < hits.len() {
                let mut j = i + 1;
                while j < hits.len() && hits[j].0 - hits[j - 1].0 <= 10 { j += 1; }
                groups.push(((j - i) as u32, hits[(i + j) / 2].0));
                i = j;
            }
            if groups.is_empty() { return (0, 0, 0, skipped, 0); }
            groups.sort_unstable_by(|a, b| b.0.cmp(&a.0));
            let top = groups[0].0;
            let mut cand: Vec<(u32, i64)> = groups.iter().take(CANDS)
                .map(|g| {
                    let strand = g.1 >> 40; let p = (g.1 & ((1 << 40) - 1)) - (1 << 32);
                    let q = if strand == 0 { &seq[..] } else { &rc[..] };
                    // beste Verschiebung um ±3 (kleine Lücken grob)
                    let mut bestmm = u32::MAX;
                    for d in -3i64..=3 {
                        let st = p + d;
                        if st < 0 || st as usize + q.len() > genome.len() { continue; }
                        let r = &genome[st as usize..st as usize + q.len()];
                        let mm = q.iter().zip(r).filter(|(a, b)| a != b).count() as u32;
                        bestmm = bestmm.min(mm);
                    }
                    (bestmm, g.1)
                }).collect();
            cand.sort_unstable();
            let second = cand.get(1).map_or(u32::MAX, |c| c.0);
            (cand[0].1, cand[0].0.saturating_add(1), second, skipped, groups[0].0)
        }).collect();
    let secs = t.elapsed().as_secs_f64();
    let (mut ok, mut wrong, mut none, mut uniq_ok, mut uniq_wrong) = (0, 0, 0, 0, 0);
    for ((h, _), (best, n, second, _, _)) in reads.iter().zip(&res) {
        let f: Vec<&str> = h[1..].split(':').collect();
        let pos: i64 = starts[f[0]] + f[1].parse::<i64>().unwrap() - 1;
        if *n == 0 { none += 1; continue; }
        let p = (best & ((1 << 40) - 1)) - (1 << 32);
        let good = (p - pos).abs() <= 30;
        if good { ok += 1 } else { wrong += 1; if wrong < 8 { eprintln!("{h} -> {p} n={n} second={second} strand={}", best >> 40); } }
        if *second > *n + 3 { if good { uniq_ok += 1 } else { uniq_wrong += 1 } }
    }
    {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(format!("strobe{}.tsv", fi + 1)).unwrap());
        for ((h, _), (best, n, second, skipped, top)) in reads.iter().zip(&res) {
            let parts: Vec<&str> = h[1..].split(':').collect();
            let p = (best & ((1 << 40) - 1)) - (1 << 32) - starts[parts[0]] + 1;
            writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}", &h[1..], p, best >> 40, n.wrapping_sub(1), second, skipped, top).unwrap();
        }
    }
    let n = reads.len();
    println!("{n} Reads in {secs:.2} s = {:.2} µs je Read ({} Kerne) | richtig {ok} falsch {wrong} ohne {none} | eindeutig richtig {uniq_ok} falsch {uniq_wrong}",
        secs * 1e6 / n as f64 * rayon::current_num_threads() as f64, rayon::current_num_threads());
    }
}
