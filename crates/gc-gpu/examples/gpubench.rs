//! Grafikeinheit gegen Prozessor für die lokalen Ausrichtungen der Partnersuche.
//!
//! `cargo run --release -p gc-gpu --example gpubench -- [Fenster] [Länge des Fensters]`
//!
//! Erzeugt zufällige Fenster (ein Drittel mit dem Read darin, mit Abweichungen und Lücken),
//! rechnet sie auf der Grafikeinheit und auf allen Kernen des Prozessors, vergleicht jedes
//! Ergebnis und gibt die Zeiten aus. Beendet sich mit Fehler, wenn ein Ergebnis abweicht.
use gc_align::Scoring;
use gc_align::local::{from_rows, local_align};

fn main() {
    let mut args = std::env::args().skip(1);
    let n: usize = args.next().map_or(20_000, |a| a.parse().unwrap());
    let tlen: usize = args.next().map_or(1320, |a| a.parse().unwrap());
    let mut s = 0x5eed_0bad_cafe_f00du64;
    let mut rng = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    let mut jobs: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(n);
    for k in 0..n {
        let t: Vec<u8> = (0..tlen).map(|_| (rng() % 4) as u8).collect();
        let mut q: Vec<u8> = if k % 3 == 0 {
            let at = (rng() % (tlen as u64 - 160)) as usize;
            t[at..at + 148].to_vec()
        } else {
            (0..148).map(|_| (rng() % 4) as u8).collect()
        };
        for _ in 0..rng() % 6 {
            let p = (rng() % q.len() as u64) as usize;
            match rng() % 4 {
                0 => q[p] = 4,
                1 => {
                    q.remove(p);
                }
                2 => q.insert(p, (rng() % 4) as u8),
                _ => q[p] = (q[p] + 1) & 3,
            }
        }
        jobs.push((q, t));
    }
    let sc = Scoring::default();
    let ks = gc_gpu::KernelScoring {
        matching: sc.matching as u32,
        mismatch: sc.mismatch as u32,
        ambiguous: sc.ambiguous as u32,
        o_del: sc.o_del as u32,
        e_del: sc.e_del as u32,
        o_ins: sc.o_ins as u32,
        e_ins: sc.e_ins as u32,
    };
    let gpu = gc_gpu::Gpu::open().expect("Grafikeinheit");
    println!("Grafikeinheit: {}", gpu.name);
    let refs: Vec<(&[u8], &[u8])> = jobs.iter().map(|(q, t)| (&q[..], &t[..])).collect();
    // Einmal zum Aufwärmen (Übersetzen des Kerns), dann gemessen.
    let _ = gpu
        .forward(&refs[..refs.len().min(256)], &ks)
        .expect("Lauf");
    let t0 = std::time::Instant::now();
    let fw = gpu.forward(&refs, &ks).expect("Lauf");
    let gpu_s = t0.elapsed().as_secs_f64();
    let cores = std::thread::available_parallelism().map_or(1, |c| c.get());
    let t0 = std::time::Instant::now();
    let chunk = n.div_ceil(cores);
    let cpu: Vec<_> = std::thread::scope(|sc_| {
        let hs: Vec<_> = jobs
            .chunks(chunk)
            .map(|c| {
                sc_.spawn(move || {
                    c.iter()
                        .map(|(q, t)| local_align(q, t, &Scoring::default(), true, 19, false))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    let cpu_s = t0.elapsed().as_secs_f64();
    let mut bad = 0;
    for (k, ((g, te, qe, rows), want)) in fw.iter().zip(&cpu).enumerate() {
        let got = from_rows(*g, *te, *qe, rows, &sc, 19);
        if got != *want {
            if bad < 5 {
                eprintln!("Fenster {k}: Grafikeinheit {got:?}, Prozessor {want:?}");
            }
            bad += 1;
        }
    }
    println!(
        "{n} Fenster zu {tlen} Basen: Grafikeinheit {:.3} s ({:.1} µs je Fenster), Prozessor {cores} Kerne {:.3} s ({:.1} µs je Fenster), abweichend {bad}",
        gpu_s,
        gpu_s / n as f64 * 1e6,
        cpu_s,
        cpu_s / n as f64 * 1e6
    );
    if bad > 0 {
        std::process::exit(1);
    }
}
