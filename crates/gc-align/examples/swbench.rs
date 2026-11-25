use gc_align::extend::Scoring;
use gc_align::local::local_align;
fn main() {
    let mut s = 88172645463325252u64;
    let mut r = |n: usize| -> Vec<u8> {
        (0..n)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                (s % 4) as u8
            })
            .collect()
    };
    let sc = Scoring::default();
    let tlen: usize = std::env::args().nth(1).map_or(1000, |a| a.parse().unwrap());
    let jobs: Vec<(Vec<u8>, Vec<u8>)> = (0..2000).map(|_| (r(148), r(tlen))).collect();
    let t = std::time::Instant::now();
    let mut acc = 0i64;
    for (q, tt) in &jobs {
        acc += local_align(q, tt, &sc, true, 19, true).score as i64;
    }
    let e = t.elapsed().as_secs_f64();
    println!(
        "{:.1} us/call, {:.2} ns/cell, acc {}",
        e / 2000.0 * 1e6,
        e / 2000.0 / (148.0 * tlen as f64) * 1e9,
        acc
    );
}
