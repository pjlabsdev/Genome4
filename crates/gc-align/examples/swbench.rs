use gc_align::extend::Scoring;
use gc_align::local::{local_align, local_align_many, local_align_pair};
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
    let mode = std::env::var("PAIR").unwrap_or_default();
    if mode == "4" || mode == "3" || mode == "2m" {
        let n = if mode == "4" {
            4
        } else if mode == "3" {
            3
        } else {
            2
        };
        for p in jobs.chunks(n) {
            let s: i64 = match p {
                [a, b, c, d] => local_align_many(
                    [
                        (&a.0[..], &a.1[..]),
                        (&b.0, &b.1),
                        (&c.0, &c.1),
                        (&d.0, &d.1),
                    ],
                    &sc,
                    true,
                    19,
                    true,
                )
                .iter()
                .map(|r| r.score as i64)
                .sum(),
                [a, b, c] => local_align_many(
                    [(&a.0[..], &a.1[..]), (&b.0, &b.1), (&c.0, &c.1)],
                    &sc,
                    true,
                    19,
                    true,
                )
                .iter()
                .map(|r| r.score as i64)
                .sum(),
                [a, b] => {
                    local_align_many([(&a.0[..], &a.1[..]), (&b.0, &b.1)], &sc, true, 19, true)
                        .iter()
                        .map(|r| r.score as i64)
                        .sum()
                }
                _ => p
                    .iter()
                    .map(|(q, t)| local_align(q, t, &sc, true, 19, true).score as i64)
                    .sum(),
            };
            acc += s;
        }
    } else if mode == "1" {
        for p in jobs.chunks(2) {
            let (a, b) =
                local_align_pair((&p[0].0, &p[0].1), (&p[1].0, &p[1].1), &sc, true, 19, true);
            acc += a.score as i64 + b.score as i64;
        }
    } else {
        for (q, tt) in &jobs {
            acc += local_align(q, tt, &sc, true, 19, true).score as i64;
        }
    }
    let e = t.elapsed().as_secs_f64();
    println!(
        "{:.1} us/call, {:.2} ns/cell, acc {}",
        e / 2000.0 * 1e6,
        e / 2000.0 / (148.0 * tlen as f64) * 1e9,
        acc
    );
}
