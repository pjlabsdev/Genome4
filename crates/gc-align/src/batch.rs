//! Einen Block von Reads auf mehreren Kernen ausrichten.
//!
//! Jeder Read und jedes Paar wird für sich gerechnet; die Ergebnisse kommen in der Reihenfolge
//! der Eingabe zurück. Die Ausgabe hängt deshalb nicht von der Zahl der Kerne ab.

use crate::finalize::mark_primary;
use crate::pair::{
    InsertStats, PairedRead, estimate_insert, pair_lines_with, planned_rescues, rescue_align,
    rescue_align_two,
};
use crate::sam::{AlignOptions, Line, read_lines, regions_from_intervals};
use crate::seed::collect_intervals_many;
use gc_index::FmIndex;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Ein Read der Eingabe.
#[derive(Debug, Clone, Copy)]
pub struct ReadInput<'a> {
    pub name: &'a str,
    /// Basen als Codes 0..3, sonst N.
    pub seq: &'a [u8],
    /// Qualitäten als ASCII; leer heißt keine.
    pub qual: &'a [u8],
}

/// Rechnet `f(0), f(1), ..., f(n - 1)` auf `threads` Kernen und liefert die Ergebnisse in
/// dieser Reihenfolge. Die Arbeit wird in kleinen Portionen verteilt, damit ein teurer
/// Abschnitt nicht an einem Kern hängt.
pub fn parallel_map<T: Send>(n: usize, threads: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let threads = threads.clamp(1, n.max(1));
    if threads == 1 {
        return (0..n).map(f).collect();
    }
    // Anfangs große Portionen (wenig Abstimmung), gegen Ende kleinere: Die Portion ist ein
    // Achtel des Rests je Kern, zwischen MIN_PORTION und MAX_PORTION. So warten am Ende eines
    // Blocks nicht alle Kerne auf den einen, der zuletzt eine große Portion genommen hat (auf
    // Rechnern mit langsamen Effizienzkernen besonders spürbar). Die Ausgabe hängt nicht davon ab.
    const MIN_PORTION: usize = 4;
    const MAX_PORTION: usize = 64;
    let next = AtomicUsize::new(0);
    let mut parts: Vec<(usize, Vec<T>)> = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| {
                    let mut mine = Vec::new();
                    loop {
                        let seen = next.load(Ordering::Relaxed);
                        let portion = (n.saturating_sub(seen) / (8 * threads))
                            .clamp(MIN_PORTION, MAX_PORTION);
                        let from = next.fetch_add(portion, Ordering::Relaxed);
                        if from >= n {
                            break;
                        }
                        let to = (from + portion).min(n);
                        mine.push((from, (from..to).map(&f).collect::<Vec<T>>()));
                    }
                    mine
                })
            })
            .collect();
        workers
            .into_iter()
            .flat_map(|w| match w.join() {
                Ok(part) => part,
                Err(panic) => std::panic::resume_unwind(panic),
            })
            .collect()
    });
    parts.sort_unstable_by_key(|p| p.0);
    parts.into_iter().flat_map(|p| p.1).collect()
}

/// Richtet einzelne Reads aus. `first_id` ist die laufende Nummer des ersten Reads.
/// Je Read seine Ausgabezeilen.
pub fn align_reads(
    index: &FmIndex,
    reads: &[ReadInput<'_>],
    first_id: u64,
    opt: &AlignOptions,
    threads: usize,
) -> Vec<Vec<Line>> {
    let found = regions_of_all(index, reads, opt, threads);
    found
        .into_iter()
        .enumerate()
        .map(|(i, (mut regions, frac_rep))| {
            let r = &reads[i];
            mark_primary(
                &mut regions,
                first_id + i as u64,
                &opt.region.scoring,
                &opt.finalize,
            );
            read_lines(index, r.name, r.seq, r.qual, &regions, frac_rep, opt)
        })
        .collect()
}

/// Reads je Gruppe, deren Seed-Suchen verschränkt laufen ([`collect_intervals_many`]).
const SEED_GROUP: usize = 32;

/// Schnellweg, dann Strobemer-Weg (wenn ein Strobe-Index da ist); `None` heißt voller Weg.
fn quick_region(
    index: &FmIndex,
    read: &[u8],
    opt: &AlignOptions,
) -> Option<Vec<crate::region::Region>> {
    use crate::fast::{FastMode, fast_region};
    let fast = match opt.fast {
        FastMode::On | FastMode::Simple => fast_region(index, read, opt, opt.fast == FastMode::On),
        _ => None,
    };
    fast.map(|x| vec![x]).or_else(|| {
        let s = opt.strobe_index.as_ref()?;
        let r = crate::strobe::strobe_regions(index, s, read, opt);
        if r.is_none()
            && strobe_seed_mode()
            && crate::strobe::LAST_CLOSE.with(|c| c.get())
            && let Some(seeds) = crate::strobe::strobe_seeds(index, s, read, opt)
        {
            let regs = crate::sam::regions_from_seeds(index, read, &seeds, opt);
            return (!regs.is_empty()).then_some(regs);
        }
        r
    })
}

/// Versuch Strobe-Saat: Reads mit zu nahem Zweitem (ohne gesperrte Randstrobes) bekommen
/// ihre Seeds aus den Strobe-Treffern statt aus der SMEM-Suche (`GC_STROBE_SAAT=1`).
fn strobe_seed_mode() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_SAAT").is_ok_and(|v| v == "1"))
}

/// Bereiche aller Reads ([`crate::sam::find_regions`]), auf `threads` Kernen, die Seed-Suche in
/// Gruppen verschränkt. Das Ergebnis ist dasselbe wie Read für Read.
fn regions_of_all(
    index: &FmIndex,
    reads: &[ReadInput<'_>],
    opt: &AlignOptions,
    threads: usize,
) -> Vec<(Vec<crate::region::Region>, f32)> {
    use crate::fast::FastMode;
    // Prüfmodus und ausgeschaltetes Verschränken: Read für Read wie bisher.
    if !opt.seed.interleave || opt.fast == FastMode::Check {
        return parallel_map(reads.len(), threads, |i| {
            crate::sam::find_regions(index, reads[i].seq, opt)
        });
    }
    if opt.seed.regroup {
        // Erst der Schnellweg für alle Reads, dann die übrigen in vollen Gruppen.
        let mut out: Vec<Option<(Vec<crate::region::Region>, f32)>> =
            parallel_map(reads.len(), threads, |i| {
                quick_region(index, reads[i].seq, opt).map(|x| (x, 0.0))
            });
        let rest: Vec<usize> = (0..reads.len()).filter(|&i| out[i].is_none()).collect();
        let found = parallel_map(rest.len().div_ceil(SEED_GROUP), threads, |g| {
            let ids = &rest[g * SEED_GROUP..((g + 1) * SEED_GROUP).min(rest.len())];
            let seqs: Vec<&[u8]> = ids.iter().map(|&i| reads[i].seq).collect();
            let intervals = collect_intervals_many(index, &seqs, &opt.seed);
            ids.iter()
                .zip(&intervals)
                .map(|(&i, iv)| regions_from_intervals(index, reads[i].seq, iv, opt))
                .collect::<Vec<_>>()
        });
        for (&i, r) in rest.iter().zip(found.into_iter().flatten()) {
            out[i] = Some(r);
        }
        return out.into_iter().map(Option::unwrap).collect();
    }
    let groups = parallel_map(reads.len().div_ceil(SEED_GROUP), threads, |g| {
        let group = &reads[g * SEED_GROUP..((g + 1) * SEED_GROUP).min(reads.len())];
        // Erst der Schnellweg; nur die übrigen Reads gehen gemeinsam in die Seed-Suche.
        let mut out: Vec<Option<(Vec<crate::region::Region>, f32)>> = group
            .iter()
            .map(|r| quick_region(index, r.seq, opt).map(|x| (x, 0.0)))
            .collect();
        let rest: Vec<usize> = (0..group.len()).filter(|&i| out[i].is_none()).collect();
        let seqs: Vec<&[u8]> = rest.iter().map(|&i| group[i].seq).collect();
        let intervals = collect_intervals_many(index, &seqs, &opt.seed);
        for (&i, iv) in rest.iter().zip(&intervals) {
            out[i] = Some(regions_from_intervals(index, group[i].seq, iv, opt));
        }
        out.into_iter().map(Option::unwrap).collect::<Vec<_>>()
    });
    groups.into_iter().flatten().collect()
}

/// Paar geht den bisherigen Weg (sonst [`crate::paarweg::Kind`] als Zahl).
const NORMAL: u8 = 255;
/// Voller Weg, aber ohne Partnersuche.
const FULL_NO_RESCUE: u8 = 254;

/// Versuch Strobe-Paarweg: erst Schnell- und Strobeweg je Read; fehlt einem Read eines Paares
/// der Bereich, werden beide aus ihren Strobe-Ankern gemeinsam verortet
/// ([`crate::paarweg::pair_regions`]). Diese Paare (zweiter Wert `true`) bekommen keine
/// Partnersuche. Was übrig bleibt, geht den vollen Weg.
#[allow(clippy::type_complexity)]
fn regions_by_pairs(
    index: &FmIndex,
    reads: &[ReadInput<'_>],
    opt: &AlignOptions,
    threads: usize,
) -> (Vec<(Vec<crate::region::Region>, f32)>, Vec<u8>) {
    let strobes = opt.strobe_index.as_ref().expect("Strobe-Index");
    let n_pairs = reads.len() / 2;
    let per_pair = parallel_map(n_pairs, threads, |i| {
        let q = [0, 1].map(|k| quick_region(index, reads[2 * i + k].seq, opt));
        if q[0].is_some() && q[1].is_some() {
            return (q, NORMAL, [0.0; 2]);
        }
        match crate::paarweg::pair_regions(
            index,
            strobes,
            [reads[2 * i].seq, reads[2 * i + 1].seq],
            [q[0].is_some(), q[1].is_some()],
            opt,
        ) {
            Some((_, crate::paarweg::Kind::Unsure, _)) if crate::paarweg::unsure_full() => {
                (q, FULL_NO_RESCUE, [0.0; 2])
            }
            Some((p, kind, frac)) => {
                // Ein schon verorteter Read behält seinen Bereich, dazu die Orte seiner Anker
                // (statt der Partnersuche dort).
                let mut k = 0;
                let both = p.map(|mut pk| {
                    if let Some(qk) = &q[k] {
                        let mut m = qk.clone();
                        if crate::paarweg::merge_mode() {
                            pk.retain(|r| {
                                !m.iter().any(|o| {
                                    let ov = o.re.min(r.re) as i64 - o.rb.max(r.rb) as i64;
                                    ov * 100 >= (o.re - o.rb).min(r.re - r.rb) as i64 * 95
                                })
                            });
                            m.extend(pk);
                            m.sort_by_key(|r| std::cmp::Reverse(r.score));
                        }
                        pk = m;
                    }
                    k += 1;
                    Some(pk)
                });
                (both, kind as u8, frac)
            }
            None => (q, NORMAL, [0.0; 2]),
        }
    });
    let mut out: Vec<Option<(Vec<crate::region::Region>, f32)>> = Vec::with_capacity(reads.len());
    let mut no_rescue = Vec::with_capacity(n_pairs);
    for (q, nr, frac) in per_pair {
        no_rescue.push(nr);
        for (r, f) in q.into_iter().zip(frac) {
            out.push(r.map(|x| (x, f)));
        }
    }
    let rest: Vec<usize> = (0..reads.len()).filter(|&i| out[i].is_none()).collect();
    let found = parallel_map(rest.len().div_ceil(SEED_GROUP), threads, |g| {
        let ids = &rest[g * SEED_GROUP..((g + 1) * SEED_GROUP).min(rest.len())];
        let seqs: Vec<&[u8]> = ids.iter().map(|&i| reads[i].seq).collect();
        let intervals = collect_intervals_many(index, &seqs, &opt.seed);
        ids.iter()
            .zip(&intervals)
            .map(|(&i, iv)| regions_from_intervals(index, reads[i].seq, iv, opt))
            .collect::<Vec<_>>()
    });
    for (&i, r) in rest.iter().zip(found.into_iter().flatten()) {
        out[i] = Some(r);
    }
    (out.into_iter().map(Option::unwrap).collect(), no_rescue)
}

/// Die gesammelten Fenster der Partnersuche mit einem Rechenkern: Vorwärtslauf dort, Anfang des
/// Alignments auf dem Prozessor. Fenster, die der Kern nicht genau so rechnen kann (zu lange
/// Reads), und alle, wenn der Kern ausfällt, rechnet der Prozessor wie bisher.
fn rescue_with_kernel(
    index: &FmIndex,
    reads: &[ReadInput<'_>],
    jobs: &[(usize, usize, crate::pair::RescueWindow)],
    opt: &AlignOptions,
    threads: usize,
    kernel: &dyn crate::local::RescueKernel,
) -> Vec<(crate::pair::RescueWindow, crate::local::LocalResult)> {
    use crate::local::{from_rows, rows_exact, with_start};
    let sc = &opt.region.scoring;
    let min_sub = opt.finalize.min_seed_len * sc.matching;
    let prepared = parallel_map(jobs.len(), threads, |j| {
        let (i, k, w) = jobs[j];
        let mate = reads[2 * i + k].seq;
        let q = crate::pair::rescue_query(mate, w.is_rev);
        let t = index.text_range(w.rb as u64, w.re as u64);
        let on_kernel = rows_exact(q.len(), sc);
        (q, t, on_kernel)
    });
    let picked: Vec<usize> = (0..jobs.len()).filter(|&j| prepared[j].2).collect();
    let refs: Vec<(&[u8], &[u8])> = picked
        .iter()
        .map(|&j| (&prepared[j].0[..], &prepared[j].1[..]))
        .collect();
    let forward = kernel.forward(&refs, sc);
    let mut slot = vec![usize::MAX; jobs.len()];
    if forward.is_some() {
        for (n, &j) in picked.iter().enumerate() {
            slot[j] = n;
        }
    }
    parallel_map(jobs.len(), threads, |j| {
        let (i, k, w) = jobs[j];
        let (q, t, _) = &prepared[j];
        let r = match (&forward, slot[j]) {
            (Some(f), n) if n != usize::MAX => {
                let (gmax, te, qe, rows) = &f[n];
                with_start(
                    q,
                    t,
                    sc,
                    true,
                    min_sub,
                    from_rows(*gmax, *te, *qe, rows, sc, min_sub),
                )
            }
            _ => rescue_align(index, reads[2 * i + k].seq, &w, opt),
        };
        (w, r)
    })
}

/// Richtet einen Block von Paaren aus. `reads` enthält Read 1, Read 2, Read 1, ...;
/// `first_pair_id` ist die laufende Nummer des ersten Paares. Je Paar seine Ausgabezeilen,
/// dazu die für den Block geschätzte Insertgröße.
pub fn align_pairs(
    index: &FmIndex,
    reads: &[ReadInput<'_>],
    first_pair_id: u64,
    opt: &AlignOptions,
    threads: usize,
) -> (Vec<Vec<Line>>, [InsertStats; 4]) {
    assert!(
        reads.len().is_multiple_of(2),
        "Paare brauchen eine gerade Zahl von Reads"
    );
    // Stufe 1: Bereiche je Read (im Versuch Strobe-Paarweg teils je Paar, dann ohne
    // Partnersuche).
    let (found, no_rescue) = if crate::paarweg::mode() && opt.strobe_index.is_some() {
        regions_by_pairs(index, reads, opt, threads)
    } else {
        (regions_of_all(index, reads, opt, threads), Vec::new())
    };
    let (regions, frac_rep): (Vec<_>, Vec<_>) = found.into_iter().unzip();
    use crate::paarweg::Kind;
    let mut opt_nr = opt.clone();
    opt_nr.pair.max_matesw = 0;
    let mut opt_un = opt_nr.clone();
    opt_un.pair.mapq_cap = crate::paarweg::unsure_cap().unwrap_or(255);
    let kind = |i: usize| no_rescue.get(i).copied().unwrap_or(NORMAL);
    let skip = |i: usize| {
        kind(i) == Kind::Sure as u8 || kind(i) == Kind::Unsure as u8 || kind(i) == FULL_NO_RESCUE
    };
    // Stufe 2: Insertgröße aus dem ganzen Block.
    let pes = estimate_insert(index, &regions, opt);
    // Stufe 3: je Paar Partnersuche, Paarbewertung, Zeilen. Die lokalen Ausrichtungen der
    // Partnersuche werden vorab für den ganzen Block gesammelt und gemeinsam gerechnet.
    let paired = |i: usize| {
        [0, 1].map(|k| {
            let r = &reads[2 * i + k];
            PairedRead {
                name: r.name,
                read: r.seq,
                qual: r.qual,
                regions: regions[2 * i + k].clone(),
                frac_rep: frac_rep[2 * i + k],
            }
        })
    };
    let n_pairs = reads.len() / 2;
    let planned = if opt.pair.batch_rescue {
        parallel_map(n_pairs, threads, |i| {
            if skip(i) {
                [Vec::new(), Vec::new()]
            } else {
                planned_rescues(index, &pes, &paired(i), opt)
            }
        })
    } else {
        Vec::new()
    };
    // Alle Fenster des Blocks in einer Liste: (Paar, gesuchter Read, Fenster).
    let jobs: Vec<(usize, usize, crate::pair::RescueWindow)> = planned
        .iter()
        .enumerate()
        .flat_map(|(i, p)| (0..2).flat_map(move |k| p[k].iter().map(move |w| (i, k, *w))))
        .collect();
    let results = match &opt.rescue_kernel {
        Some(kernel) => rescue_with_kernel(index, reads, &jobs, opt, threads, kernel.as_ref()),
        None if opt.pair.rescue_two => {
            let two = parallel_map(jobs.len().div_ceil(2), threads, |g| {
                let (i, k, w) = jobs[2 * g];
                let a = reads[2 * i + k].seq;
                match jobs.get(2 * g + 1) {
                    Some(&(i2, k2, w2)) => {
                        let b = reads[2 * i2 + k2].seq;
                        let (ra, rb) = rescue_align_two(index, (a, &w), (b, &w2), opt);
                        vec![(w, ra), (w2, rb)]
                    }
                    None => vec![(w, rescue_align(index, a, &w, opt))],
                }
            });
            two.into_iter().flatten().collect()
        }
        None => parallel_map(jobs.len(), threads, |j| {
            let (i, k, w) = jobs[j];
            (w, rescue_align(index, reads[2 * i + k].seq, &w, opt))
        }),
    };
    // Je Paar und Read der Abschnitt seiner Ergebnisse.
    let mut spans = vec![[(0usize, 0usize); 2]; if planned.is_empty() { 0 } else { n_pairs }];
    let mut at = 0;
    for (i, p) in planned.iter().enumerate() {
        for k in 0..2 {
            spans[i][k] = (at, at + p[k].len());
            at += p[k].len();
        }
    }
    let lines = parallel_map(n_pairs, threads, |i| {
        let mut s = paired(i);
        let done = if spans.is_empty() {
            [&[][..], &[][..]]
        } else {
            [0, 1].map(|k| &results[spans[i][k].0..spans[i][k].1])
        };
        let o = match kind(i) {
            k if k == Kind::Sure as u8 || k == FULL_NO_RESCUE => &opt_nr,
            k if k == Kind::Unsure as u8 => &opt_un,
            _ => opt,
        };
        pair_lines_with(index, &pes, first_pair_id + i as u64, &mut s, o, done).0
    });
    (lines, pes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn results_come_back_in_input_order_for_any_thread_count() {
        let want: Vec<usize> = (0..1000).map(|i| i * i).collect();
        for threads in [1, 2, 3, 8, 64] {
            assert_eq!(parallel_map(1000, threads, |i| i * i), want);
        }
        assert!(parallel_map(0, 4, |i| i).is_empty());
    }

    #[test]
    #[should_panic(expected = "absichtlich")]
    fn a_panic_in_a_worker_reaches_the_caller() {
        parallel_map(500, 4, |i| {
            assert!(i != 321, "absichtlich");
            i
        });
    }
}
