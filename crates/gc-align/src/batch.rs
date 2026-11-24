//! Einen Block von Reads auf mehreren Kernen ausrichten.
//!
//! Jeder Read und jedes Paar wird für sich gerechnet; die Ergebnisse kommen in der Reihenfolge
//! der Eingabe zurück. Die Ausgabe hängt deshalb nicht von der Zahl der Kerne ab.

use crate::pair::{InsertStats, PairedRead, estimate_insert, pair_lines};
use crate::sam::{AlignOptions, Line, align_read, find_regions, read_lines};
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
    parallel_map(reads.len(), threads, |i| {
        let r = &reads[i];
        let (regions, frac_rep) = align_read(index, r.seq, first_id + i as u64, opt);
        read_lines(index, r.name, r.seq, r.qual, &regions, frac_rep, opt)
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
    // Stufe 1: Bereiche je Read.
    let found = parallel_map(reads.len(), threads, |i| {
        find_regions(index, reads[i].seq, opt)
    });
    let (regions, frac_rep): (Vec<_>, Vec<_>) = found.into_iter().unzip();
    // Stufe 2: Insertgröße aus dem ganzen Block.
    let pes = estimate_insert(index, &regions, opt);
    // Stufe 3: je Paar Partnersuche, Paarbewertung, Zeilen.
    let lines = parallel_map(reads.len() / 2, threads, |i| {
        let mut s = [0, 1].map(|k| {
            let r = &reads[2 * i + k];
            PairedRead {
                name: r.name,
                read: r.seq,
                qual: r.qual,
                regions: regions[2 * i + k].clone(),
                frac_rep: frac_rep[2 * i + k],
            }
        });
        pair_lines(index, &pes, first_pair_id + i as u64, &mut s, opt).0
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
