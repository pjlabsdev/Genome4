//! Blockweiser Bau des FM-Index mit festem RAM-Bedarf, für das ganze Genom.
//!
//! Das volle Suffix-Array des Humangenoms (beide Stränge, 6,2 Mrd. Zeilen) passt nicht in den
//! Speicher. Deshalb werden die Suffixe nach ihren ersten `k` Basen in Fächer geteilt, benachbarte
//! Fächer zu Blöcken fester Größe zusammengefasst und Block für Block sortiert. Die Blöcke kommen
//! in lexikographischer Reihenfolge, die BWT entsteht also fortlaufend.
//!
//! Damit das Sortieren auch in langen Wiederholungen schnell bleibt, werden zwei Suffixe höchstens
//! `v` Basen weit direkt verglichen. Sind sie dann noch gleich, entscheidet eine vorab berechnete
//! Rangliste über eine Stichprobe der Positionen (Differenzüberdeckung nach Kärkkäinen und
//! Burkhardt 2003; für den blockweisen Bau Kärkkäinen 2007).

use std::cmp::Ordering;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

use serde_json::json;

use crate::fm::{BuildError, Builder, FmIndex, SENTINEL, Variant};
use crate::reference::Reference;

/// Einstellungen für den blockweisen Bau.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// Kerne (0 = alle).
    pub threads: usize,
    /// Höchstzahl Suffixe je Block; bestimmt den RAM-Bedarf (16 Bytes je Suffix).
    pub block_suffixes: u64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            threads: 0,
            block_suffixes: 96_000_000,
        }
    }
}

/// Innere Stellgrößen; in Tests klein gewählt, damit alle Wege durchlaufen werden.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Params {
    /// Periode der Differenzüberdeckung und zugleich die Tiefe des direkten Vergleichs.
    pub v: u64,
    /// Länge des Präfixes, nach dem in Fächer geteilt wird (höchstens 12).
    pub k: u32,
    pub block_suffixes: u64,
    pub threads: usize,
}

/// Beide Stränge als 2-Bit-Folge, erste Base eines Wortes in den höchsten Bits. So ordnet der
/// Vergleich zweier Wörter 32 Basen auf einmal lexikographisch.
struct Text {
    n: u64,
    words: Vec<u64>,
}

impl Text {
    fn from_reference(reference: &Reference) -> Text {
        let l = reference.len();
        let n = 2 * l;
        let mut words = vec![0u64; (n / 32) as usize + 2];
        let mut put = |i: u64, code: u8| {
            words[(i / 32) as usize] |= u64::from(code) << (62 - 2 * (i % 32));
        };
        for i in 0..l {
            let c = reference.forward.get(i);
            put(i, c);
            put(n - 1 - i, 3 - c);
        }
        Text { n, words }
    }

    #[inline]
    fn get(&self, i: u64) -> u8 {
        (self.words[(i / 32) as usize] >> (62 - 2 * (i % 32)) & 3) as u8
    }

    /// 32 Basen ab Position `i`; hinter dem Textende stehen Nullen.
    #[inline]
    fn word(&self, i: u64) -> u64 {
        let k = (i / 32) as usize;
        let r = (i % 32) as u32;
        if r == 0 {
            self.words[k]
        } else {
            self.words[k] << (2 * r) | self.words[k + 1] >> (64 - 2 * r)
        }
    }

    /// Die ersten `k` Basen ab `i` als Zahl (hinter dem Textende mit A aufgefüllt).
    #[inline]
    fn kmer(&self, i: u64, k: u32) -> usize {
        (self.word(i) >> (64 - 2 * k)) as usize
    }

    /// Vergleicht die Suffixe ab `i` und `j` höchstens `depth` Basen weit. `Equal` heißt: beide
    /// haben mindestens `depth` Basen und stimmen darin überein. Endet ein Suffix vorher, ist das
    /// kürzere das kleinere.
    #[inline]
    fn cmp_prefix(&self, i: u64, j: u64, depth: u64) -> Ordering {
        let limit = (self.n - i).min(self.n - j).min(depth);
        let mut off = 0;
        while off < limit {
            let (mut a, mut b) = (self.word(i + off), self.word(j + off));
            let rem = limit - off;
            if rem < 32 {
                let s = 64 - 2 * rem as u32;
                a >>= s;
                b >>= s;
            }
            if a != b {
                return a.cmp(&b);
            }
            off += 32;
        }
        if limit == depth {
            Ordering::Equal
        } else {
            (self.n - i).cmp(&(self.n - j))
        }
    }
}

/// Differenzüberdeckung modulo `v`: eine kleine Menge D, sodass es zu jedem Abstand r ein a in D
/// gibt mit a + r (mod v) ebenfalls in D.
struct Cover {
    v: u64,
    /// Elemente von D, aufsteigend.
    d: Vec<u64>,
    /// Rest → Index in `d` oder `u32::MAX`.
    index_of: Vec<u32>,
    /// Abstand r → ein a aus D mit (a + r) mod v in D.
    shift: Vec<u32>,
}

impl Cover {
    fn new(v: u64) -> Cover {
        // Konstruktion nach Colbourn und Ling: 6r + 4 Elemente decken alle Abstände modulo
        // 24r² + 36r + 13. Falls sie für dieses v nicht aufgeht, wird gierig ergänzt.
        let mut r = 0u64;
        while 24 * r * r + 36 * r + 13 < v {
            r += 1;
        }
        let steps = [
            (r, 1),
            (1, r + 1),
            (r, 2 * r + 1),
            (2 * r + 1, 4 * r + 3),
            (r + 1, 2 * r + 2),
            (r, 1),
        ];
        let mut d = vec![0u64];
        let mut at = 0;
        for (times, step) in steps {
            for _ in 0..times {
                at += step;
                d.push(at % v);
            }
        }
        d.sort_unstable();
        d.dedup();
        let covered = |d: &[u64]| -> Vec<bool> {
            let mut seen = vec![false; v as usize];
            for &a in d {
                for &b in d {
                    seen[((b + v - a) % v) as usize] = true;
                }
            }
            seen
        };
        let mut seen = covered(&d);
        while let Some(missing) = seen.iter().position(|s| !s) {
            // gierig: ein Element ergänzen, das den fehlenden Abstand zu 0 herstellt
            d.push(missing as u64 % v);
            d.sort_unstable();
            d.dedup();
            seen = covered(&d);
        }
        let mut index_of = vec![u32::MAX; v as usize];
        for (i, &a) in d.iter().enumerate() {
            index_of[a as usize] = i as u32;
        }
        let mut shift = vec![u32::MAX; v as usize];
        for &a in &d {
            for &b in &d {
                let diff = ((b + v - a) % v) as usize;
                if shift[diff] == u32::MAX {
                    shift[diff] = a as u32;
                }
            }
        }
        Cover {
            v,
            d,
            index_of,
            shift,
        }
    }

    fn size(&self) -> u64 {
        self.d.len() as u64
    }

    /// Platz der Stichprobenposition `pos` in der Rangliste.
    #[inline]
    fn slot(&self, pos: u64) -> usize {
        (pos / self.v * self.size() + u64::from(self.index_of[(pos % self.v) as usize])) as usize
    }

    #[inline]
    fn pos_of_slot(&self, slot: u64) -> u64 {
        slot / self.size() * self.v + self.d[(slot % self.size()) as usize]
    }

    /// Kleinster Versatz d < v, bei dem `i + d` und `j + d` beide Stichprobenpositionen sind.
    #[inline]
    fn common_offset(&self, i: u64, j: u64) -> u64 {
        let v = self.v;
        let diff = (j % v + v - i % v) % v;
        let a = u64::from(self.shift[diff as usize]);
        (a + v - i % v) % v
    }
}

fn thread_count(requested: usize) -> usize {
    if requested > 0 {
        requested
    } else {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    }
}

/// Teilt `0..len` in höchstens `parts` Bereiche, deren Grenzen auf Fachgrenzen (`starts`) liegen.
/// Liefert Fachindizes `[b0, b1, ...]`.
fn bucket_cuts(starts: &[usize], parts: usize) -> Vec<usize> {
    let buckets = starts.len() - 1;
    let total = starts[buckets];
    let mut cuts = vec![0usize];
    for p in 1..parts {
        let target = total / parts * p;
        let b = starts.partition_point(|&s| s < target).min(buckets);
        if b > *cuts.last().expect("nicht leer") {
            cuts.push(b);
        }
    }
    if *cuts.last().expect("nicht leer") != buckets {
        cuts.push(buckets);
    }
    cuts
}

/// Sortiert jedes Fach von `data` für sich, verteilt auf mehrere Kerne.
fn sort_buckets<T: Send>(
    data: &mut [T],
    starts: &[usize],
    threads: usize,
    cmp: impl Fn(&T, &T) -> Ordering + Sync,
) {
    let cuts = bucket_cuts(starts, threads);
    std::thread::scope(|scope| {
        let mut rest = data;
        let mut consumed = 0;
        for w in cuts.windows(2) {
            let (b0, b1) = (w[0], w[1]);
            let (mine, tail) = rest.split_at_mut(starts[b1] - consumed);
            rest = tail;
            let base = consumed;
            consumed = starts[b1];
            let cmp = &cmp;
            scope.spawn(move || {
                for b in b0..b1 {
                    let slice = &mut mine[starts[b] - base..starts[b + 1] - base];
                    if slice.len() > 1 {
                        slice.sort_unstable_by(cmp);
                    }
                }
            });
        }
    });
}

/// Rangliste aller Stichproben-Suffixe: `rank[slot]` ist für jede gültige Stichprobenposition
/// eindeutig und ordnet wie die Suffixe; 0 steht für "hinter dem Textende".
fn sample_ranks(text: &Text, cover: &Cover, k: u32, threads: usize) -> Vec<u32> {
    let (v, m, n) = (cover.v, cover.size(), text.n);
    let slots = n.div_ceil(v) * m;
    assert!(
        slots < u64::from(u32::MAX),
        "Text zu lang für diese Stichprobendichte"
    );

    // Nach den ersten k Basen in Fächer legen (Zählsortieren).
    let n_buckets = 1usize << (2 * k);
    let mut starts = vec![0usize; n_buckets + 1];
    let valid = |slot: u64| cover.pos_of_slot(slot) < n;
    for slot in (0..slots).filter(|&s| valid(s)) {
        starts[text.kmer(cover.pos_of_slot(slot), k) + 1] += 1;
    }
    for b in 0..n_buckets {
        starts[b + 1] += starts[b];
    }
    let mut order = vec![0u32; starts[n_buckets]];
    let mut next = starts.clone();
    for slot in (0..slots).filter(|&s| valid(s)) {
        let b = text.kmer(cover.pos_of_slot(slot), k);
        order[next[b]] = slot as u32;
        next[b] += 1;
    }
    drop(next);

    // Jedes Fach nach den ersten v Basen sortieren.
    sort_buckets(&mut order, &starts, threads, |&a, &b| {
        text.cmp_prefix(
            cover.pos_of_slot(u64::from(a)),
            cover.pos_of_slot(u64::from(b)),
            v,
        )
    });
    drop(starts);

    // Ränge vergeben; Gruppen mit gleichen ersten v Basen bleiben vorerst gleichrangig.
    let mut rank = vec![0u32; slots as usize];
    let mut groups: Vec<(usize, usize)> = Vec::new();
    let mut group_start = 0;
    for i in 0..order.len() {
        if i > 0 {
            let (a, b) = (u64::from(order[i - 1]), u64::from(order[i]));
            if text.cmp_prefix(cover.pos_of_slot(a), cover.pos_of_slot(b), v) != Ordering::Equal {
                if i - group_start > 1 {
                    groups.push((group_start, i));
                }
                group_start = i;
            }
        }
        rank[order[i] as usize] = group_start as u32 + 1;
    }
    if order.len() - group_start > 1 {
        groups.push((group_start, order.len()));
    }

    // Verdoppeln (Larsson und Sadakane): eine Gruppe ist nach h Basen gleich; der Rang der
    // Position h weiter (wieder eine Stichprobe, weil h ein Vielfaches von v ist) trennt sie
    // nach 2h Basen.
    let mut h_slots = m;
    let mut round = 0;
    let mut keyed: Vec<(u32, u32)> = Vec::new();
    while !groups.is_empty() {
        round += 1;
        gc_log::debug(
            "INDEX.STICHPROBE.RUNDE",
            "Stichproben-Ränge verfeinern",
            json!({ "runde": round, "offene_gruppen": groups.len(), "tiefe_basen": h_slots / m * v }),
        );
        let mut next_groups = Vec::new();
        for &(s, e) in &groups {
            keyed.clear();
            keyed.extend(order[s..e].iter().map(|&slot| {
                let ahead = u64::from(slot) + h_slots;
                (
                    if ahead < slots {
                        rank[ahead as usize]
                    } else {
                        0
                    },
                    slot,
                )
            }));
            keyed.sort_unstable();
            let mut sub_start = s;
            for (i, &(key, slot)) in keyed.iter().enumerate() {
                let at = s + i;
                if i > 0 && key != keyed[i - 1].0 {
                    if at - sub_start > 1 {
                        next_groups.push((sub_start, at));
                    }
                    sub_start = at;
                }
                order[at] = slot;
                rank[slot as usize] = sub_start as u32 + 1;
            }
            if e - sub_start > 1 {
                next_groups.push((sub_start, e));
            }
        }
        groups = next_groups;
        h_slots *= 2;
    }
    rank
}

#[derive(Clone, Copy)]
#[repr(C)]
struct Item {
    /// Die 32 Basen nach dem Fach-Präfix: entscheidet die meisten Vergleiche ohne Textzugriff.
    key: u64,
    pos_lo: u32,
    pos_hi: u8,
    /// Base vor dem Suffix (oder `SENTINEL`).
    prev: u8,
}

impl Item {
    #[inline]
    fn pos(&self) -> u64 {
        u64::from(self.pos_lo) | u64::from(self.pos_hi) << 32
    }
}

/// Zeiger auf den Block, in den mehrere Kerne an getrennten Stellen schreiben.
struct SharedItems(*mut Item);
// SAFETY: Jeder Kern schreibt nur in Plätze, die ihm über die Zählung vorab exklusiv
// zugeteilt wurden (siehe `fill_block`); kein Platz wird zweimal vergeben.
unsafe impl Sync for SharedItems {}

/// Legt alle Suffixe, deren Präfix in den Fächern `b0..b1` liegt, nach Fach geordnet in `items`.
/// `offsets[b - b0]` ist der Anfang des Fachs im Block.
fn fill_block(
    text: &Text,
    k: u32,
    b0: usize,
    b1: usize,
    offsets: &[usize],
    items: &mut Vec<Item>,
    threads: usize,
) {
    let n = text.n;
    let total = offsets[b1 - b0];
    items.clear();
    items.resize(
        total,
        Item {
            key: 0,
            pos_lo: 0,
            pos_hi: 0,
            prev: 0,
        },
    );
    let ranges: Vec<(u64, u64)> = (0..threads as u64)
        .map(|t| {
            (
                n / threads as u64 * t,
                if t + 1 == threads as u64 {
                    n
                } else {
                    n / threads as u64 * (t + 1)
                },
            )
        })
        .collect();

    // Durchgang 1: je Kern zählen, wie viele seiner Positionen in jedes Fach fallen.
    let counts: Vec<Vec<u32>> = std::thread::scope(|scope| {
        let handles: Vec<_> = ranges
            .iter()
            .map(|&(from, to)| {
                scope.spawn(move || {
                    let mut local = vec![0u32; b1 - b0];
                    for pos in from..to {
                        let b = text.kmer(pos, k);
                        if b >= b0 && b < b1 {
                            local[b - b0] += 1;
                        }
                    }
                    local
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("Zählkern"))
            .collect()
    });

    // Plätze zuteilen: innerhalb eines Fachs kommt Kern 0 zuerst, dann Kern 1 und so weiter.
    let mut cursors: Vec<Vec<usize>> = Vec::with_capacity(threads);
    let mut running: Vec<usize> = offsets[..b1 - b0].to_vec();
    for local in &counts {
        cursors.push(running.clone());
        for (r, &c) in running.iter_mut().zip(local) {
            *r += c as usize;
        }
    }
    debug_assert!(running.iter().zip(&offsets[1..]).all(|(a, b)| a == b));
    drop(counts);

    // Durchgang 2: eintragen.
    let shared = SharedItems(items.as_mut_ptr());
    let shared = &shared;
    std::thread::scope(|scope| {
        for (&(from, to), mut cursor) in ranges.iter().zip(cursors) {
            scope.spawn(move || {
                for pos in from..to {
                    let b = text.kmer(pos, k);
                    if b >= b0 && b < b1 {
                        let at = cursor[b - b0];
                        cursor[b - b0] += 1;
                        let item = Item {
                            key: if pos + u64::from(k) < n {
                                text.word(pos + u64::from(k))
                            } else {
                                0
                            },
                            pos_lo: pos as u32,
                            pos_hi: (pos >> 32) as u8,
                            prev: if pos == 0 {
                                SENTINEL
                            } else {
                                text.get(pos - 1)
                            },
                        };
                        debug_assert!(at < total);
                        // SAFETY: `at` liegt im Block und gehört nach der Zuteilung oben nur diesem Kern.
                        unsafe { shared.0.add(at).write(item) };
                    }
                }
            });
        }
    });
}

pub(crate) fn build_with(
    reference: Reference,
    variant: Variant,
    params: Params,
) -> Result<FmIndex, BuildError> {
    let l = reference.len();
    if l == 0 {
        return Err(BuildError::Empty);
    }
    let max = (1u64 << 39) - 1;
    if l > max {
        return Err(BuildError::TooLong { bases: l, max });
    }
    let (occ_interval, sa_interval) = variant.parameters();
    let Params {
        v,
        k,
        block_suffixes,
        threads,
    } = params;
    assert!((1..=12).contains(&k) && v >= 4);

    let st = gc_log::stage("index.text", json!({ "verfahren": "blockweise" }));
    let text = Text::from_reference(&reference);
    let n = text.n;
    st.finish(json!({ "zeichen": n, "bytes": text.words.len() * 8 }));

    let st = gc_log::stage("index.stichprobe", json!({ "periode": v }));
    let cover = Cover::new(v);
    let rank = sample_ranks(&text, &cover, k, threads);
    st.finish(json!({
        "ueberdeckung": cover.size(), "plaetze": rank.len(), "bytes": rank.len() * 4,
    }));

    // Fachgrößen über alle Suffixe.
    let st = gc_log::stage("index.faecher", json!({ "praefix_basen": k }));
    let n_buckets = 1usize << (2 * k);
    let per_thread: Vec<Vec<u32>> = std::thread::scope(|scope| {
        let text = &text;
        let handles: Vec<_> = (0..threads as u64)
            .map(|t| {
                let from = n / threads as u64 * t;
                let to = if t + 1 == threads as u64 {
                    n
                } else {
                    n / threads as u64 * (t + 1)
                };
                scope.spawn(move || {
                    let mut local = vec![0u32; n_buckets];
                    for pos in from..to {
                        local[text.kmer(pos, k)] += 1;
                    }
                    local
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("Zählkern"))
            .collect()
    });
    let mut sizes = vec![0u64; n_buckets];
    for local in &per_thread {
        for (s, &c) in sizes.iter_mut().zip(local) {
            *s += u64::from(c);
        }
    }
    drop(per_thread);

    // Benachbarte Fächer zu Blöcken zusammenfassen.
    let mut blocks: Vec<(usize, usize)> = Vec::new();
    let (mut b0, mut in_block) = (0usize, 0u64);
    for (b, &size) in sizes.iter().enumerate() {
        if in_block > 0 && in_block + size > block_suffixes {
            blocks.push((b0, b));
            b0 = b;
            in_block = 0;
        }
        in_block += size;
    }
    blocks.push((b0, n_buckets));
    let largest = sizes.iter().copied().max().unwrap_or(0);
    st.finish(json!({
        "faecher": n_buckets, "groesstes_fach": largest, "bloecke": blocks.len(),
        "suffixe_je_block_hoechstens": block_suffixes,
    }));

    let mut st = gc_log::stage(
        "index.bloecke",
        json!({ "bloecke": blocks.len(), "kerne": threads }),
    );
    let mut builder = Builder::new(n + 1, occ_interval, sa_interval);
    builder.push(text.get(n - 1), n); // Zeile 0: das Suffix, das nur aus dem Endzeichen besteht
    let mut items: Vec<Item> = Vec::new();
    let mut offsets: Vec<usize> = Vec::new();
    let mut done = 0u64;
    let deep_ties = AtomicUsize::new(0);
    for (bi, &(b0, b1)) in blocks.iter().enumerate() {
        offsets.clear();
        offsets.push(0);
        for &size in &sizes[b0..b1] {
            offsets.push(offsets.last().expect("nicht leer") + size as usize);
        }
        fill_block(&text, k, b0, b1, &offsets, &mut items, threads);
        sort_buckets(&mut items, &offsets, threads, |a, b| {
            a.key.cmp(&b.key).then_with(|| {
                let (i, j) = (a.pos(), b.pos());
                text.cmp_prefix(i, j, v).then_with(|| {
                    // Nach v Basen noch gleich: die Rangliste der Stichprobe entscheidet.
                    deep_ties.fetch_add(1, AtomicOrdering::Relaxed);
                    let d = cover.common_offset(i, j);
                    rank[cover.slot(i + d)].cmp(&rank[cover.slot(j + d)])
                })
            })
        });
        for item in &items {
            builder.push(item.prev, item.pos());
        }
        done += items.len() as u64;
        st.heartbeat(10, || {
            json!({ "block": bi + 1, "bloecke": blocks.len(), "suffixe_fertig": done, "anteil": done as f64 / n as f64 })
        });
    }
    debug_assert_eq!(done, n);
    let index = builder.finish(reference);
    st.finish(json!({
        "zeilen": index.rows(), "vergleiche_ueber_stichprobe": deep_ties.load(AtomicOrdering::Relaxed),
    }));
    Ok(index)
}

/// Baut den Index blockweise mit festem RAM-Bedarf.
pub fn build(
    reference: Reference,
    variant: Variant,
    options: Options,
) -> Result<FmIndex, BuildError> {
    let threads = thread_count(options.threads);
    gc_log::info(
        "INDEX.BAU.PARAMETER",
        "Index-Variante",
        json!({
            "variante": variant.name(), "verfahren": "blockweise", "kerne": threads,
            "suffixe_je_block": options.block_suffixes, "referenz_basen": reference.len(),
            "contigs": reference.contigs.len(), "n_strecken": reference.holes.len(),
        }),
    );
    build_with(
        reference,
        variant,
        Params {
            v: 1024,
            k: 12,
            block_suffixes: options.block_suffixes.max(1),
            threads,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fm::Rng;

    #[test]
    fn covers_are_complete() {
        for v in [4u64, 5, 7, 8, 13, 16, 31, 32, 64, 100, 256, 1024, 2048] {
            let c = Cover::new(v);
            assert!(c.shift.iter().all(|&s| s != u32::MAX), "v = {v}");
            for i in 0..3 * v {
                for j in [i + 1, i + 2, i + v - 1, i + v, i + 3 * v + 5] {
                    let d = c.common_offset(i, j);
                    assert!(d < v);
                    assert_ne!(c.index_of[((i + d) % v) as usize], u32::MAX);
                    assert_ne!(c.index_of[((j + d) % v) as usize], u32::MAX);
                }
            }
            if v == 1024 {
                assert_eq!(c.size(), 40);
            }
        }
    }

    fn genome(seed: u64, len: usize, style: u64) -> Vec<u8> {
        let mut rng = Rng(seed | 1);
        let mut s: Vec<u8> = (0..len)
            .map(|_| b"ACGT"[(rng.next() & 3) as usize])
            .collect();
        match style % 5 {
            0 => {}
            1 => {
                // lange exakte Wiederholungen, auch als reverses Komplement
                let unit = s[..len / 5].to_vec();
                let rc: Vec<u8> = unit
                    .iter()
                    .rev()
                    .map(|b| match b {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        _ => b'A',
                    })
                    .collect();
                s[len / 4..len / 4 + unit.len()].copy_from_slice(&unit);
                s[len / 2..len / 2 + rc.len()].copy_from_slice(&rc);
            }
            2 => {
                // langes Tandem mit kurzer Periode
                let period = 1 + (rng.next() % 7) as usize;
                for i in len / 3..len - len / 6 {
                    s[i] = s[len / 3 + (i - len / 3) % period];
                }
            }
            3 => s.fill(b'A'),
            _ => {
                for (i, b) in s.iter_mut().enumerate() {
                    *b = b"AC"[i % 2];
                }
                if len > 10 {
                    s[len / 2] = b'N';
                }
            }
        }
        s
    }

    #[test]
    fn default_parameters_on_small_genome() {
        let seq = genome(5, 5000, 1);
        let mut a = Reference::new();
        a.push_contig("x", &seq);
        let mut b = Reference::new();
        b.push_contig("x", &seq);
        let (want, _) = FmIndex::build(a, Variant::Kompakt, 0).expect("direkter Bau");
        let got = build(
            b,
            Variant::Kompakt,
            Options {
                threads: 2,
                block_suffixes: 3000,
            },
        )
        .expect("blockweise");
        assert!(got == want);
    }

    #[test]
    fn equals_direct_build() {
        let mut rng = Rng(77);
        for round in 0..60u64 {
            let mut reference = Reference::new();
            let mut direct = Reference::new();
            let contigs = 1 + round % 3;
            for c in 0..contigs {
                let len = 1 + (rng.next() % 900) as usize;
                let seq = genome(rng.next(), len, round + c);
                reference.push_contig(&format!("c{c}"), &seq);
                direct.push_contig(&format!("c{c}"), &seq);
            }
            let variant = if round % 2 == 0 {
                Variant::Kompakt
            } else {
                Variant::Schnell
            };
            let (want, _) = FmIndex::build(direct, variant, 200).expect("direkter Bau");
            let params = Params {
                v: [4, 8, 16, 64, 1024][(round % 5) as usize],
                k: [1, 2, 3, 6, 9][(round / 5 % 5) as usize],
                block_suffixes: [1, 7, 100, 1_000_000][(round / 3 % 4) as usize],
                threads: 1 + (round % 4) as usize,
            };
            let got = build_with(reference, variant, params).expect("blockweiser Bau");
            assert!(got == want, "Runde {round}: {params:?}");
        }
    }
}
