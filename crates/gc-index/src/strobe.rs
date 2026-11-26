//! Randstrobe-Index für den Strobemer-Weg von `fq2bam`.
//!
//! Seeds wie bei strobealign (Sahlin 2022): Aus dem Text werden offene Syncmere gewählt
//! (k = 20, s = 16, das kleinste s-mer steht an Stelle 2), und jedes Syncmer wird mit einem
//! zweiten aus den 2. bis 9. folgenden zu einem Randstrobe verbunden (das mit der kleinsten
//! Bitdifferenz der Hashes). Ein Randstrobe überspannt so 25 bis etwa 70 Basen und bleibt trotz
//! einzelner Abweichungen zwischen den beiden Strobes erhalten.
//!
//! Der Index enthält nur den Vorwärtsstrang. Ein Read wird mit beiden Orientierungen
//! nachgeschlagen. Je Randstrobe ein Eintrag von 8 Bytes: 32 Bit des Hashes und die Position im
//! Vorwärtsstrang. Die obersten 24 Bit des Hashes wählen das Fach, die nächsten 32 stehen im
//! Eintrag; zusammen unterscheiden 56 Bit. Für GRCh38 sind das rund 587 Mio. Einträge (4,7 GB)
//! und 64 MB Fächer. N-Strecken (`holes`) liefern keine Syncmere.

use crate::reference::Reference;
use std::sync::atomic::{AtomicUsize, Ordering};

pub const K: usize = 20;
pub const S: usize = 16;
const T: usize = 3;
pub const W_MIN: usize = 2;
pub const W_MAX: usize = 9;
const BUCKET_BITS: u32 = 24;

#[inline]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^= x >> 33;
    x
}

/// Offene Syncmere einer Folge von Codes (0 bis 3; größer = unbekannt): (Anfang, Hash des k-mers).
pub fn syncmers(seq: &[u8], out: &mut Vec<(u32, u64)>, scratch: &mut Vec<u64>) {
    out.clear();
    if seq.len() < K {
        return;
    }
    let smask = (1u64 << (2 * S)) - 1;
    let kmask = (1u64 << (2 * K)) - 1;
    let (mut sv, mut kv) = (0u64, 0u64);
    let mut valid = 0usize;
    scratch.clear();
    scratch.resize(seq.len(), 0);
    for (i, &c) in seq.iter().enumerate() {
        if c > 3 {
            valid = 0;
            sv = 0;
            kv = 0;
            continue;
        }
        sv = ((sv << 2) | u64::from(c)) & smask;
        kv = ((kv << 2) | u64::from(c)) & kmask;
        valid += 1;
        if valid >= S {
            scratch[i + 1 - S] = mix(sv);
        }
        if valid >= K {
            let start = i + 1 - K;
            let w = &scratch[start..start + K - S + 1];
            let mut m = 0;
            for j in 1..w.len() {
                if w[j] < w[m] {
                    m = j;
                }
            }
            if m == T - 1 {
                out.push((start as u32, mix(kv ^ 0x9e37_79b9_7f4a_7c15)));
            }
        }
    }
}

/// Randstrobes aus Syncmeren: (Hash, Anfang des ersten Strobes).
pub fn randstrobes(sy: &[(u32, u64)], out: &mut Vec<(u64, u32)>) {
    out.clear();
    for i in 0..sy.len() {
        let lo = i + W_MIN;
        if lo >= sy.len() {
            break;
        }
        let hi = (i + W_MAX).min(sy.len() - 1);
        let h1 = sy[i].1;
        let (mut best, mut bv) = (lo, u32::MAX);
        for (j, s) in sy.iter().enumerate().take(hi + 1).skip(lo) {
            let v = (h1 ^ s.1).count_ones();
            if v < bv {
                bv = v;
                best = j;
            }
        }
        out.push(((h1 >> 1).wrapping_add(sy[best].1 / 3), sy[i].0));
    }
}

/// Der Index.
pub struct StrobeIndex {
    /// Anfang jedes Fachs in `entries`, ein Eintrag mehr als Fächer.
    buckets: Vec<u32>,
    /// Je Fach nach Schlüssel sortiert: Schlüssel (32 Bit) << 32 | Position.
    entries: Vec<u64>,
}

impl std::fmt::Debug for StrobeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "StrobeIndex({} Einträge)", self.entries.len())
    }
}

/// Stücke von 4 Mb mit 10 kb Überhang: Ein Randstrobe braucht nur die folgenden 9 Syncmere.
const CHUNK: u64 = 4 << 20;
const OVER: u64 = 10_000;

impl StrobeIndex {
    /// Baut den Index über den Vorwärtsstrang der Referenz, auf `threads` Kernen.
    pub fn build(reference: &Reference, threads: usize) -> Self {
        let total = reference.len();
        assert!(
            total < u64::from(u32::MAX),
            "Referenz zu lang für 32-Bit-Positionen"
        );
        // Stücke je Contig: (Anfang, Ende der Anfänge, Ende der Folge).
        let mut pieces = Vec::new();
        for c in &reference.contigs {
            let e = c.offset + c.len;
            let mut a = c.offset;
            while a < e {
                let b = (a + CHUNK).min(e);
                pieces.push((a, b, (b + OVER).min(e)));
                a = b;
            }
        }
        let strobes_of = |&(a, b, e): &(u64, u64, u64), out: &mut Vec<(u64, u32)>| {
            let mut seq = Vec::with_capacity((e - a) as usize);
            reference.forward.extend_range(a, e, &mut seq);
            // N-Strecken unbekannt machen
            let first = reference.holes.partition_point(|h| h.offset + h.len <= a);
            for h in reference.holes[first..].iter().take_while(|h| h.offset < e) {
                for p in h.offset.max(a)..(h.offset + h.len).min(e) {
                    seq[(p - a) as usize] = 4;
                }
            }
            let (mut sy, mut scratch) = (Vec::new(), Vec::new());
            syncmers(&seq, &mut sy, &mut scratch);
            randstrobes(&sy, out);
            out.retain(|x| u64::from(x.1) < b - a);
            for x in out.iter_mut() {
                x.1 += a as u32;
            }
        };
        let threads = threads.max(1);
        let nb = 1usize << BUCKET_BITS;
        // 1. Zählen
        let next = AtomicUsize::new(0);
        let counts: Vec<Vec<u32>> = std::thread::scope(|s| {
            let hs: Vec<_> = (0..threads)
                .map(|_| {
                    s.spawn(|| {
                        let mut c = vec![0u32; nb];
                        let mut st = Vec::new();
                        loop {
                            let i = next.fetch_add(1, Ordering::Relaxed);
                            let Some(p) = pieces.get(i) else { break };
                            strobes_of(p, &mut st);
                            for &(h, _) in &st {
                                c[(h >> (64 - BUCKET_BITS)) as usize] += 1;
                            }
                        }
                        c
                    })
                })
                .collect();
            hs.into_iter()
                .map(|h| h.join().expect("Zähl-Thread"))
                .collect()
        });
        let mut buckets = vec![0u32; nb + 1];
        for c in &counts {
            for (i, &v) in c.iter().enumerate() {
                buckets[i + 1] += v;
            }
        }
        drop(counts);
        for i in 1..buckets.len() {
            buckets[i] += buckets[i - 1];
        }
        // 2. Füllen, Stück für Stück in der Reihenfolge der Stücke (gleiches Ergebnis bei jeder
        // Zahl von Kernen); gerechnet werden die Strobes parallel in Gruppen.
        let mut entries = vec![0u64; buckets[nb] as usize];
        let mut fill: Vec<u32> = buckets[..nb].to_vec();
        for group in pieces.chunks(threads * 4) {
            let lists: Vec<Vec<(u64, u32)>> = std::thread::scope(|s| {
                let hs: Vec<_> = group
                    .chunks(group.len().div_ceil(threads))
                    .map(|part| {
                        s.spawn(move || {
                            part.iter()
                                .map(|p| {
                                    let mut st = Vec::new();
                                    strobes_of(p, &mut st);
                                    st
                                })
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                hs.into_iter()
                    .flat_map(|h| h.join().expect("Strobe-Thread"))
                    .collect()
            });
            for list in lists {
                for (h, pos) in list {
                    let bi = (h >> (64 - BUCKET_BITS)) as usize;
                    entries[fill[bi] as usize] = (((h >> 8) & 0xffff_ffff) << 32) | u64::from(pos);
                    fill[bi] += 1;
                }
            }
        }
        drop(fill);
        // 3. Je Fach sortieren
        {
            let mut slices: Vec<&mut [u64]> = Vec::with_capacity(nb);
            let mut rest: &mut [u64] = &mut entries;
            for i in 0..nb {
                let (a, b) = rest.split_at_mut((buckets[i + 1] - buckets[i]) as usize);
                slices.push(a);
                rest = b;
            }
            let per = slices.len().div_ceil(threads);
            std::thread::scope(|s| {
                for part in slices.chunks_mut(per) {
                    s.spawn(move || {
                        for sl in part.iter_mut() {
                            sl.sort_unstable();
                        }
                    });
                }
            });
        }
        StrobeIndex { buckets, entries }
    }

    /// Zahl der Einträge.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Bytes im Speicher.
    pub fn bytes(&self) -> usize {
        self.entries.len() * 8 + self.buckets.len() * 4
    }

    /// Positionen (Vorwärtsstrang, Anfang des ersten Strobes) eines Randstrobe-Hashes.
    /// `None`, wenn es mehr als `max` gibt.
    pub fn lookup(&self, h: u64, max: usize) -> Option<impl Iterator<Item = u64> + '_> {
        let bi = (h >> (64 - BUCKET_BITS)) as usize;
        let e = &self.entries[self.buckets[bi] as usize..self.buckets[bi + 1] as usize];
        let key = (h >> 8) & 0xffff_ffff;
        let a = e.partition_point(|&x| (x >> 32) < key);
        let n = e[a..]
            .iter()
            .take(max + 1)
            .take_while(|&&x| (x >> 32) == key)
            .count();
        (n <= max).then(|| e[a..a + n].iter().map(|&x| x & 0xffff_ffff))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_strobe_of_the_reference_is_found_at_its_place() {
        let mut state = 0x1234_5678_9abc_def1u64;
        let mut seq = Vec::new();
        for i in 0..30_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            seq.push(if (12_000..12_050).contains(&i) {
                b'N'
            } else {
                b"ACGT"[(state & 3) as usize]
            });
        }
        let mut reference = Reference::new();
        reference.push_contig("a", &seq[..20_000]);
        reference.push_contig("b", &seq[20_000..]);
        let idx = StrobeIndex::build(&reference, 3);
        // Dieselben Strobes, Contig für Contig ohne Stücke gerechnet
        let mut want = 0;
        for c in &reference.contigs {
            let mut codes = Vec::new();
            reference
                .forward
                .extend_range(c.offset, c.offset + c.len, &mut codes);
            for h in &reference.holes {
                for p in h.offset.max(c.offset)..(h.offset + h.len).min(c.offset + c.len) {
                    codes[(p - c.offset) as usize] = 4;
                }
            }
            let (mut sy, mut sc, mut st) = (Vec::new(), Vec::new(), Vec::new());
            syncmers(&codes, &mut sy, &mut sc);
            randstrobes(&sy, &mut st);
            for &(h, p) in &st {
                let at: Vec<u64> = idx.lookup(h, 1000).unwrap().collect();
                assert!(at.contains(&(c.offset + u64::from(p))));
            }
            want += st.len();
        }
        assert_eq!(idx.len(), want);
        assert!(want > 3000);
        // mit einem und mit vielen Kernen gleich
        let one = StrobeIndex::build(&reference, 1);
        assert_eq!(one.entries, idx.entries);
        assert_eq!(one.buckets, idx.buckets);
    }
}
