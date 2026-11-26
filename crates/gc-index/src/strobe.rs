//! Randstrobe-Index für den Strobemer-Weg von `fq2bam`.
//!
//! Seeds wie bei strobealign (Sahlin 2022): Aus dem Text werden offene Syncmere gewählt
//! (k = 20, s = 16, das kleinste s-mer steht an Stelle 2), und jedes Syncmer wird mit einem
//! zweiten aus den 2. bis 9. folgenden zu einem Randstrobe verbunden (das mit der kleinsten
//! Bitdifferenz der Hashes). Ein Randstrobe überspannt so 25 bis etwa 70 Basen und bleibt trotz
//! einzelner Abweichungen zwischen den beiden Strobes erhalten.
//!
//! Der Index enthält nur den Vorwärtsstrang. Ein Read wird mit beiden Orientierungen
//! nachgeschlagen. Je Ort 6 Bytes: 16 Bit des Hashes und die Position im Vorwärtsstrang; die
//! obersten 24 Bit des Hashes wählen das Fach, zusammen unterscheiden 40 Bit. Randstrobes mit mehr
//! als [`KEEP_MAX`] Orten stehen nur in einer Sperrliste. Für GRCh38 rund 577 Mio. Orte (3,5 GB)
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
    /// Anfang jedes Fachs in `keys` und `pos`, ein Eintrag mehr als Fächer.
    buckets: Vec<u32>,
    /// Je Fach aufsteigend: 16 weitere Bit des Hashes, dazu die Position.
    keys: Vec<u16>,
    pos: Vec<u32>,
    /// Hashes (Fach << 16 | Schlüssel) mit mehr als [`KEEP_MAX`] Orten; ihre Orte sind nicht
    /// gespeichert, sie gelten als gesperrt.
    blocked: Vec<u64>,
}

/// Höchstens so viele Orte je Randstrobe werden gespeichert.
pub const KEEP_MAX: usize = 1000;

impl std::fmt::Debug for StrobeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "StrobeIndex({} Einträge, {} gesperrt)",
            self.keys.len(),
            self.blocked.len()
        )
    }
}

/// Stücke von 4 Mb mit 10 kb Überhang: Ein Randstrobe braucht nur die folgenden 9 Syncmere.
const CHUNK: u64 = 4 << 20;
const OVER: u64 = 10_000;
/// Der Bau füllt die Fächer in so vielen Durchgängen; so braucht er neben dem fertigen Index
/// nur ein Viertel der Einträge zusätzlich.
const PASSES: usize = 4;

#[inline]
fn bucket_of(h: u64) -> usize {
    (h >> (64 - BUCKET_BITS)) as usize
}

#[inline]
fn key_of(h: u64) -> u16 {
    (h >> (64 - BUCKET_BITS - 16)) as u16
}

impl StrobeIndex {
    /// Baut den Index über den Vorwärtsstrang der Referenz, auf `threads` Kernen. Das Ergebnis
    /// hängt nicht von der Zahl der Kerne ab.
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
        // Alle Stücke auf allen Kernen, je Stück `f` mit seinen Strobes; die Ergebnisse je Kern.
        let over_pieces = |f: &(dyn Fn(&[(u64, u32)], &mut Vec<(u64, u32)>) + Sync)| {
            let next = AtomicUsize::new(0);
            std::thread::scope(|s| {
                let hs: Vec<_> = (0..threads)
                    .map(|_| {
                        s.spawn(|| {
                            let (mut st, mut kept) = (Vec::new(), Vec::new());
                            let mut mine = Vec::new();
                            loop {
                                let i = next.fetch_add(1, Ordering::Relaxed);
                                let Some(p) = pieces.get(i) else { break };
                                strobes_of(p, &mut st);
                                kept.clear();
                                f(&st, &mut kept);
                                mine.extend_from_slice(&kept);
                            }
                            mine
                        })
                    })
                    .collect();
                hs.into_iter()
                    .map(|h| h.join().expect("Strobe-Thread"))
                    .collect::<Vec<Vec<(u64, u32)>>>()
            })
        };
        // 1. Zählen je Fach
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
                                c[bucket_of(h)] += 1;
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
        let mut count = vec![0u64; nb];
        for c in &counts {
            for (i, &v) in c.iter().enumerate() {
                count[i] += u64::from(v);
            }
        }
        drop(counts);
        let all: u64 = count.iter().sum();
        let mut keys: Vec<u16> = Vec::with_capacity(all as usize);
        let mut pos: Vec<u32> = Vec::with_capacity(all as usize);
        let mut blocked = Vec::new();
        let mut buckets = vec![0u32; nb + 1];
        // 2. In Durchgängen über Bereiche von Fächern: sammeln, je Fach sortieren, verdichten.
        let per_pass = nb.div_ceil(PASSES);
        for pass in 0..PASSES {
            let (b0, b1) = (pass * per_pass, ((pass + 1) * per_pass).min(nb));
            let mut start = vec![0usize; b1 - b0 + 1];
            for i in b0..b1 {
                start[i - b0 + 1] = start[i - b0] + count[i] as usize;
            }
            let mut tmp = vec![0u64; start[b1 - b0]];
            let mut fill: Vec<usize> = start[..b1 - b0].to_vec();
            let lists = over_pieces(&|st, kept| {
                kept.extend(st.iter().filter(|x| (b0..b1).contains(&bucket_of(x.0))));
            });
            for list in lists {
                for (h, p) in list {
                    let bi = bucket_of(h) - b0;
                    tmp[fill[bi]] = (u64::from(key_of(h)) << 32) | u64::from(p);
                    fill[bi] += 1;
                }
            }
            drop(fill);
            {
                let mut slices: Vec<&mut [u64]> = Vec::with_capacity(b1 - b0);
                let mut rest: &mut [u64] = &mut tmp;
                for i in 0..b1 - b0 {
                    let (a, b) = rest.split_at_mut(start[i + 1] - start[i]);
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
            for i in 0..b1 - b0 {
                let sl = &tmp[start[i]..start[i + 1]];
                let mut a = 0;
                while a < sl.len() {
                    let k = sl[a] >> 32;
                    let mut e = a + 1;
                    while e < sl.len() && sl[e] >> 32 == k {
                        e += 1;
                    }
                    if e - a > KEEP_MAX {
                        blocked.push(((b0 + i) as u64) << 16 | k);
                    } else {
                        for &x in &sl[a..e] {
                            keys.push(k as u16);
                            pos.push(x as u32);
                        }
                    }
                    a = e;
                }
                buckets[b0 + i + 1] = keys.len() as u32;
            }
        }
        StrobeIndex {
            buckets,
            keys,
            pos,
            blocked,
        }
    }

    /// Zahl der gespeicherten Orte.
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Zahl der gesperrten Randstrobes (mehr als [`KEEP_MAX`] Orte).
    pub fn blocked(&self) -> usize {
        self.blocked.len()
    }

    /// Bytes im Speicher.
    pub fn bytes(&self) -> usize {
        self.keys.len() * 6 + self.buckets.len() * 4 + self.blocked.len() * 8
    }

    /// Positionen (Vorwärtsstrang, Anfang des ersten Strobes) eines Randstrobe-Hashes.
    /// `None`, wenn es mehr als `max` gibt (oder mehr als [`KEEP_MAX`]).
    pub fn lookup(&self, h: u64, max: usize) -> Option<impl Iterator<Item = u64> + '_> {
        let bi = bucket_of(h);
        let k = key_of(h);
        if self
            .blocked
            .binary_search(&((bi as u64) << 16 | u64::from(k)))
            .is_ok()
        {
            return None;
        }
        let (lo, hi) = (self.buckets[bi] as usize, self.buckets[bi + 1] as usize);
        let ks = &self.keys[lo..hi];
        let a = ks.partition_point(|&x| x < k);
        let n = ks[a..]
            .iter()
            .take(max + 1)
            .take_while(|&&x| x == k)
            .count();
        (n <= max).then(|| self.pos[lo + a..lo + a + n].iter().map(|&p| u64::from(p)))
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
        assert_eq!(one.keys, idx.keys);
        assert_eq!(one.pos, idx.pos);
        assert_eq!(one.buckets, idx.buckets);
        assert_eq!(one.blocked, idx.blocked);
    }

    #[test]
    fn strobes_of_a_repeat_with_many_copies_are_blocked() {
        let mut state = 0x0bad_cafe_u64;
        let unit: Vec<u8> = (0..300)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect();
        let seq: Vec<u8> = unit.iter().cycle().take(300 * 1200).copied().collect();
        let mut reference = Reference::new();
        reference.push_contig("wdh", &seq);
        let idx = StrobeIndex::build(&reference, 2);
        assert!(idx.blocked() > 0);
        // drei Einheiten hintereinander, geprüft werden die Strobes der mittleren
        let codes: Vec<u8> = unit
            .iter()
            .cycle()
            .take(900)
            .map(|&b| crate::reference::encode(b).unwrap())
            .collect();
        let (mut sy, mut sc, mut st) = (Vec::new(), Vec::new(), Vec::new());
        syncmers(&codes, &mut sy, &mut sc);
        randstrobes(&sy, &mut st);
        st.retain(|x| (300..600).contains(&x.1));
        assert!(!st.is_empty());
        for &(h, _) in &st {
            assert!(idx.lookup(h, usize::MAX).is_none());
        }
    }
}
