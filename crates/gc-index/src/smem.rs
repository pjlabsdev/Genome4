//! Seed-Suche: Suche in beide Richtungen auf dem FM-Index und längste exakte Treffer (SMEM).
//!
//! Der Indextext ist der Vorwärtsstrang gefolgt von seinem reversen Komplement. Er ist damit
//! sein eigenes reverses Komplement, und jedes Vorkommen eines Musters `P` an Textposition `i`
//! hat ein Vorkommen von `rk(P)` an Position `2L - i - |P|`. Deshalb reicht eine einzige BWT, um
//! ein Muster nach links und nach rechts zu verlängern (FMD-Index, Li 2012): Man führt die
//! Zeilenbereiche von `P` und von `rk(P)` gemeinsam.

use crate::fm::FmIndex;

/// Zeilenbereiche eines Musters `P` und seines reversen Komplements, beide gleich groß.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BiInterval {
    /// Erste Zeile der Suffixe, die mit `P` beginnen.
    pub k: u64,
    /// Erste Zeile der Suffixe, die mit `rk(P)` beginnen.
    pub kr: u64,
    /// Anzahl der Zeilen, also Vorkommen von `P` im Indextext (beide Stränge).
    pub size: u64,
}

/// Ein längster exakter Treffer eines Reads: `read[start..end]` kommt `interval.size`-mal vor
/// und lässt sich in keine Richtung verlängern, ohne unter die verlangte Häufigkeit zu fallen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Smem {
    pub start: u32,
    pub end: u32,
    pub interval: BiInterval,
}

impl Smem {
    pub fn len(&self) -> u32 {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }
}

impl FmIndex {
    /// Bereich eines einzelnen Zeichens (Code 0..3).
    pub fn bi_start(&self, c: u8) -> BiInterval {
        let end = |c: u8| {
            if c == 3 {
                self.rows()
            } else {
                self.first_row(c + 1)
            }
        };
        BiInterval {
            k: self.first_row(c),
            kr: self.first_row(3 - c),
            size: end(c) - self.first_row(c),
        }
    }

    /// Verlängert das Muster nach links um jede der vier Basen: Ergebnis `[cP für c in ACGT]`.
    pub fn extend_left_all(&self, iv: BiInterval) -> [BiInterval; 4] {
        let mut out = [BiInterval {
            k: 0,
            kr: 0,
            size: 0,
        }; 4];
        let mut total = 0;
        for c in 0..4u8 {
            let lo = self.occ(c, iv.k);
            let hi = self.occ(c, iv.k + iv.size);
            out[c as usize].k = self.first_row(c) + lo;
            out[c as usize].size = hi - lo;
            total += hi - lo;
        }
        // Im Bereich von rk(P) stehen die Fortsetzungen in der Reihenfolge Endzeichen, A, C, G, T.
        // rk(P)·x kommt so oft vor wie komp(x)·P. Das Endzeichen folgt höchstens einmal.
        let mut kr = iv.kr + (iv.size - total);
        for c in (0..4usize).rev() {
            out[c].kr = kr;
            kr += out[c].size;
        }
        out
    }

    /// Verlängert das Muster nach links um die Base `c`.
    pub fn extend_left(&self, iv: BiInterval, c: u8) -> BiInterval {
        self.extend_left_all(iv)[c as usize]
    }

    /// Verlängert das Muster nach rechts um die Base `c`: Linksverlängerung von `rk(P)` um
    /// das Komplement.
    pub fn extend_right(&self, iv: BiInterval, c: u8) -> BiInterval {
        let swapped = BiInterval {
            k: iv.kr,
            kr: iv.k,
            size: iv.size,
        };
        let r = self.extend_left(swapped, 3 - c);
        BiInterval {
            k: r.kr,
            kr: r.k,
            size: r.size,
        }
    }

    /// Alle längsten Treffer, die die Read-Position `x` überdecken und mindestens `min_size`-mal
    /// vorkommen. Hängt sie an `out` an, längster Anfang zuletzt, und gibt die Position zurück,
    /// an der der am weitesten nach rechts reichende Treffer endet. `read` enthält Codes 0..3,
    /// jeder andere Wert gilt als N und beendet jeden Treffer.
    ///
    /// Ablauf wie in BWA-MEM (Li 2013): erst von `x` nach rechts verlängern und jeden Bereich
    /// merken, bevor er kleiner wird; dann alle gemerkten Bereiche gemeinsam nach links
    /// verlängern, die weitest reichenden zuerst.
    pub fn smems_at(&self, read: &[u8], x: usize, min_size: u64, out: &mut Vec<Smem>) -> usize {
        let min_size = min_size.max(1);
        if x >= read.len() || read[x] > 3 {
            return x + 1;
        }
        // Nach rechts.
        let mut forward: Vec<(usize, BiInterval)> = Vec::new();
        let mut iv = self.bi_start(read[x]);
        let mut end = x + 1;
        while end < read.len() && read[end] <= 3 {
            let next = self.extend_right(iv, read[end]);
            if next.size != iv.size {
                forward.push((end, iv));
                if next.size < min_size {
                    break;
                }
            }
            iv = next;
            end += 1;
        }
        if forward.last().is_none_or(|&(e, _)| e != end) && iv.size >= min_size {
            forward.push((end, iv));
        }
        if forward.is_empty() {
            return x + 1;
        }
        let reach = forward.last().map_or(x + 1, |&(e, _)| e);
        forward.reverse();

        // Nach links, alle gemeinsam; `forward` ist nach fallendem Ende geordnet.
        let first_new = out.len();
        let mut alive = forward;
        let mut next_alive: Vec<(usize, BiInterval)> = Vec::with_capacity(alive.len());
        let mut i = x as isize - 1;
        loop {
            let c = if i >= 0 { read[i as usize] } else { 4 };
            next_alive.clear();
            for &(e, iv) in &alive {
                let ext = if c <= 3 {
                    self.extend_left(iv, c)
                } else {
                    BiInterval {
                        k: 0,
                        kr: 0,
                        size: 0,
                    }
                };
                if ext.size < min_size {
                    // Nicht weiter verlängerbar. Ein Treffer ist es nur, wenn kein weiter nach
                    // rechts reichender Bereich diesen Schritt überlebt hat und nicht schon einer
                    // mit demselben Anfang gemeldet wurde.
                    let start = (i + 1) as usize;
                    let contained = out[first_new..]
                        .last()
                        .is_some_and(|m| m.start as usize <= start);
                    if next_alive.is_empty() && !contained {
                        out.push(Smem {
                            start: start as u32,
                            end: e as u32,
                            interval: iv,
                        });
                    }
                } else if next_alive.last().is_none_or(|&(_, l)| l.size != ext.size) {
                    next_alive.push((e, ext));
                }
            }
            if next_alive.is_empty() {
                break;
            }
            std::mem::swap(&mut alive, &mut next_alive);
            i -= 1;
        }
        out[first_new..].reverse();
        reach
    }

    /// Alle längsten exakten Treffer (SMEM) eines Reads mit mindestens `min_len` Basen,
    /// nach Anfang geordnet.
    pub fn smems(&self, read: &[u8], min_len: u32) -> Vec<Smem> {
        let mut out = Vec::new();
        let mut x = 0;
        while x < read.len() {
            x = self.smems_at(read, x, 1, &mut out);
        }
        out.retain(|m| m.len() >= min_len);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fm::{Rng, Variant};
    use crate::reference::Reference;

    /// Indextext ohne Endzeichen als Codes 0..3.
    fn text_of(index: &FmIndex) -> Vec<u8> {
        let t = index.reference().index_text();
        t[..t.len() - 1].iter().map(|c| c - 1).collect()
    }

    fn occurrences(text: &[u8], p: &[u8]) -> u64 {
        if p.is_empty() || p.len() > text.len() {
            return 0;
        }
        text.windows(p.len()).filter(|w| *w == p).count() as u64
    }

    fn genome(seed: u64, len: usize) -> FmIndex {
        let mut rng = Rng(seed);
        let mut s: Vec<u8> = (0..len)
            .map(|_| b"ACGT"[(rng.next() & 3) as usize])
            .collect();
        // Wiederholungen in beiden Richtungen, damit Bereiche mit mehreren Zeilen vorkommen.
        let copy = s[10..90].to_vec();
        s[150..230].copy_from_slice(&copy);
        let rc: Vec<u8> = s[300..360]
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        s[400..460].copy_from_slice(&rc);
        for b in &mut s[500..520] {
            *b = b'N';
        }
        let mut r = Reference::new();
        r.push_contig("a", &s[..len / 2]);
        r.push_contig("b", &s[len / 2..]);
        FmIndex::build(r, Variant::Kompakt, 0).unwrap().0
    }

    /// Ein Read: Stück des Textes mit Fehlern und gelegentlich einem N.
    fn read_from(text: &[u8], rng: &mut Rng) -> Vec<u8> {
        let len = 20 + (rng.next() % 80) as usize;
        let start = (rng.next() as usize) % (text.len() - len);
        let mut r = text[start..start + len].to_vec();
        for _ in 0..rng.next() % 5 {
            let i = (rng.next() as usize) % len;
            r[i] = (rng.next() & 3) as u8;
        }
        if rng.next().is_multiple_of(4) {
            let i = (rng.next() as usize) % len;
            r[i] = 4;
        }
        r
    }

    #[test]
    fn extensions_match_plain_search() {
        let index = genome(3, 1200);
        let text = text_of(&index);
        let mut rng = Rng(11);
        for _ in 0..300 {
            let len = 1 + (rng.next() % 30) as usize;
            let start = (rng.next() as usize) % (text.len() - len);
            let mut p = text[start..start + len].to_vec();
            if rng.next().is_multiple_of(3) {
                let i = (rng.next() as usize) % len;
                p[i] = (rng.next() & 3) as u8;
            }
            let rc: Vec<u8> = p.iter().rev().map(|c| 3 - c).collect();
            let check = |iv: BiInterval, upto: &[u8]| {
                let rcu: Vec<u8> = upto.iter().rev().map(|c| 3 - c).collect();
                let (lo, hi) = index.interval(upto);
                let (rlo, rhi) = index.interval(&rcu);
                assert_eq!(iv.size, hi - lo);
                assert_eq!(iv.size, rhi - rlo);
                if iv.size > 0 {
                    assert_eq!((iv.k, iv.kr), (lo, rlo), "Muster {upto:?}");
                }
            };
            // nach rechts aufbauen
            let mut iv = index.bi_start(p[0]);
            check(iv, &p[..1]);
            for i in 1..len {
                iv = index.extend_right(iv, p[i]);
                check(iv, &p[..=i]);
                if iv.size == 0 {
                    break;
                }
            }
            // nach links aufbauen
            let mut iv = index.bi_start(p[len - 1]);
            for i in (0..len - 1).rev() {
                iv = index.extend_left(iv, p[i]);
                check(iv, &p[i..]);
                if iv.size == 0 {
                    break;
                }
            }
            let _ = rc;
        }
    }

    /// Unabhängige Wahrheit: alle Abschnitte, die x überdecken, mindestens `min_size`-mal
    /// vorkommen und in keine Richtung verlängerbar sind.
    fn naive_at(text: &[u8], read: &[u8], x: usize, min_size: u64) -> Vec<(usize, usize, u64)> {
        let ok = |a: usize, b: usize| -> u64 {
            if read[a..b].iter().any(|&c| c > 3) {
                0
            } else {
                occurrences(text, &read[a..b])
            }
        };
        let mut out = Vec::new();
        for a in 0..=x {
            for b in x + 1..=read.len() {
                let n = ok(a, b);
                if n < min_size {
                    continue;
                }
                let left = a > 0 && ok(a - 1, b) >= min_size;
                let right = b < read.len() && ok(a, b + 1) >= min_size;
                if !left && !right {
                    out.push((a, b, n));
                }
            }
        }
        out
    }

    #[test]
    fn smems_at_match_naive_search() {
        let index = genome(5, 1200);
        let text = text_of(&index);
        let mut rng = Rng(21);
        let mut found = 0;
        for round in 0..120 {
            let read = read_from(&text, &mut rng);
            let x = (rng.next() as usize) % read.len();
            let min_size = 1 + round % 3;
            let mut got = Vec::new();
            let reach = index.smems_at(&read, x, min_size, &mut got);
            let got: Vec<(usize, usize, u64)> = got
                .iter()
                .map(|m| (m.start as usize, m.end as usize, m.interval.size))
                .collect();
            let want = naive_at(&text, &read, x, min_size);
            assert_eq!(
                got, want,
                "Runde {round} x={x} min={min_size} Read {read:?}"
            );
            let far = want.iter().map(|w| w.1).max().unwrap_or(x + 1);
            assert_eq!(reach, far);
            found += want.len();
        }
        assert!(found > 150);
    }

    #[test]
    fn smems_of_read_match_naive_search() {
        let index = genome(8, 1200);
        let text = text_of(&index);
        let mut rng = Rng(31);
        for round in 0..150 {
            let read = read_from(&text, &mut rng);
            // Längste Übereinstimmung ab jeder Position.
            let longest: Vec<usize> = (0..read.len())
                .map(|i| {
                    let mut l = 0;
                    while i + l < read.len()
                        && read[i + l] <= 3
                        && occurrences(&text, &read[i..=i + l]) > 0
                    {
                        l += 1;
                    }
                    l
                })
                .collect();
            let mut want = Vec::new();
            let mut prev_end = 0;
            for (i, &l) in longest.iter().enumerate() {
                if l > 0 && i + l > prev_end {
                    want.push((i, i + l));
                }
                prev_end = prev_end.max(i + l);
            }
            let got = index.smems(&read, 1);
            let pairs: Vec<(usize, usize)> = got
                .iter()
                .map(|m| (m.start as usize, m.end as usize))
                .collect();
            assert_eq!(pairs, want, "Runde {round} Read {read:?}");
            for m in &got {
                let p = &read[m.start as usize..m.end as usize];
                assert_eq!(
                    index.interval(p),
                    (m.interval.k, m.interval.k + m.interval.size)
                );
            }
            // Mindestlänge filtert nur.
            let long: Vec<Smem> = got.iter().copied().filter(|m| m.len() >= 12).collect();
            assert_eq!(index.smems(&read, 12), long);
        }
    }
}
