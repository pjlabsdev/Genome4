//! Suffix-Array in linearer Zeit nach SA-IS (Nong, Zhang, Chan 2009).
//!
//! Voraussetzung: das letzte Zeichen des Textes ist ein Endzeichen, das genau einmal vorkommt
//! und kleiner ist als alle anderen. Positionen sind `u32`, der Text darf also höchstens
//! `u32::MAX - 1` Zeichen haben.

const EMPTY: u32 = u32::MAX;

/// Größte Textlänge, die diese Umsetzung verarbeitet.
pub const MAX_LEN: usize = (u32::MAX - 1) as usize;

struct Bits(Vec<u64>);

impl Bits {
    fn new(n: usize) -> Self {
        Bits(vec![0; n.div_ceil(64)])
    }
    #[inline]
    fn get(&self, i: usize) -> bool {
        self.0[i >> 6] >> (i & 63) & 1 == 1
    }
    #[inline]
    fn set(&mut self, i: usize) {
        self.0[i >> 6] |= 1 << (i & 63);
    }
}

fn bucket_sizes<T: Copy + Into<u32>>(s: &[T], k: usize) -> Vec<u32> {
    let mut sizes = vec![0u32; k];
    for &c in s {
        sizes[c.into() as usize] += 1;
    }
    sizes
}

fn bucket_heads(sizes: &[u32], heads: &mut [u32]) {
    let mut sum = 0;
    for (h, &size) in heads.iter_mut().zip(sizes) {
        *h = sum;
        sum += size;
    }
}

fn bucket_tails(sizes: &[u32], tails: &mut [u32]) {
    let mut sum = 0;
    for (t, &size) in tails.iter_mut().zip(sizes) {
        sum += size;
        *t = sum;
    }
}

/// Aus den einsortierten LMS-Suffixen erst alle L-, dann alle S-Suffixe ableiten.
fn induce<T: Copy + Into<u32>>(
    s: &[T],
    sa: &mut [u32],
    is_s: &Bits,
    sizes: &[u32],
    bkt: &mut [u32],
) {
    let n = s.len();
    bucket_heads(sizes, bkt);
    for i in 0..n {
        let j = sa[i];
        if j != EMPTY && j > 0 {
            let p = (j - 1) as usize;
            if !is_s.get(p) {
                let c = s[p].into() as usize;
                sa[bkt[c] as usize] = p as u32;
                bkt[c] += 1;
            }
        }
    }
    bucket_tails(sizes, bkt);
    for i in (0..n).rev() {
        let j = sa[i];
        if j != EMPTY && j > 0 {
            let p = (j - 1) as usize;
            if is_s.get(p) {
                let c = s[p].into() as usize;
                bkt[c] -= 1;
                sa[bkt[c] as usize] = p as u32;
            }
        }
    }
}

/// Suffix-Array von `s` über das Alphabet `0..k`.
pub fn suffix_array<T: Copy + Into<u32>>(s: &[T], k: usize) -> Vec<u32> {
    let n = s.len();
    assert!(n <= MAX_LEN, "Text zu lang für 32-Bit-Positionen");
    let mut sa = vec![EMPTY; n];
    if n == 0 {
        return sa;
    }
    if n == 1 {
        sa[0] = 0;
        return sa;
    }
    let ch = |i: usize| -> u32 { s[i].into() };

    // Typen: S, wenn das Suffix kleiner ist als das folgende.
    let mut is_s = Bits::new(n);
    is_s.set(n - 1);
    for i in (0..n - 1).rev() {
        let (a, b) = (ch(i), ch(i + 1));
        if a < b || (a == b && is_s.get(i + 1)) {
            is_s.set(i);
        }
    }
    let is_lms = |i: usize| i > 0 && is_s.get(i) && !is_s.get(i - 1);

    let sizes = bucket_sizes(s, k);
    let mut bkt = vec![0u32; k];

    // Schritt 1: LMS-Suffixe ans Ende ihrer Fächer, dann ableiten. Danach sind die
    // LMS-Teilstücke sortiert.
    bucket_tails(&sizes, &mut bkt);
    let mut n_lms = 0usize;
    for i in (1..n).rev() {
        if is_lms(i) {
            let c = ch(i) as usize;
            bkt[c] -= 1;
            sa[bkt[c] as usize] = i as u32;
            n_lms += 1;
        }
    }
    induce(s, &mut sa, &is_s, &sizes, &mut bkt);

    // Schritt 2: LMS-Teilstücke benennen (gleiche Stücke, gleicher Name).
    let mut names = vec![EMPTY; n / 2 + 1];
    let mut name = 0u32;
    let mut prev = usize::MAX;
    for &entry in sa.iter() {
        let j = entry as usize;
        if !is_lms(j) {
            continue;
        }
        if prev != usize::MAX {
            let mut d = 0;
            let differ = loop {
                let (a, b) = (j + d, prev + d);
                if ch(a) != ch(b) || is_s.get(a) != is_s.get(b) {
                    break true;
                }
                if d > 0 && (is_lms(a) || is_lms(b)) {
                    break !(is_lms(a) && is_lms(b));
                }
                d += 1;
            };
            if differ {
                name += 1;
            }
        }
        names[j / 2] = name;
        prev = j;
    }
    let n_names = name as usize + 1;

    // Verkürzter Text aus den Namen in Textreihenfolge, dazu die LMS-Positionen.
    let mut lms_pos = Vec::with_capacity(n_lms);
    let mut s1 = Vec::with_capacity(n_lms);
    for i in 1..n {
        if is_lms(i) {
            lms_pos.push(i as u32);
            s1.push(names[i / 2]);
        }
    }
    drop(names);

    // Schritt 3: Reihenfolge der LMS-Suffixe, direkt oder über Rekursion.
    let sa1 = if n_names < n_lms {
        suffix_array(&s1, n_names)
    } else {
        let mut direct = vec![0u32; n_lms];
        for (i, &nm) in s1.iter().enumerate() {
            direct[nm as usize] = i as u32;
        }
        direct
    };
    drop(s1);

    // Schritt 4: LMS-Suffixe in richtiger Reihenfolge einsortieren und endgültig ableiten.
    sa.fill(EMPTY);
    bucket_tails(&sizes, &mut bkt);
    for &r in sa1.iter().rev() {
        let j = lms_pos[r as usize];
        let c = ch(j as usize) as usize;
        bkt[c] -= 1;
        sa[bkt[c] as usize] = j;
    }
    drop(sa1);
    drop(lms_pos);
    induce(s, &mut sa, &is_s, &sizes, &mut bkt);
    sa
}

/// Prüft in linearer Zeit, ob `sa` das Suffix-Array von `s` ist (gleiche Voraussetzung wie oben).
/// Liefert die erste fehlerhafte Zeile.
pub fn verify<T: Copy + Into<u32>>(s: &[T], sa: &[u32]) -> Result<(), usize> {
    let n = s.len();
    if sa.len() != n {
        return Err(0);
    }
    let mut rank = vec![EMPTY; n + 1];
    for (i, &p) in sa.iter().enumerate() {
        if p as usize >= n || rank[p as usize] != EMPTY {
            return Err(i);
        }
        rank[p as usize] = i as u32 + 1;
    }
    rank[n] = 0; // das leere Suffix ist das kleinste
    for i in 1..n {
        let (a, b) = (sa[i - 1] as usize, sa[i] as usize);
        let (ca, cb): (u32, u32) = (s[a].into(), s[b].into());
        if ca > cb || (ca == cb && rank[a + 1] >= rank[b + 1]) {
            return Err(i);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive(s: &[u8]) -> Vec<u32> {
        let mut sa: Vec<u32> = (0..s.len() as u32).collect();
        sa.sort_by(|&a, &b| s[a as usize..].cmp(&s[b as usize..]));
        sa
    }

    fn with_sentinel(body: &[u8]) -> Vec<u8> {
        let mut s: Vec<u8> = body.iter().map(|c| c + 1).collect();
        s.push(0);
        s
    }

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
    }

    #[test]
    fn small_known_cases() {
        for body in [
            &b""[..],
            b"a",
            b"aa",
            b"ab",
            b"ba",
            b"banana",
            b"mississippi",
            b"aaaaaaaa",
            b"abababab",
            b"abcabcabcabc",
        ] {
            let s = with_sentinel(body);
            let sa = suffix_array(&s, 256);
            assert_eq!(sa, naive(&s), "{:?}", std::str::from_utf8(body));
            assert_eq!(verify(&s, &sa), Ok(()));
        }
    }

    #[test]
    fn random_texts_match_naive_sort() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        for round in 0..400 {
            let len = (rng.next() % 300) as usize;
            let sigma = [1, 2, 4, 4, 20][round % 5];
            // Wiederholungen erzwingen, damit die Rekursion mehrfach greift.
            let unit: Vec<u8> = (0..1 + rng.next() % 7)
                .map(|_| (rng.next() % sigma) as u8)
                .collect();
            let body: Vec<u8> = (0..len)
                .map(|i| {
                    if rng.next().is_multiple_of(4) {
                        (rng.next() % sigma) as u8
                    } else {
                        unit[i % unit.len()]
                    }
                })
                .collect();
            let s = with_sentinel(&body);
            let sa = suffix_array(&s, sigma as usize + 1);
            assert_eq!(sa, naive(&s), "Runde {round}");
            assert_eq!(verify(&s, &sa), Ok(()));
        }
    }

    #[test]
    fn verify_detects_errors() {
        let s = with_sentinel(b"mississippi");
        let mut sa = suffix_array(&s, 256);
        sa.swap(3, 4);
        assert!(verify(&s, &sa).is_err());
        sa.swap(3, 4);
        sa[5] = sa[6];
        assert!(verify(&s, &sa).is_err());
    }

    #[test]
    fn larger_dna_like_text() {
        let mut rng = Rng(42);
        let mut body: Vec<u8> = (0..200_000).map(|_| (rng.next() % 4) as u8).collect();
        // lange exakte Wiederholung und Homopolymer
        let copy = body[1000..31_000].to_vec();
        body.extend_from_slice(&copy);
        body.extend(std::iter::repeat_n(0u8, 5000));
        let s = with_sentinel(&body);
        let sa = suffix_array(&s, 5);
        assert_eq!(verify(&s, &sa), Ok(()));
    }
}
