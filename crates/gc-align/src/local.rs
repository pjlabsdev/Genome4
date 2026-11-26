//! Lokales Alignment (Smith-Waterman) eines Reads gegen ein Stück Referenz, mit bestem und
//! zweitbestem Score und den Enden (BWA: `ksw_align2` mit `ksw_u8` und `ksw_i16`).
//!
//! bwa rechnet das mit Vektorbefehlen in "gestreifter" Anordnung (Farrar 2007): Der Read wird
//! auf 16 oder 8 Spuren verteilt, die gleichzeitig laufen. Diese Fassung rechnet dieselben
//! Schritte in derselben Anordnung und liefert auch in den Feinheiten dasselbe: Eine Deletion
//! direkt nach einer Insertion gibt es nicht, und was die Spuren erst in der Nachkorrektur
//! erfahren, fließt nicht mehr in die offene Deletion ein.
//!
//! Die Spuren sind feste Felder (16 Bytes oder 8 Werte mit 16 Bit) und die Schritte gesättigte
//! Additionen, Subtraktionen und Maxima. Daraus macht der Übersetzer auf jedem System
//! Vektorbefehle; systemnaher Code ist nicht nötig, und das Ergebnis ist überall dasselbe.

use crate::extend::Scoring;
use gc_sys::simd::{I16x8, Lanes, U8x16};

/// Ergebnis. Unbekannte Werte sind -1.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocalResult {
    pub score: i32,
    /// Letzte ausgerichtete Position in der Referenz und im Read (einschließlich).
    pub te: i32,
    pub qe: i32,
    /// Bester Score, der an einer anderen Stelle der Referenz endet, und wo.
    pub score2: i32,
    pub te2: i32,
    /// Erste ausgerichtete Position in Referenz und Read; nur mit `want_start`.
    pub tb: i32,
    pub qb: i32,
}

const NONE: LocalResult = LocalResult {
    score: 0,
    te: -1,
    qe: -1,
    score2: -1,
    te2: -1,
    tb: -1,
    qb: -1,
};

/// Ein Wert einer Spur: ein Byte ohne Vorzeichen (16 Spuren) oder 16 Bit (8 Spuren). Die
/// Rechenschritte sind so geschrieben, dass der Übersetzer aus den Schleifen über die Spuren
/// Vektorbefehle macht (SSE2, NEON), ohne dass hier systemnaher Code steht.
trait Cell: Copy + Ord + Default + Send + Sync {
    /// Ob die Punkte verschoben gespeichert sind (Byte-Lauf).
    const BYTE: bool;
    fn from_i32(v: i32) -> Self;
    fn to_i32(self) -> i32;
}

impl Cell for u8 {
    const BYTE: bool = true;
    #[inline(always)]
    fn from_i32(v: i32) -> Self {
        v as u8
    }
    #[inline(always)]
    fn to_i32(self) -> i32 {
        i32::from(self)
    }
}

impl Cell for i16 {
    const BYTE: bool = false;
    #[inline(always)]
    fn from_i32(v: i32) -> Self {
        v as i16
    }
    #[inline(always)]
    fn to_i32(self) -> i32 {
        i32::from(self)
    }
}

/// Ein Durchlauf mit `L` Spuren: 16 (Scores in einem Byte, höchstens 255) oder 8 (16 Bit).
/// `min_sub`: ab diesem Score merkt sich der Lauf Zeilen für den zweitbesten Score.
/// `stop_at`: Abbruch, sobald dieser Score erreicht ist.
// Die Schleifen über die Spuren laufen absichtlich über den Index: So erkennt der Übersetzer
// die feste Länge und setzt Vektorbefehle ein.
#[allow(clippy::needless_range_loop)]
fn striped<T: Cell, V: Lanes<Elem = T>, const L: usize>(
    query: &[u8],
    target: &[u8],
    sc: &Scoring,
    min_sub: i32,
    stop_at: i32,
) -> LocalResult {
    let qlen = query.len();
    let slen = qlen.div_ceil(L);
    // Im Byte-Lauf sind alle Punkte um den größten Abzug verschoben, damit nichts negativ wird.
    let shift = if T::BYTE {
        sc.mismatch.max(sc.ambiguous)
    } else {
        0
    };
    let best_pair = sc.matching.max(-sc.ambiguous);
    let vz = V::zero();
    // Ein einziger Speicherblock für Punkte und Arbeitszeilen (vorher fünf je Aufruf).
    let mut block = vec![vz; 9 * slen];
    let (profile, rest) = block.split_at_mut(5 * slen);
    let (h0, rest) = rest.split_at_mut(slen);
    let (h1, rest) = rest.split_at_mut(slen);
    let (e_row, h_max) = rest.split_at_mut(slen);
    let (mut h0, mut h1) = (h0, h1);
    // Punkte je Referenzbase, in gestreifter Anordnung: Vektor j, Spur k gehört zu Read-Position
    // j + k * slen. Hinter dem Read-Ende stehen 0 Punkte.
    for a in 0..5u8 {
        for j in 0..slen {
            let mut lanes = [T::default(); L];
            for (k, lane) in lanes.iter_mut().enumerate() {
                let p = j + k * slen;
                let s = if p < qlen { sc.score(a, query[p]) } else { 0 };
                *lane = T::from_i32(s + shift);
            }
            profile[usize::from(a) * slen + j] = V::load(&lanes);
        }
    }
    let profile = &*profile;
    let shift_v = V::splat(T::from_i32(shift));
    let (e_del, e_ins) = (
        V::splat(T::from_i32(sc.e_del)),
        V::splat(T::from_i32(sc.e_ins)),
    );
    let oe_del = V::splat(T::from_i32(sc.o_del + sc.e_del));
    let oe_ins = V::splat(T::from_i32(sc.o_ins + sc.e_ins));
    let add = |h: V, s: V| {
        if T::BYTE {
            h.adds(s).subs(shift_v)
        } else {
            h.adds(s)
        }
    };
    let mut marks: Vec<(i32, i32)> = Vec::new(); // (Score, Zeile)
    let (mut gmax, mut te) = (0i32, -1i32);
    for (i, &base) in target.iter().enumerate() {
        let s = &profile[usize::from(base.min(4)) * slen..][..slen];
        // h = H der vorigen Zeile, letzte Read-Position jeder Spur, um eine Spur verschoben.
        let mut h = h0[slen - 1].shift_in_zero();
        let mut f = vz;
        let mut row_max = vz;
        for j in 0..slen {
            // `t` hängt nicht von `f` ab. Weil das begrenzte Abziehen mit dem Maximum
            // vertauscht und `oe_ins >= e_ins` ist, gilt
            // max(f - e_ins, max(t, f) - oe_ins) = max(f - e_ins, t - oe_ins): Die Kette von
            // Spalte zu Spalte ist so nur zwei Schritte lang statt drei, bei gleichem Ergebnis.
            let t = add(h, s[j]).max(e_row[j]);
            let v = t.max(f);
            row_max = row_max.max(v);
            h1[j] = v;
            e_row[j] = e_row[j].subs(e_del).max(v.subs(oe_del));
            f = f.subs(e_ins).max(t.subs(oe_ins));
            h = h0[j];
        }
        // Nachkorrektur: die offene Insertion von Spur zu Spur weiterreichen.
        'lazy: for _ in 0..16 {
            f = f.shift_in_zero();
            for h1j in h1.iter_mut() {
                let v = h1j.max(f);
                *h1j = v;
                f = f.subs(e_ins);
                if !f.any_gt(v.subs(oe_ins)) {
                    break 'lazy;
                }
            }
        }
        let imax = row_max.hmax().to_i32();
        let i = i as i32;
        if imax >= min_sub {
            match marks.last_mut() {
                Some(last) if last.1 + 1 == i => {
                    if last.0 < imax {
                        *last = (imax, i);
                    }
                }
                _ => marks.push((imax, i)),
            }
        }
        if imax > gmax {
            gmax = imax;
            te = i;
            h_max.copy_from_slice(h1);
            if (T::BYTE && gmax + shift >= 255) || gmax >= stop_at {
                break;
            }
        }
        std::mem::swap(&mut h0, &mut h1);
    }
    let mut r = NONE;
    r.score = if T::BYTE && gmax + shift >= 255 {
        255
    } else {
        gmax
    };
    r.te = te;
    if !(T::BYTE && r.score == 255) {
        let mut max = -1;
        let mut lanes = [T::default(); 16];
        for (j, hv) in h_max.iter().enumerate() {
            hv.store(&mut lanes);
            for (k, &v) in lanes[..L].iter().enumerate() {
                let (v, p) = (v.to_i32(), (j + k * slen) as i32);
                if v > max {
                    max = v;
                    r.qe = p;
                } else if v == max && p < r.qe {
                    r.qe = p;
                }
            }
        }
        if !marks.is_empty() {
            let d = (r.score + best_pair - 1) / best_pair;
            let (low, high) = (te - d, te + d);
            for &(score, row) in &marks {
                if (row < low || row > high) && score > r.score2 {
                    r.score2 = score;
                    r.te2 = row;
                }
            }
        }
    }
    r
}

/// Ein Durchlauf mit Scores in einem Byte (16 Spuren) oder in 16 Bit (8 Spuren).
fn run(
    query: &[u8],
    target: &[u8],
    sc: &Scoring,
    byte: bool,
    min_sub: i32,
    stop_at: i32,
) -> LocalResult {
    use gc_sys::simd::plain;
    match (byte, vectors_on()) {
        (true, true) => striped::<u8, U8x16, 16>(query, target, sc, min_sub, stop_at),
        (false, true) => striped::<i16, I16x8, 8>(query, target, sc, min_sub, stop_at),
        (true, false) => striped::<u8, plain::U8x16, 16>(query, target, sc, min_sub, stop_at),
        (false, false) => striped::<i16, plain::I16x8, 8>(query, target, sc, min_sub, stop_at),
    }
}

/// Schalter `GENOME_SIMD=0`: Vektoren als gewöhnliche Felder rechnen (Stand vor den
/// Vektorbefehlen). Das Ergebnis ist in beiden Fällen dasselbe.
pub(crate) fn vectors_on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GENOME_SIMD").map_or(true, |v| v != "0"))
}

/// Dieselbe Rechnung Wert für Wert mit 32 Bit, ohne jede Rücksicht auf Tempo. Sie war die erste
/// Fassung und ist gegen bwa geprüft; die schnelle Fassung muss in jedem Wert dasselbe liefern.
#[cfg(test)]
fn striped_reference(
    query: &[u8],
    target: &[u8],
    sc: &Scoring,
    lanes: usize,
    min_sub: i32,
    stop_at: i32,
) -> LocalResult {
    let byte = lanes == 16;
    let qlen = query.len();
    let slen = qlen.div_ceil(lanes);
    // Im Byte-Lauf sind alle Punkte um den größten Abzug verschoben, damit nichts negativ wird.
    let shift = if byte {
        sc.mismatch.max(sc.ambiguous)
    } else {
        0
    };
    let cap = if byte { 255 } else { i32::from(i16::MAX) };
    let best_pair = sc.matching.max(-sc.ambiguous);
    // Punkte je Referenzbase, in gestreifter Anordnung: Vektor j, Spur k gehört zu Read-Position
    // j + k * slen. Hinter dem Read-Ende stehen 0 Punkte.
    let mut profile = vec![0i32; 5 * slen * lanes];
    for a in 0..5u8 {
        for j in 0..slen {
            for k in 0..lanes {
                let p = j + k * slen;
                let s = if p < qlen { sc.score(a, query[p]) } else { 0 };
                profile[(usize::from(a) * slen + j) * lanes + k] = s + shift;
            }
        }
    }
    let n = slen * lanes;
    let (mut h0, mut h1) = (vec![0i32; n], vec![0i32; n]);
    let (mut e_row, mut h_max) = (vec![0i32; n], vec![0i32; n]);
    let (oe_del, oe_ins) = (sc.o_del + sc.e_del, sc.o_ins + sc.e_ins);
    let sub = |a: i32, b: i32| (a - b).max(0);
    let mut marks: Vec<(i32, i32)> = Vec::new(); // (Score, Zeile)
    let (mut gmax, mut te) = (0i32, -1i32);
    let mut h = vec![0i32; lanes];
    let mut f = vec![0i32; lanes];
    let mut row_max = vec![0i32; lanes];
    for (i, &base) in target.iter().enumerate() {
        let s = &profile[usize::from(base.min(4)) * n..][..n];
        // h = H der vorigen Zeile, letzte Read-Position jeder Spur, um eine Spur verschoben.
        h[0] = 0;
        for k in 1..lanes {
            h[k] = h0[(slen - 1) * lanes + k - 1];
        }
        f.fill(0);
        row_max.fill(0);
        for j in 0..slen {
            for k in 0..lanes {
                let at = j * lanes + k;
                let mut v = if byte {
                    sub((h[k] + s[at]).min(cap), shift)
                } else {
                    (h[k] + s[at]).min(cap)
                };
                let mut e = e_row[at];
                v = v.max(e).max(f[k]);
                row_max[k] = row_max[k].max(v);
                h1[at] = v;
                e = sub(e, sc.e_del).max(sub(v, oe_del));
                e_row[at] = e;
                f[k] = sub(f[k], sc.e_ins).max(sub(v, oe_ins));
                h[k] = h0[at];
            }
        }
        // Nachkorrektur: die offene Insertion von Spur zu Spur weiterreichen.
        'lazy: for _ in 0..16 {
            for k in (1..lanes).rev() {
                f[k] = f[k - 1];
            }
            f[0] = 0;
            for j in 0..slen {
                let mut settled = true;
                for (k, fk) in f.iter_mut().enumerate() {
                    let at = j * lanes + k;
                    let v = h1[at].max(*fk);
                    h1[at] = v;
                    *fk = sub(*fk, sc.e_ins);
                    settled &= *fk <= sub(v, oe_ins);
                }
                if settled {
                    break 'lazy;
                }
            }
        }
        let imax = row_max.iter().copied().max().unwrap_or(0);
        let i = i as i32;
        if imax >= min_sub {
            match marks.last_mut() {
                Some(last) if last.1 + 1 == i => {
                    if last.0 < imax {
                        *last = (imax, i);
                    }
                }
                _ => marks.push((imax, i)),
            }
        }
        if imax > gmax {
            gmax = imax;
            te = i;
            h_max.copy_from_slice(&h1);
            if (byte && gmax + shift >= 255) || gmax >= stop_at {
                break;
            }
        }
        std::mem::swap(&mut h0, &mut h1);
    }
    let mut r = NONE;
    r.score = if byte && gmax + shift >= 255 {
        255
    } else {
        gmax
    };
    r.te = te;
    if !(byte && r.score == 255) {
        let mut max = -1;
        for (idx, &v) in h_max.iter().enumerate() {
            let p = (idx / lanes + idx % lanes * slen) as i32;
            if v > max {
                max = v;
                r.qe = p;
            } else if v == max && p < r.qe {
                r.qe = p;
            }
        }
        if !marks.is_empty() {
            let d = (r.score + best_pair - 1) / best_pair;
            let (low, high) = (te - d, te + d);
            for &(score, row) in &marks {
                if (row < low || row > high) && score > r.score2 {
                    r.score2 = score;
                    r.te2 = row;
                }
            }
        }
    }
    r
}

/// Lokales Alignment von `query` gegen `target`. `byte`: mit Scores in einem Byte rechnen
/// (bwa tut das, wenn der höchste mögliche Score unter 250 liegt). `min_sub`: Schwelle für den
/// zweitbesten Score; liegt der beste darunter, wird der Anfang nicht gesucht. Mit
/// `want_start` wird das Alignment vom Ende aus rückwärts wiederholt, um den Anfang zu finden.
pub fn local_align(
    query: &[u8],
    target: &[u8],
    sc: &Scoring,
    byte: bool,
    min_sub: i32,
    want_start: bool,
) -> LocalResult {
    if query.is_empty() {
        return NONE;
    }
    let mut r = run(query, target, sc, byte, min_sub, 0x10000);
    if !want_start || r.score < min_sub || r.qe < 0 || r.te < 0 {
        return r;
    }
    // Rückwärts: Read bis qe und Referenz bis te umdrehen. Wie in bwa läuft der zweite Durchlauf
    // über die ganze Länge der Referenz; er bricht ab, sobald der Score erreicht ist.
    let q: Vec<u8> = query[..=r.qe as usize].iter().rev().copied().collect();
    let mut t = target.to_vec();
    t[..=r.te as usize].reverse();
    let rr = run(&q, &t, sc, byte, 0x10000, r.score);
    if rr.score == r.score {
        r.tb = r.te - rr.te;
        r.qb = r.qe - rr.qe;
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn random(state: &mut u64, n: usize) -> Vec<u8> {
        (0..n).map(|_| (xorshift(state) % 4) as u8).collect()
    }

    /// Lehrbuch-Rechnung nach Gotoh, lokal. Obere Schranke für den Score.
    fn textbook(q: &[u8], t: &[u8], sc: &Scoring) -> i32 {
        let (n, m) = (t.len(), q.len());
        let mut h = vec![vec![0i32; m + 1]; n + 1];
        let mut e = vec![vec![0i32; m + 1]; n + 1];
        let mut f = vec![vec![0i32; m + 1]; n + 1];
        let mut best = 0;
        for i in 1..=n {
            for j in 1..=m {
                e[i][j] = (e[i - 1][j] - sc.e_del).max(h[i - 1][j] - sc.o_del - sc.e_del);
                f[i][j] = (f[i][j - 1] - sc.e_ins).max(h[i][j - 1] - sc.o_ins - sc.e_ins);
                h[i][j] = (h[i - 1][j - 1] + sc.score(t[i - 1], q[j - 1]))
                    .max(e[i][j])
                    .max(f[i][j])
                    .max(0);
                best = best.max(h[i][j]);
            }
        }
        best
    }

    #[test]
    fn known_cases_give_score_and_both_ends() {
        let sc = Scoring::default();
        let mut state = 0x0123_4567_89AB_CDEFu64;
        let t = random(&mut state, 400);
        for byte in [true, false] {
            // Read genau aus der Referenz
            let q = t[150..250].to_vec();
            let r = local_align(&q, &t, &sc, byte, 19, true);
            assert_eq!((r.score, r.tb, r.te, r.qb, r.qe), (100, 150, 249, 0, 99));
            // ein Austausch in der Mitte
            let mut q1 = q.clone();
            q1[50] = (q1[50] + 1) & 3;
            let r = local_align(&q1, &t, &sc, byte, 19, true);
            assert_eq!((r.score, r.tb, r.te, r.qb, r.qe), (95, 150, 249, 0, 99));
            // 3 Basen der Referenz fehlen im Read
            let q2: Vec<u8> = [&t[150..200], &t[203..253]].concat();
            let r = local_align(&q2, &t, &sc, byte, 19, true);
            assert_eq!(
                (r.score, r.tb, r.te, r.qb, r.qe),
                (100 - 9, 150, 252, 0, 99)
            );
            // zehn fremde Basen vorn: das Alignment beginnt erst dahinter
            let mut q3 = q.clone();
            for b in q3.iter_mut().take(10) {
                *b = (*b + 2) & 3;
            }
            let r = local_align(&q3, &t, &sc, byte, 19, true);
            assert_eq!((r.score, r.tb, r.te, r.qb, r.qe), (90, 160, 249, 10, 99));
        }
    }

    #[test]
    fn second_best_is_reported_when_the_read_fits_twice() {
        let sc = Scoring::default();
        let mut state = 0x0F0F_1234_5678_9ABCu64;
        let mut t = random(&mut state, 600);
        let copy = t[100..180].to_vec();
        t[400..480].copy_from_slice(&copy);
        t[440] = (t[440] + 1) & 3; // zweite Kopie mit einem Austausch
        let r = local_align(&copy, &t, &sc, true, 19, true);
        assert_eq!((r.score, r.te, r.tb), (80, 179, 100));
        assert_eq!((r.score2, r.te2), (75, 479));
        // ohne zweite Stelle bleibt der zweitbeste Score unter der Schwelle unbekannt
        let r = local_align(&t[200..300], &t[150..350], &sc, true, 60, false);
        assert_eq!((r.score, r.score2), (100, -1));
    }

    /// Die schnelle Fassung gegen die Rechnung Wert für Wert: jeder Wert des Ergebnisses, beide
    /// Breiten, mit und ohne vorzeitigen Abbruch, auch mit N und mit Scores über 255.
    #[test]
    fn fast_version_equals_the_reference_in_every_value() {
        let sc = Scoring::default();
        let mut state = 0x5151_DEAD_BEEF_0001u64;
        for round in 0..3000 {
            let tlen = 1 + (xorshift(&mut state) % 700) as usize;
            let mut t = random(&mut state, tlen);
            let qlen = 1 + (xorshift(&mut state) % 320) as usize;
            let start = (xorshift(&mut state) as usize) % tlen;
            let mut q: Vec<u8> = t[start..(start + qlen).min(tlen)].to_vec();
            if round % 7 == 0 {
                q = random(&mut state, qlen);
            }
            for _ in 0..(xorshift(&mut state) % 10) {
                let p = (xorshift(&mut state) as usize) % q.len();
                match xorshift(&mut state) % 4 {
                    0 => q[p] = (xorshift(&mut state) % 5) as u8,
                    1 if q.len() > 1 => {
                        q.remove(p);
                    }
                    2 => q.insert(p, (xorshift(&mut state) % 4) as u8),
                    _ => {
                        let at = (xorshift(&mut state) as usize) % tlen;
                        t[at] = 4;
                    }
                }
            }
            let min_sub = [0, 19, 40, 0x10000][(xorshift(&mut state) % 4) as usize];
            let stop_at = [0x10000, 30, 77][(xorshift(&mut state) % 3) as usize];
            for byte in [true, false] {
                let lanes = if byte { 16 } else { 8 };
                assert_eq!(
                    run(&q, &t, &sc, byte, min_sub, stop_at),
                    striped_reference(&q, &t, &sc, lanes, min_sub, stop_at),
                    "Runde {round}, Byte {byte}"
                );
            }
        }
    }

    #[test]
    fn both_widths_agree_and_stay_between_ungapped_and_textbook() {
        let sc = Scoring::default();
        let mut state = 0x7777_AAAA_1111_3333u64;
        for round in 0..300 {
            let tlen = 60 + (xorshift(&mut state) % 200) as usize;
            let t = random(&mut state, tlen);
            let qlen = 20 + (xorshift(&mut state) % 120) as usize;
            let start = (xorshift(&mut state) as usize) % tlen;
            let mut q: Vec<u8> = t[start..(start + qlen).min(tlen)].to_vec();
            for _ in 0..(xorshift(&mut state) % 6) {
                if q.len() < 5 {
                    break;
                }
                let p = (xorshift(&mut state) as usize) % q.len();
                match xorshift(&mut state) % 3 {
                    0 => q[p] = (xorshift(&mut state) % 5) as u8,
                    1 => {
                        q.remove(p);
                    }
                    _ => q.insert(p, (xorshift(&mut state) % 4) as u8),
                }
            }
            let a = local_align(&q, &t, &sc, true, 19, true);
            let b = local_align(&q, &t, &sc, false, 19, true);
            assert_eq!(a, b, "Runde {round}: Byte und 16 Bit");
            assert!(a.score <= textbook(&q, &t, &sc), "Runde {round}");
            if a.score >= 19 {
                assert!(
                    a.tb >= 0 && a.qb >= 0 && a.tb <= a.te && a.qb <= a.qe,
                    "Runde {round}"
                );
            }
        }
    }
}
