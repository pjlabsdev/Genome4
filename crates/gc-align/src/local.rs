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

/// Zustand eines gestreiften Laufs für [`striped_pair`]: dieselben Schritte wie [`striped`],
/// nur in Teile zerlegt, damit zwei Läufe Spalte für Spalte ineinander verzahnt rechnen können.
struct Lane2<'a, T: Cell, V: Lanes<Elem = T>> {
    target: &'a [u8],
    slen: usize,
    shift: i32,
    profile: Vec<V>,
    h0: Vec<V>,
    h1: Vec<V>,
    e_row: Vec<V>,
    h_max: Vec<V>,
    marks: Vec<(i32, i32)>,
    gmax: i32,
    te: i32,
    /// Nächste Zeile der Referenz; `done`, wenn der Lauf fertig ist (Ende oder Abbruch).
    row: usize,
    done: bool,
    qlen: usize,
}

impl<'a, T: Cell, V: Lanes<Elem = T>> Lane2<'a, T, V> {
    #[allow(clippy::needless_range_loop)]
    fn new<const L: usize>(query: &[u8], target: &'a [u8], sc: &Scoring) -> Self {
        let qlen = query.len();
        let slen = qlen.div_ceil(L);
        let shift = if T::BYTE {
            sc.mismatch.max(sc.ambiguous)
        } else {
            0
        };
        let vz = V::zero();
        let mut profile = vec![vz; 5 * slen];
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
        Lane2 {
            target,
            slen,
            shift,
            profile,
            h0: vec![vz; slen],
            h1: vec![vz; slen],
            e_row: vec![vz; slen],
            h_max: vec![vz; slen],
            marks: Vec::new(),
            gmax: 0,
            te: -1,
            row: 0,
            done: target.is_empty(),
            qlen,
        }
    }

    /// Abschluss einer Zeile nach der Hauptschleife: Nachkorrektur, Merkzeilen, Bestwert.
    #[inline(always)]
    fn close_row(&mut self, mut f: V, row_max: V, e_ins: V, oe_ins: V, min_sub: i32, stop_at: i32) {
        'lazy: for _ in 0..16 {
            f = f.shift_in_zero();
            for h1j in self.h1.iter_mut() {
                let v = h1j.max(f);
                *h1j = v;
                f = f.subs(e_ins);
                if !f.any_gt(v.subs(oe_ins)) {
                    break 'lazy;
                }
            }
        }
        let imax = row_max.hmax().to_i32();
        let i = self.row as i32;
        if imax >= min_sub {
            match self.marks.last_mut() {
                Some(last) if last.1 + 1 == i => {
                    if last.0 < imax {
                        *last = (imax, i);
                    }
                }
                _ => self.marks.push((imax, i)),
            }
        }
        self.row += 1;
        if imax > self.gmax {
            self.gmax = imax;
            self.te = i;
            self.h_max.copy_from_slice(&self.h1);
            if (T::BYTE && self.gmax + self.shift >= 255) || self.gmax >= stop_at {
                self.done = true;
                return;
            }
        }
        std::mem::swap(&mut self.h0, &mut self.h1);
        if self.row >= self.target.len() {
            self.done = true;
        }
    }

    fn result<const L: usize>(&self, sc: &Scoring) -> LocalResult {
        let best_pair = sc.matching.max(-sc.ambiguous);
        let mut r = NONE;
        r.score = if T::BYTE && self.gmax + self.shift >= 255 {
            255
        } else {
            self.gmax
        };
        r.te = self.te;
        if !(T::BYTE && r.score == 255) {
            let mut max = -1;
            let mut lanes = [T::default(); 16];
            for (j, hv) in self.h_max.iter().enumerate() {
                hv.store(&mut lanes);
                for (k, &v) in lanes[..L].iter().enumerate() {
                    let (v, p) = (v.to_i32(), (j + k * self.slen) as i32);
                    if v > max {
                        max = v;
                        r.qe = p;
                    } else if v == max && p < r.qe {
                        r.qe = p;
                    }
                }
            }
            if !self.marks.is_empty() {
                let d = (r.score + best_pair - 1) / best_pair;
                let (low, high) = (self.te - d, self.te + d);
                for &(score, row) in &self.marks {
                    if (row < low || row > high) && score > r.score2 {
                        r.score2 = score;
                        r.te2 = row;
                    }
                }
            }
        }
        let _ = self.qlen;
        r
    }
}

/// Zwei Läufe von [`striped`] gleichzeitig, Zeile für Zeile verzahnt. Die Hauptschleife eines
/// Laufs ist eine Kette (jeder Schritt wartet auf den vorigen); zwei unabhängige Ketten
/// nebeneinander halten die Rechenwerke des Prozessors besser beschäftigt. Jeder Lauf rechnet
/// genau dieselben Schritte wie allein, die Ergebnisse sind dieselben wie zweimal [`striped`].
fn striped_pair<T: Cell, V: Lanes<Elem = T>, const L: usize>(
    a: (&[u8], &[u8]),
    b: (&[u8], &[u8]),
    sc: &Scoring,
    min_sub: i32,
    stop_at: i32,
) -> (LocalResult, LocalResult) {
    let mut x = Lane2::<T, V>::new::<L>(a.0, a.1, sc);
    let mut y = Lane2::<T, V>::new::<L>(b.0, b.1, sc);
    let shift_v = V::splat(T::from_i32(x.shift));
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
    let vz = V::zero();
    // Eine Zeile eines Laufs allein.
    let one = |z: &mut Lane2<'_, T, V>| {
        let slen = z.slen;
        let base = z.target[z.row];
        let s = &z.profile[usize::from(base.min(4)) * slen..][..slen];
        let mut h = z.h0[slen - 1].shift_in_zero();
        let (mut f, mut row_max) = (vz, vz);
        for (j, &sj) in s.iter().enumerate() {
            let t = add(h, sj).max(z.e_row[j]);
            let v = t.max(f);
            row_max = row_max.max(v);
            z.h1[j] = v;
            z.e_row[j] = z.e_row[j].subs(e_del).max(v.subs(oe_del));
            f = f.subs(e_ins).max(t.subs(oe_ins));
            h = z.h0[j];
        }
        z.close_row(f, row_max, e_ins, oe_ins, min_sub, stop_at);
    };
    while !x.done && !y.done {
        if x.slen != y.slen {
            one(&mut x);
            one(&mut y);
            continue;
        }
        let slen = x.slen;
        let sx = &x.profile[usize::from(x.target[x.row].min(4)) * slen..][..slen];
        let sy = &y.profile[usize::from(y.target[y.row].min(4)) * slen..][..slen];
        let mut hx = x.h0[slen - 1].shift_in_zero();
        let mut hy = y.h0[slen - 1].shift_in_zero();
        let (mut fx, mut mx, mut fy, mut my) = (vz, vz, vz, vz);
        {
            let (xh0, xh1, xe) = (&x.h0[..slen], &mut x.h1[..slen], &mut x.e_row[..slen]);
            let (yh0, yh1, ye) = (&y.h0[..slen], &mut y.h1[..slen], &mut y.e_row[..slen]);
            let (sx, sy) = (&sx[..slen], &sy[..slen]);
            for j in 0..slen {
                let tx = add(hx, sx[j]).max(xe[j]);
                let ty = add(hy, sy[j]).max(ye[j]);
                let vx = tx.max(fx);
                let vy = ty.max(fy);
                mx = mx.max(vx);
                my = my.max(vy);
                xh1[j] = vx;
                yh1[j] = vy;
                xe[j] = xe[j].subs(e_del).max(vx.subs(oe_del));
                ye[j] = ye[j].subs(e_del).max(vy.subs(oe_del));
                fx = fx.subs(e_ins).max(tx.subs(oe_ins));
                fy = fy.subs(e_ins).max(ty.subs(oe_ins));
                hx = xh0[j];
                hy = yh0[j];
            }
        }
        x.close_row(fx, mx, e_ins, oe_ins, min_sub, stop_at);
        y.close_row(fy, my, e_ins, oe_ins, min_sub, stop_at);
    }
    while !x.done {
        one(&mut x);
    }
    while !y.done {
        one(&mut y);
    }
    (x.result::<L>(sc), y.result::<L>(sc))
}

/// Zwei Vorwärtsläufe von [`run`] auf einmal (siehe [`striped_pair`]). Gleiche Ergebnisse wie
/// zweimal [`run`] mit `stop_at` = 0x10000.
fn run_pair(
    a: (&[u8], &[u8]),
    b: (&[u8], &[u8]),
    sc: &Scoring,
    byte: bool,
    min_sub: i32,
) -> (LocalResult, LocalResult) {
    use gc_sys::simd::plain;
    let stop = 0x10000;
    match (byte, vectors_on()) {
        (true, true) => striped_pair::<u8, U8x16, 16>(a, b, sc, min_sub, stop),
        (false, true) => striped_pair::<i16, I16x8, 8>(a, b, sc, min_sub, stop),
        (true, false) => striped_pair::<u8, plain::U8x16, 16>(a, b, sc, min_sub, stop),
        (false, false) => striped_pair::<i16, plain::I16x8, 8>(a, b, sc, min_sub, stop),
    }
}

/// Zwei Aufrufe von [`local_align`] auf einmal, mit gleichem `byte`, `min_sub` und
/// `want_start`; die Vorwärtsläufe rechnen verzahnt. Gleiche Ergebnisse wie zwei Aufrufe.
pub fn local_align_pair(
    a: (&[u8], &[u8]),
    b: (&[u8], &[u8]),
    sc: &Scoring,
    byte: bool,
    min_sub: i32,
    want_start: bool,
) -> (LocalResult, LocalResult) {
    if a.0.is_empty() || b.0.is_empty() {
        return (
            local_align(a.0, a.1, sc, byte, min_sub, want_start),
            local_align(b.0, b.1, sc, byte, min_sub, want_start),
        );
    }
    let (ra, rb) = run_pair(a, b, sc, byte, min_sub);
    if !want_start {
        return (ra, rb);
    }
    (
        with_start(a.0, a.1, sc, byte, min_sub, ra),
        with_start(b.0, b.1, sc, byte, min_sub, rb),
    )
}

/// Wie [`striped_pair`], mit `N` Läufen nebeneinander. Läufe mit anderer Read-Länge als der
/// erste und Läufe, die übrig bleiben, wenn andere fertig sind, rechnen einzeln weiter.
fn striped_many<T: Cell, V: Lanes<Elem = T>, const L: usize, const N: usize>(
    jobs: [(&[u8], &[u8]); N],
    sc: &Scoring,
    min_sub: i32,
    stop_at: i32,
) -> [LocalResult; N] {
    let mut z: [Lane2<'_, T, V>; N] = jobs.map(|(q, t)| Lane2::<T, V>::new::<L>(q, t, sc));
    let shift_v = V::splat(T::from_i32(z[0].shift));
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
    let vz = V::zero();
    let one = |z: &mut Lane2<'_, T, V>| {
        let slen = z.slen;
        let base = z.target[z.row];
        let s = &z.profile[usize::from(base.min(4)) * slen..][..slen];
        let mut h = z.h0[slen - 1].shift_in_zero();
        let (mut f, mut row_max) = (vz, vz);
        for (j, &sj) in s.iter().enumerate() {
            let t = add(h, sj).max(z.e_row[j]);
            let v = t.max(f);
            row_max = row_max.max(v);
            z.h1[j] = v;
            z.e_row[j] = z.e_row[j].subs(e_del).max(v.subs(oe_del));
            f = f.subs(e_ins).max(t.subs(oe_ins));
            h = z.h0[j];
        }
        z.close_row(f, row_max, e_ins, oe_ins, min_sub, stop_at);
    };
    let slen = z[0].slen;
    let same = z.iter().all(|x| x.slen == slen);
    while same && z.iter().all(|x| !x.done) {
        let mut base = [0usize; N];
        let mut h = [vz; N];
        for n in 0..N {
            base[n] = usize::from(z[n].target[z[n].row].min(4)) * slen;
            h[n] = z[n].h0[slen - 1].shift_in_zero();
        }
        let (mut f, mut m) = ([vz; N], [vz; N]);
        for j in 0..slen {
            for n in 0..N {
                let x = &mut z[n];
                let t = add(h[n], x.profile[base[n] + j]).max(x.e_row[j]);
                let v = t.max(f[n]);
                m[n] = m[n].max(v);
                x.h1[j] = v;
                x.e_row[j] = x.e_row[j].subs(e_del).max(v.subs(oe_del));
                f[n] = f[n].subs(e_ins).max(t.subs(oe_ins));
                h[n] = x.h0[j];
            }
        }
        for n in 0..N {
            z[n].close_row(f[n], m[n], e_ins, oe_ins, min_sub, stop_at);
        }
    }
    for x in z.iter_mut() {
        while !x.done {
            one(x);
        }
    }
    std::array::from_fn(|n| z[n].result::<L>(sc))
}

/// `N` Aufrufe von [`local_align`] auf einmal (siehe [`striped_many`]), mit gleichem `byte`,
/// `min_sub` und `want_start`. Gleiche Ergebnisse wie `N` einzelne Aufrufe.
pub fn local_align_many<const N: usize>(
    jobs: [(&[u8], &[u8]); N],
    sc: &Scoring,
    byte: bool,
    min_sub: i32,
    want_start: bool,
) -> [LocalResult; N] {
    use gc_sys::simd::plain;
    if jobs.iter().any(|j| j.0.is_empty()) {
        return jobs.map(|(q, t)| local_align(q, t, sc, byte, min_sub, want_start));
    }
    let stop = 0x10000;
    let r = match (byte, vectors_on()) {
        (true, true) => striped_many::<u8, U8x16, 16, N>(jobs, sc, min_sub, stop),
        (false, true) => striped_many::<i16, I16x8, 8, N>(jobs, sc, min_sub, stop),
        (true, false) => striped_many::<u8, plain::U8x16, 16, N>(jobs, sc, min_sub, stop),
        (false, false) => striped_many::<i16, plain::I16x8, 8, N>(jobs, sc, min_sub, stop),
    };
    if !want_start {
        return r;
    }
    std::array::from_fn(|n| with_start(jobs[n].0, jobs[n].1, sc, byte, min_sub, r[n]))
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
    let r = run(query, target, sc, byte, min_sub, 0x10000);
    if !want_start {
        return r;
    }
    with_start(query, target, sc, byte, min_sub, r)
}

/// Ergänzt einen Vorwärtslauf um den Anfang des Alignments (zweiter Teil von [`local_align`]).
pub fn with_start(
    query: &[u8],
    target: &[u8],
    sc: &Scoring,
    byte: bool,
    min_sub: i32,
    mut r: LocalResult,
) -> LocalResult {
    if r.score < min_sub || r.qe < 0 || r.te < 0 {
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

/// Vorwärtslauf im Byte-Format, wie ihn eine Grafikeinheit liefert: bester Score `gmax`, seine
/// Zeile `te` und Read-Position `qe`, und je Zeile der Referenz das Zeilenmaximum (vor der
/// Nachkorrektur, wie im gestreiften Lauf). Daraus wird dasselbe [`LocalResult`] wie aus dem
/// Byte-Lauf von [`local_align`] ohne Anfang.
pub fn from_rows(
    gmax: i32,
    te: i32,
    qe: i32,
    row_max: &[u8],
    sc: &Scoring,
    min_sub: i32,
) -> LocalResult {
    let best_pair = sc.matching.max(-sc.ambiguous);
    let mut marks: Vec<(i32, i32)> = Vec::new();
    for (i, &m) in row_max.iter().enumerate() {
        let (imax, i) = (i32::from(m), i as i32);
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
        // Wie im gestreiften Lauf: Nach dem besten Score kommen nur Zeilen bis zu ihm vor,
        // wenn der Lauf dort abbrach; im Byte-Lauf der Partnersuche bricht er nie ab.
    }
    let mut r = NONE;
    r.score = gmax;
    r.te = te;
    r.qe = qe;
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
    r
}

/// Ein Rechenkern für viele Vorwärtsläufe auf einmal (Byte-Lauf, wie [`kernel_reference`]).
/// Liefert je Paar (Read, Referenzstück) `(gmax, te, qe, Zeilenmaxima)` oder `None`, wenn er
/// gerade nicht rechnen kann; dann rechnet der Prozessor.
/// Ergebnis eines Vorwärtslaufs: `(gmax, te, qe, Zeilenmaxima)`.
pub type Forward = (i32, i32, i32, Vec<u8>);

pub trait RescueKernel: Send + Sync + std::fmt::Debug {
    fn forward(&self, jobs: &[(&[u8], &[u8])], sc: &Scoring) -> Option<Vec<Forward>>;
}

/// Ob ein Byte-Lauf für diesen Read sicher nicht sättigt und nie vorzeitig abbricht: Dann gibt
/// [`from_rows`] mit den Werten einer Grafikeinheit dasselbe wie der Lauf auf dem Prozessor.
pub fn rows_exact(query_len: usize, sc: &Scoring) -> bool {
    let shift = sc.mismatch.max(sc.ambiguous);
    query_len as i32 * sc.matching + shift < 250 && query_len <= 256 && query_len > 0
}

/// Der Rechenkern der Grafikeinheit, Schritt für Schritt in Rust: gestreifter Byte-Lauf mit 16
/// Spuren, jede Spur ein eigener Wert. Nur zum Prüfen, dass die Übertragung in die Sprache der
/// Grafikeinheit stimmt (gleiche Schleifen, gleiche Reihenfolge).
// Die Schleifen über die Spuren laufen wie im Rechenkern über den Index.
#[allow(clippy::needless_range_loop)]
pub fn kernel_reference(query: &[u8], target: &[u8], sc: &Scoring) -> Forward {
    const L: usize = 16;
    let qlen = query.len();
    let slen = qlen.div_ceil(L);
    let shift = sc.mismatch.max(sc.ambiguous) as u32;
    let sat_add = |a: u32, b: u32| (a + b).min(255);
    let sat_sub = |a: u32, b: u32| a.saturating_sub(b);
    let score = |a: u8, p: usize| -> u32 {
        if p < qlen {
            (sc.score(a, query[p]) + shift as i32) as u32
        } else {
            shift
        }
    };
    let (e_del, e_ins) = (sc.e_del as u32, sc.e_ins as u32);
    let oe_del = (sc.o_del + sc.e_del) as u32;
    let oe_ins = (sc.o_ins + sc.e_ins) as u32;
    let mut h0 = vec![0u32; slen * L];
    let mut h1 = vec![0u32; slen * L];
    let mut e_row = vec![0u32; slen * L];
    let mut h_max = vec![0u32; slen * L];
    let mut rows = Vec::with_capacity(target.len());
    let (mut gmax, mut te) = (0u32, -1i32);
    for (i, &base) in target.iter().enumerate() {
        let b = base.min(4);
        let mut h = [0u32; L];
        for k in 1..L {
            h[k] = h0[(slen - 1) * L + k - 1];
        }
        let mut f = [0u32; L];
        let mut row_max = [0u32; L];
        for j in 0..slen {
            for k in 0..L {
                let x = j * L + k;
                let t = sat_sub(sat_add(h[k], score(b, j + k * slen)), shift).max(e_row[x]);
                let v = t.max(f[k]);
                row_max[k] = row_max[k].max(v);
                h1[x] = v;
                e_row[x] = sat_sub(e_row[x], e_del).max(sat_sub(v, oe_del));
                f[k] = sat_sub(f[k], e_ins).max(sat_sub(t, oe_ins));
                h[k] = h0[x];
            }
        }
        'lazy: for _ in 0..16 {
            for k in (1..L).rev() {
                f[k] = f[k - 1];
            }
            f[0] = 0;
            for j in 0..slen {
                let mut any = false;
                for k in 0..L {
                    let x = j * L + k;
                    let v = h1[x].max(f[k]);
                    h1[x] = v;
                    f[k] = sat_sub(f[k], e_ins);
                    any |= f[k] > sat_sub(v, oe_ins);
                }
                if !any {
                    break 'lazy;
                }
            }
        }
        let imax = row_max.iter().copied().max().unwrap_or(0);
        rows.push(imax as u8);
        if imax > gmax {
            gmax = imax;
            te = i as i32;
            h_max.copy_from_slice(&h1);
        }
        std::mem::swap(&mut h0, &mut h1);
    }
    let (mut max, mut qe) = (-1i32, -1i32);
    for j in 0..slen {
        for k in 0..L {
            let (v, p) = (h_max[j * L + k] as i32, (j + k * slen) as i32);
            if v > max || (v == max && p < qe) {
                max = v;
                qe = p;
            }
        }
    }
    (gmax as i32, te, qe, rows)
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
    fn the_gpu_kernel_steps_give_the_same_result() {
        let sc = Scoring::default();
        let mut state = 0x7777_1234_5678_9ABCu64;
        for round in 0..300 {
            let tlen = 200 + (xorshift(&mut state) % 1400) as usize;
            let t = random(&mut state, tlen);
            let qlen = 30 + (xorshift(&mut state) % 220) as usize;
            // Ein Teil der Reads stammt aus der Referenz, mit Abweichungen, Lücken und N.
            let mut q = if round % 3 == 0 || qlen + 10 > tlen {
                random(&mut state, qlen)
            } else {
                let at = (xorshift(&mut state) % (tlen - qlen) as u64) as usize;
                t[at..at + qlen].to_vec()
            };
            for _ in 0..(xorshift(&mut state) % 8) {
                let p = (xorshift(&mut state) % q.len() as u64) as usize;
                match xorshift(&mut state) % 4 {
                    0 => q[p] = 4,
                    1 => {
                        q.remove(p);
                    }
                    2 => q.insert(p, (xorshift(&mut state) % 4) as u8),
                    _ => q[p] = (q[p] + 1) & 3,
                }
            }
            if !rows_exact(q.len(), &sc) {
                continue;
            }
            let want = run(&q, &t, &sc, true, 19, 0x10000);
            let (gmax, te, qe, rows) = kernel_reference(&q, &t, &sc);
            let got = from_rows(gmax, te, qe, &rows, &sc, 19);
            assert_eq!(got, want, "Runde {round}");
            assert_eq!(
                with_start(&q, &t, &sc, true, 19, got),
                local_align(&q, &t, &sc, true, 19, true)
            );
        }
    }

    #[test]
    fn two_runs_side_by_side_equal_two_single_runs() {
        let sc = Scoring::default();
        let mut state = 0x2468_ACE0_1357_9BDFu64;
        for round in 0..400 {
            let n = |s: &mut u64, lo: u64, hi: u64| (lo + xorshift(s) % (hi - lo)) as usize;
            let qa = {
                let len = n(&mut state, 1, 200);
                random(&mut state, len)
            };
            // Gleich lange und verschieden lange Reads, gleich und verschieden lange Fenster,
            // auch leere Fenster und Reads, die ganz im Fenster stehen.
            let qb = if round % 3 == 0 {
                random(&mut state, qa.len())
            } else {
                {
                    let len = n(&mut state, 1, 200);
                    random(&mut state, len)
                }
            };
            let mut ta = {
                let len = n(&mut state, 0, 900);
                random(&mut state, len)
            };
            let tb = {
                let len = n(&mut state, 0, 900);
                random(&mut state, len)
            };
            if round % 5 == 0 && ta.len() > qa.len() {
                let at = ta.len() - qa.len();
                ta[at..].copy_from_slice(&qa);
            }
            for byte in [true, false] {
                for want_start in [true, false] {
                    let pair = local_align_pair((&qa, &ta), (&qb, &tb), &sc, byte, 19, want_start);
                    let single = (
                        local_align(&qa, &ta, &sc, byte, 19, want_start),
                        local_align(&qb, &tb, &sc, byte, 19, want_start),
                    );
                    assert_eq!(pair, single, "Runde {round}, byte {byte}");
                    let many = local_align_many(
                        [(&qa, &ta), (&qb, &tb), (&qb, &ta), (&qa, &tb)],
                        &sc,
                        byte,
                        19,
                        want_start,
                    );
                    let want = [
                        single.0,
                        single.1,
                        local_align(&qb, &ta, &sc, byte, 19, want_start),
                        local_align(&qa, &tb, &sc, byte, 19, want_start),
                    ];
                    assert_eq!(many, want, "Runde {round}, byte {byte}, vier");
                }
            }
        }
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
