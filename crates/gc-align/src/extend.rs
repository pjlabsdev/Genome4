//! Verlängern eines Seeds mit Smith-Waterman und affinen Lücken (BWA: `ksw_extend2`).
//!
//! Das Alignment ist am Seed verankert: Es beginnt mit dem Score `h0` des Seeds und läuft von
//! dort in eine Richtung. Gesucht ist die Stelle mit dem höchsten Score. Ein Pfad endet, sobald
//! sein Score auf 0 fällt. Lücken öffnen nur aus einem Treffer oder Fehltreffer heraus, nie
//! direkt aus einer anderen Lücke; ein Alignment wie `100M3I3D20M` entsteht so nicht.

use gc_sys::simd::Signed16;

/// Bewertung. Die Vorgaben sind die von BWA-MEM.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Scoring {
    /// Punkte für eine gleiche Base (`-A`).
    pub matching: i32,
    /// Abzug für eine ungleiche Base (`-B`), positiv angegeben.
    pub mismatch: i32,
    /// Abzug, wenn eine der beiden Basen N ist.
    pub ambiguous: i32,
    /// Öffnen und Verlängern einer Deletion (Basen der Referenz fehlen im Read; `-O`, `-E`).
    pub o_del: i32,
    pub e_del: i32,
    /// Öffnen und Verlängern einer Insertion (zusätzliche Basen im Read).
    pub o_ins: i32,
    pub e_ins: i32,
    /// Bandbreite (`-w`).
    pub band_width: i32,
    /// Abbruch, wenn der Score so weit unter das bisherige Maximum fällt (`-d`). 0 schaltet ab.
    pub zdrop: i32,
}

impl Default for Scoring {
    fn default() -> Self {
        Self {
            matching: 1,
            mismatch: 4,
            ambiguous: 1,
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            band_width: 100,
            zdrop: 100,
        }
    }
}

impl Scoring {
    /// Punkte für das Basenpaar (Codes 0..3, alles andere N).
    #[inline]
    pub fn score(&self, a: u8, b: u8) -> i32 {
        if a > 3 || b > 3 {
            -self.ambiguous
        } else if a == b {
            self.matching
        } else {
            -self.mismatch
        }
    }
}

/// Ergebnis einer Verlängerung.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Extension {
    /// Höchster Score (mindestens `h0`).
    pub score: i32,
    /// Verbrauchte Basen des Reads und der Referenz an der Stelle des höchsten Scores.
    pub query_len: i32,
    pub target_len: i32,
    /// Höchster Score, bei dem der Read bis zu seinem Ende ausgerichtet ist, und die dabei
    /// verbrauchten Basen der Referenz. -1 und 0, wenn das Ende nicht erreicht wird.
    pub global_score: i32,
    pub global_target_len: i32,
    /// Größter Abstand von der Diagonale, an dem ein neues Maximum lag.
    pub max_off: i32,
}

/// Verlängert ab dem Anker mit Score `h0` über `query` (Read) und `target` (Referenz).
/// `end_bonus` geht nur in die Begrenzung der Bandbreite ein, wie in bwa.
///
/// Die Rechnung endet, sobald feststeht, dass keine weitere Zeile der Referenz das Ergebnis noch
/// ändern kann (siehe `nothing_left_to_gain`). bwa rechnet in dem Fall weiter, bis der Score um
/// `zdrop` gefallen ist, und kommt zum selben Ergebnis.
///
/// Gerechnet wird, wenn die Zahlen in 16 Bit passen (bei Reads bis einige tausend Basen immer),
/// zeilenweise in Durchgängen, aus denen der Übersetzer Vektorbefehle macht (`extend_rows`);
/// sonst Zelle für Zelle (`extend_cells`). Beide liefern dasselbe.
pub fn extend(query: &[u8], target: &[u8], sc: &Scoring, end_bonus: i32, h0: i32) -> Extension {
    if fits_16_bit(query.len(), sc, h0) {
        extend_rows(query, target, sc, end_bonus, h0, true)
    } else {
        extend_cells(query, target, sc, end_bonus, h0, true)
    }
}

/// Zellen je Block in `extend_rows`: acht Zahlen zu 16 Bit füllen ein Vektorregister von 128 Bit,
/// das jeder unterstützte Prozessor hat.
const LANES: usize = 8;

/// Bleiben alle Zwischenwerte von `extend_rows` sicher in 16 Bit? Der größte ist der höchste
/// mögliche Score plus der Ausgleich für die Insertionsfolge, der kleinste liegt um dieselbe
/// Größe unter dem Startwert der Folge (-16384).
fn fits_16_bit(qlen: usize, sc: &Scoring, h0: i32) -> bool {
    let costs = [
        sc.matching,
        sc.mismatch,
        sc.ambiguous,
        sc.o_del,
        sc.e_del,
        sc.o_ins,
        sc.e_ins,
    ];
    if costs.iter().any(|&c| !(0..=1000).contains(&c)) || !(0..=16_000).contains(&h0) {
        return false;
    }
    let q = qlen as i64;
    let span = i64::from(h0)
        + (q + LANES as i64) * i64::from(sc.matching + sc.e_ins)
        + i64::from(sc.o_ins + sc.e_ins + sc.o_del + sc.e_del + sc.mismatch + sc.ambiguous);
    span < 16_000
}

/// Wie [`extend_cells`], aber jede Zeile in getrennten, einfachen Durchgängen über Felder von
/// 16-Bit-Zahlen, damit der Übersetzer Vektorbefehle daraus macht. Die einzige Größe, die in
/// einer Zeile von Zelle zu Zelle weiterläuft, ist die Insertion `F`: Sie ist das Maximum über
/// alle Zellen links davon, je Schritt um `e_ins` gesunken. Mit dem Ausgleich `ramp[k] = k *
/// e_ins` wird daraus ein laufendes Maximum, ein einziger Rechenschritt je Zelle; alles andere
/// hängt nur von der vorigen Zeile ab.
fn extend_rows(
    query: &[u8],
    target: &[u8],
    sc: &Scoring,
    end_bonus: i32,
    h0: i32,
    stop_early: bool,
) -> Extension {
    if crate::local::vectors_on() {
        extend_rows_v::<gc_sys::simd::I16x8>(query, target, sc, end_bonus, h0, stop_early)
    } else {
        extend_rows_v::<gc_sys::simd::plain::I16x8>(query, target, sc, end_bonus, h0, stop_early)
    }
}

fn extend_rows_v<V: Signed16>(
    query: &[u8],
    target: &[u8],
    sc: &Scoring,
    end_bonus: i32,
    h0: i32,
    stop_early: bool,
) -> Extension {
    let n_q = query.len();
    let qlen = n_q as i32;
    let tlen = target.len() as i32;
    let (e_del, e_ins) = (sc.e_del as i16, sc.e_ins as i16);
    let oe_del = (sc.o_del + sc.e_del) as i16;
    let oe_ins = (sc.o_ins + sc.e_ins) as i16;
    // Je Spalte j: h = H(i-1, j-1), e = E(i, j).
    let mut h = vec![0i16; n_q + 2 + LANES];
    let mut e = vec![0i16; n_q + 2 + LANES];
    h[0] = h0 as i16;
    if qlen >= 1 {
        h[1] = (h0 - i32::from(oe_ins)).max(0) as i16;
    }
    let mut j = 2;
    while j <= n_q && h[j - 1] > e_ins {
        h[j] = h[j - 1] - e_ins;
        j += 1;
    }
    let limit = |open: i32, ext: i32| -> i32 {
        ((f64::from(qlen * sc.matching + end_bonus - open) / f64::from(ext) + 1.0) as i32).max(1)
    };
    let w = sc
        .band_width
        .min(limit(sc.o_ins, sc.e_ins))
        .min(limit(sc.o_del, sc.e_del));

    // Alle Felder haben `LANES` Plätze Luft, damit ganze Blöcke gelesen werden können.
    let stride = n_q + LANES;
    let mut profile = vec![0i16; 5 * stride];
    for (a, row) in profile.chunks_exact_mut(stride).enumerate() {
        for (p, &q) in row.iter_mut().zip(query) {
            *p = sc.score(a as u8, q) as i16;
        }
    }
    let ramp: Vec<i16> = (0..=stride).map(|k| (k as i32 * sc.e_ins) as i16).collect();
    // Arbeitsfelder einer Zeile: bester Wert ohne Insertion, neue Deletion, laufendes Maximum,
    // fertige Zelle.
    let mut top = vec![0i16; stride];
    let mut del = vec![0i16; stride];
    let mut run = vec![0i16; stride];
    let mut out = vec![0i16; stride];

    let (mut max, mut max_i, mut max_j) = (h0, -1i32, -1i32);
    let (mut max_ie, mut gscore, mut max_off) = (-1i32, -1i32, 0i32);
    let (mut beg, mut end) = (0i32, qlen);
    for i in 0..tlen {
        let t = usize::from(target[i as usize].min(4));
        beg = beg.max(i - w);
        end = end.min(i + w + 1).min(qlen);
        let first = if beg == 0 {
            (h0 - (sc.o_del + sc.e_del * (i + 1))).max(0) as i16
        } else {
            0
        };
        let (mut h1, mut m, mut mj) = (first, 0i16, -1i32);
        if beg < end {
            let (b, n) = (beg as usize, end as usize);
            let len = n - b;
            // Gerechnet wird in ganzen Blöcken von `LANES` Zellen, auch über das Ende der Zeile
            // hinaus; die überzähligen Zellen landen nur in den Arbeitsfeldern und werden nicht
            // zurückgeschrieben.
            let width = len.div_ceil(LANES) * LANES;
            let points = &profile[t * stride + b..][..width];
            // 1. Paar-Score aus der Diagonale (tot bleibt tot), bester Wert ohne Insertion,
            //    Deletion für die nächste Zeile und der ausgeglichene Wert für die Insertion.
            let (v_e_del, v_oe_del) = (V::splat(e_del), V::splat(oe_del));
            let (v_oe_ins, v_e_ins, vz) = (V::splat(oe_ins), V::splat(e_ins), V::zero());
            for c in (0..width).step_by(LANES) {
                let d = V::load(&h[b + c..]);
                let p = V::load(&points[c..]);
                let ej = V::load(&e[b + c..]);
                let up = V::load(&ramp[c..]);
                // Tot bleibt tot: 0, wo die Diagonale 0 ist.
                let v = d.add_w(p).keep_if_nonzero(d);
                v.max(ej).store(&mut top[c..]);
                ej.sub_w(v_e_del)
                    .max(v.sub_w(v_oe_del))
                    .max(vz)
                    .store(&mut del[c..]);
                v.sub_w(v_oe_ins).add_w(up).store(&mut run[c..]);
            }
            // 2. Insertion: laufendes Maximum der ausgeglichenen Werte links von der Zelle, je
            //    Block in drei Schritten (um 1, 2 und 4 Spuren verschoben) und mit dem Übertrag
            //    aus dem Block davor. Werte hinter `len` gehen nicht in die Zeile ein.
            let mut carry = V::splat(i16::MIN / 2);
            for c in (0..width).step_by(LANES) {
                let x = V::load(&run[c..]);
                let mut incl = x.max(x.shift_fill(1, carry));
                incl = incl.max(incl.shift_fill(2, carry));
                incl = incl.max(incl.shift_fill(4, carry));
                incl.shift_fill(1, carry).store(&mut run[c..]);
                carry = V::splat(incl.hmax());
            }
            // 3. Zelle fertig: H(i, j).
            for c in (0..width).step_by(LANES) {
                let (x, r, up) = (V::load(&top[c..]), V::load(&run[c..]), V::load(&ramp[c..]));
                x.max(r.sub_w(up).add_w(v_e_ins))
                    .max(vz)
                    .store(&mut out[c..]);
            }
            // 4. Zurückschreiben: H um eins versetzt als Diagonale der nächsten Zeile.
            h[b + 1..=n].copy_from_slice(&out[..len]);
            e[b..n].copy_from_slice(&del[..len]);
            h[b] = first;
            // 5. Höchster Wert der Zeile und die letzte Spalte, in der er steht.
            let row = &h[b + 1..=n];
            m = row.iter().copied().fold(0, i16::max);
            mj = (b + row.iter().rposition(|&x| x == m).unwrap_or(len - 1)) as i32;
            h1 = row[len - 1];
        } else {
            h[end as usize] = h1;
        }
        e[end as usize] = 0;
        let (h1, m) = (i32::from(h1), i32::from(m));
        if end.max(beg) == qlen {
            if gscore <= h1 {
                max_ie = i;
            }
            gscore = gscore.max(h1);
        }
        if m == 0 {
            break;
        }
        if m > max {
            max = m;
            max_i = i;
            max_j = mj;
            max_off = max_off.max((mj - i).abs());
        } else if sc.zdrop > 0 {
            let (di, dj) = (i - max_i, mj - max_j);
            let gap = if di > dj {
                (di - dj) * sc.e_del
            } else {
                (dj - di) * sc.e_ins
            };
            if max - m - gap > sc.zdrop {
                break;
            }
        }
        if stop_early && gscore >= 0 {
            let rest = (qlen - 1 - mj) * sc.matching;
            let bound = max.min(gscore - 1);
            if (rest == 0 || m + rest <= bound)
                && nothing_left_to_gain(
                    &h,
                    &e,
                    (beg as usize, end as usize),
                    n_q,
                    sc.matching,
                    (max, gscore),
                )
            {
                break;
            }
        }
        // Band für die nächste Zeile auf die noch lebenden Spalten einengen.
        let mut j = beg;
        while j < end && h[j as usize] == 0 && e[j as usize] == 0 {
            j += 1;
        }
        beg = j;
        let mut j = end;
        while j >= beg && h[j as usize] == 0 && e[j as usize] == 0 {
            j -= 1;
        }
        end = (j + 2).min(qlen);
    }
    Extension {
        score: max,
        query_len: max_j + 1,
        target_len: max_i + 1,
        global_score: gscore,
        global_target_len: max_ie + 1,
        max_off,
    }
}

/// Obergrenze für alles, was nach der fertigen Zeile noch kommen kann, gegen das schon
/// Erreichte. Jeder Pfad in spätere Zeilen kommt aus einer lebenden Zelle dieser Zeile: schräg
/// aus `h[j]` in die Spalte `j` oder senkrecht über `e[j]`. Von dort kann er höchstens noch für
/// jede übrige Base des Reads einen Treffer bekommen; Lücken und Fehltreffer kosten nur, und eine
/// Zelle, die auf 0 gefallen ist, lebt nicht wieder auf. Liegt diese Grenze für alle Zellen nicht
/// über dem bisherigen Maximum und unter dem besten Score am Read-Ende, dann ändert sich keine
/// der gemeldeten Größen mehr: Ein neues Maximum bräuchte einen höheren Score, ein neues Ende
/// einen mindestens gleichen.
#[inline]
fn nothing_left_to_gain<T: Copy + Into<i32>>(
    h: &[T],
    e: &[T],
    (beg, end): (usize, usize),
    qlen: usize,
    matching: i32,
    (max, gscore): (i32, i32),
) -> bool {
    let limit = max.min(gscore - 1);
    // h[qlen] wäre der Weg in eine Spalte hinter dem Read-Ende, die es nicht gibt.
    for j in beg..=end.min(qlen.saturating_sub(1)) {
        let left = (qlen - j) as i32;
        let (hj, ej): (i32, i32) = (h[j].into(), e[j].into());
        if hj > 0 && hj + left * matching > limit {
            return false;
        }
        if ej > 0 && ej + (left - 1) * matching > limit {
            return false;
        }
    }
    true
}

fn extend_cells(
    query: &[u8],
    target: &[u8],
    sc: &Scoring,
    end_bonus: i32,
    h0: i32,
    stop_early: bool,
) -> Extension {
    let qlen = query.len() as i32;
    let tlen = target.len() as i32;
    let oe_del = sc.o_del + sc.e_del;
    let oe_ins = sc.o_ins + sc.e_ins;
    // Je Spalte j: h = H(i-1, j-1), e = E(i, j).
    let mut h = vec![0i32; query.len() + 2];
    let mut e = vec![0i32; query.len() + 2];
    h[0] = h0;
    if qlen >= 1 {
        h[1] = (h0 - oe_ins).max(0);
    }
    let mut j = 2;
    while j <= query.len() && h[j - 1] > sc.e_ins {
        h[j] = h[j - 1] - sc.e_ins;
        j += 1;
    }
    // Das Band braucht nie breiter zu sein als die längste Lücke, die sich noch lohnen kann.
    let best = sc.matching;
    let limit = |open: i32, ext: i32| -> i32 {
        ((f64::from(qlen * best + end_bonus - open) / f64::from(ext) + 1.0) as i32).max(1)
    };
    let w = sc
        .band_width
        .min(limit(sc.o_ins, sc.e_ins))
        .min(limit(sc.o_del, sc.e_del));

    // Punkte je Referenzbase gegen jede Read-Position, einmal vorab (wie das `qp` in bwa).
    let mut profile = vec![0i32; 5 * query.len()];
    for (a, row) in profile.chunks_exact_mut(query.len().max(1)).enumerate() {
        for (p, &q) in row.iter_mut().zip(query) {
            *p = sc.score(a as u8, q);
        }
    }

    let (mut max, mut max_i, mut max_j) = (h0, -1i32, -1i32);
    let (mut max_ie, mut gscore, mut max_off) = (-1i32, -1i32, 0i32);
    let (mut beg, mut end) = (0i32, qlen);
    for i in 0..tlen {
        let t = usize::from(target[i as usize].min(4));
        let (mut f, mut m, mut mj) = (0i32, 0i32, -1i32);
        beg = beg.max(i - w);
        end = end.min(i + w + 1).min(qlen);
        let mut h1 = if beg == 0 {
            (h0 - (sc.o_del + sc.e_del * (i + 1))).max(0)
        } else {
            0
        };
        let mut j = beg;
        if beg < end {
            let (b, n) = (beg as usize, end as usize);
            let points = &profile[t * query.len()..][b..n];
            for ((hj, ej), &pt) in h[b..n].iter_mut().zip(&mut e[b..n]).zip(points) {
                let diag = *hj;
                *hj = h1;
                let mm = if diag != 0 { diag + pt } else { 0 };
                let hh = mm.max(*ej).max(f);
                h1 = hh;
                if m <= hh {
                    mj = j;
                    m = hh;
                }
                *ej = (*ej - sc.e_del).max(mm - oe_del).max(0);
                f = (f - sc.e_ins).max(mm - oe_ins).max(0);
                j += 1;
            }
        }
        h[end as usize] = h1;
        e[end as usize] = 0;
        if j == qlen {
            if gscore <= h1 {
                max_ie = i;
            }
            gscore = gscore.max(h1);
        }
        if m == 0 {
            break;
        }
        if m > max {
            max = m;
            max_i = i;
            max_j = mj;
            max_off = max_off.max((mj - i).abs());
        } else if sc.zdrop > 0 {
            let (di, dj) = (i - max_i, mj - max_j);
            let gap = if di > dj {
                (di - dj) * sc.e_del
            } else {
                (dj - di) * sc.e_ins
            };
            if max - m - gap > sc.zdrop {
                break;
            }
        }
        // Vorprüfung ohne Durchgang: Die beste Zelle der Zeile selbst könnte noch gewinnen.
        if stop_early && gscore >= 0 {
            let rest = (qlen - 1 - mj) * sc.matching;
            let limit = max.min(gscore - 1);
            if (rest == 0 || m + rest <= limit)
                && nothing_left_to_gain(
                    &h,
                    &e,
                    (beg as usize, end as usize),
                    query.len(),
                    sc.matching,
                    (max, gscore),
                )
            {
                break;
            }
        }
        // Band für die nächste Zeile auf die noch lebenden Spalten einengen.
        let mut j = beg;
        while j < end && h[j as usize] == 0 && e[j as usize] == 0 {
            j += 1;
        }
        beg = j;
        let mut j = end;
        while j >= beg && h[j as usize] == 0 && e[j as usize] == 0 {
            j -= 1;
        }
        end = (j + 2).min(qlen);
    }
    Extension {
        score: max,
        query_len: max_j + 1,
        target_len: max_i + 1,
        global_score: gscore,
        global_target_len: max_ie + 1,
        max_off,
    }
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

    /// Unabhängige Rechnung mit drei vollen Matrizen, ohne Band und ohne Abbruch.
    /// Zustände: M (Paar), E (Deletion), F (Insertion); tot ist jeder Wert bis 0.
    fn reference(query: &[u8], target: &[u8], sc: &Scoring, h0: i32) -> Extension {
        const DEAD: i32 = i32::MIN / 4;
        let (n, m) = (target.len(), query.len());
        let live = |v: i32| if v > 0 { v } else { DEAD };
        // Index +1: Zeile 0 und Spalte 0 sind der Rand vor der ersten Base.
        let mut mm = vec![vec![DEAD; m + 1]; n + 1];
        let mut ee = vec![vec![DEAD; m + 1]; n + 1];
        let mut ff = vec![vec![DEAD; m + 1]; n + 1];
        mm[0][0] = h0;
        for (j, cell) in ff[0].iter_mut().enumerate().skip(1) {
            *cell = live(h0 - sc.o_ins - sc.e_ins * j as i32);
        }
        for (i, row) in ee.iter_mut().enumerate().skip(1) {
            row[0] = live(h0 - sc.o_del - sc.e_del * i as i32);
        }
        let best3 = |a: i32, b: i32, c: i32| a.max(b).max(c);
        let (mut max, mut max_i, mut max_j) = (h0, 0usize, 0usize);
        let (mut gscore, mut gi) = (-1, 0usize);
        for i in 1..=n {
            let mut row_max = 0;
            let mut row_j = 0;
            for j in 1..=m {
                let prev = best3(mm[i - 1][j - 1], ee[i - 1][j - 1], ff[i - 1][j - 1]);
                if prev > 0 {
                    mm[i][j] = prev + sc.score(target[i - 1], query[j - 1]);
                }
                let open_e = if mm[i - 1][j] > 0 {
                    mm[i - 1][j] - sc.o_del - sc.e_del
                } else {
                    DEAD
                };
                ee[i][j] = live(open_e.max(ee[i - 1][j] - sc.e_del));
                let open_f = if mm[i][j - 1] > 0 {
                    mm[i][j - 1] - sc.o_ins - sc.e_ins
                } else {
                    DEAD
                };
                ff[i][j] = live(open_f.max(ff[i][j - 1] - sc.e_ins));
                let hh = best3(mm[i][j], ee[i][j], ff[i][j]).max(0);
                if hh >= row_max {
                    row_max = hh;
                    row_j = j;
                }
                if j == m && hh >= gscore {
                    gscore = hh;
                    gi = i;
                }
            }
            if row_max > max {
                max = row_max;
                max_i = i;
                max_j = row_j;
            }
        }
        Extension {
            score: max,
            query_len: max_j as i32,
            target_len: max_i as i32,
            global_score: gscore,
            global_target_len: gi as i32,
            max_off: 0,
        }
    }

    fn mutate(rng: &mut u64, src: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        for &b in src {
            match xorshift(rng) % 40 {
                0 => out.push((xorshift(rng) & 3) as u8),
                1 => {}
                2 => {
                    out.push(b);
                    for _ in 0..1 + xorshift(rng) % 4 {
                        out.push((xorshift(rng) & 3) as u8);
                    }
                }
                3 if xorshift(rng).is_multiple_of(4) => out.push(4),
                _ => out.push(b),
            }
        }
        out
    }

    #[test]
    fn equals_full_matrix_without_band_and_zdrop() {
        let mut rng = 0x0123_4567_89AB_CDEF_u64;
        let sc = Scoring {
            band_width: 10_000,
            zdrop: 0,
            ..Scoring::default()
        };
        let mut reached_end = 0;
        for round in 0..3000 {
            let tlen = 1 + (xorshift(&mut rng) % 70) as usize;
            let target: Vec<u8> = (0..tlen).map(|_| (xorshift(&mut rng) & 3) as u8).collect();
            let mut query = mutate(&mut rng, &target);
            if xorshift(&mut rng).is_multiple_of(3) {
                let cut = (xorshift(&mut rng) as usize) % (query.len() + 1);
                query.truncate(cut);
            }
            if query.is_empty() {
                continue;
            }
            let h0 = 1 + (xorshift(&mut rng) % 40) as i32;
            // end_bonus groß genug, dass die Begrenzung des Bandes nie greift
            let got = extend(&query, &target, &sc, 10_000, h0);
            let want = reference(&query, &target, &sc, h0);
            let got_cmp = Extension { max_off: 0, ..got };
            // Erreicht kein Pfad das Read-Ende, meldet die Bandrechnung -1 oder 0; beides heißt
            // "nicht erreicht".
            if want.global_score <= 0 {
                assert!(got.global_score <= 0, "Runde {round}");
                assert_eq!(
                    (got.score, got.query_len, got.target_len),
                    (want.score, want.query_len, want.target_len),
                    "Runde {round} q={query:?} t={target:?} h0={h0}"
                );
            } else {
                assert_eq!(
                    got_cmp, want,
                    "Runde {round} q={query:?} t={target:?} h0={h0}"
                );
                reached_end += 1;
            }
        }
        assert!(reached_end > 1000);
    }

    #[test]
    fn stopping_early_changes_nothing() {
        let mut rng = 0x5EED_1234_ABCD_u64;
        let (mut stopped_rows, mut cases) = (0, 0);
        for round in 0..20_000u32 {
            let sc = Scoring {
                band_width: [100, 200, 3, 10][(round % 4) as usize],
                zdrop: [100, 100, 0, 20][((round / 4) % 4) as usize],
                matching: 1 + (round / 16 % 2) as i32,
                ..Scoring::default()
            };
            // Referenz länger als der Read, wie im Aligner: Platz für Lücken.
            let qlen = (xorshift(&mut rng) % 90) as usize;
            let tlen = qlen + (xorshift(&mut rng) % 120) as usize;
            let target: Vec<u8> = (0..tlen).map(|_| (xorshift(&mut rng) & 3) as u8).collect();
            let mut query = match round % 3 {
                // gleich, mit wenigen Fehlern, oder mit Lücken und N
                0 => target.clone(),
                1 => {
                    let mut q = target.clone();
                    for _ in 0..xorshift(&mut rng) % 4 {
                        if !q.is_empty() {
                            let at = (xorshift(&mut rng) as usize) % q.len();
                            q[at] = (xorshift(&mut rng) & 3) as u8;
                        }
                    }
                    q
                }
                _ => mutate(&mut rng, &target),
            };
            query.truncate(qlen);
            let h0 = 1 + (xorshift(&mut rng) % 150) as i32;
            let end_bonus = [5, 0, 10_000][(round % 3) as usize];
            let plain = extend_cells(&query, &target, &sc, end_bonus, h0, false);
            let early = extend_cells(&query, &target, &sc, end_bonus, h0, true);
            assert_eq!(
                early, plain,
                "Runde {round} q={query:?} t={target:?} h0={h0} {sc:?}"
            );
            // Die Rechnung in Zeilen mit 16 Bit liefert dasselbe, mit und ohne frühes Ende.
            assert!(fits_16_bit(query.len(), &sc, h0));
            for stop in [false, true] {
                assert_eq!(
                    extend_rows(&query, &target, &sc, end_bonus, h0, stop),
                    plain,
                    "Zeilen, Runde {round} q={query:?} t={target:?} h0={h0} {sc:?}"
                );
            }
            cases += 1;
            stopped_rows += usize::from(plain.global_score > 0);
        }
        assert!(cases == 20_000 && stopped_rows > 5_000);
    }

    #[test]
    fn large_numbers_take_the_cell_by_cell_path() {
        let sc = Scoring::default();
        assert!(fits_16_bit(250, &sc, 250));
        assert!(!fits_16_bit(20_000, &sc, 19));
        assert!(!fits_16_bit(100, &sc, 20_000));
        let wide = Scoring {
            matching: 2000,
            ..Scoring::default()
        };
        assert!(!fits_16_bit(10, &wide, 19));
        // Beide Wege geben dasselbe, wo beide gehen; der lange Read geht den Weg der Zellen.
        let long: Vec<u8> = (0..9000).map(|i| (i * 7 % 4) as u8).collect();
        let x = extend(&long, &long, &sc, 5, 19);
        assert_eq!((x.score, x.query_len, x.target_len), (9019, 9000, 9000));
    }

    #[test]
    fn band_and_zdrop_only_restrict() {
        let mut rng = 0xFEED_BEEF_u64;
        let free = Scoring {
            band_width: 10_000,
            zdrop: 0,
            ..Scoring::default()
        };
        let narrow = Scoring {
            band_width: 3,
            zdrop: 20,
            ..Scoring::default()
        };
        let mut same = 0;
        for _ in 0..1000 {
            let tlen = 20 + (xorshift(&mut rng) % 100) as usize;
            let target: Vec<u8> = (0..tlen).map(|_| (xorshift(&mut rng) & 3) as u8).collect();
            let query = mutate(&mut rng, &target);
            let a = extend(&query, &target, &free, 10_000, 20);
            let b = extend(&query, &target, &narrow, 5, 20);
            assert!(b.score <= a.score);
            assert!(b.score >= 20);
            same += usize::from(a.score == b.score);
        }
        assert!(same > 500);
    }

    #[test]
    fn simple_cases() {
        let sc = Scoring::default();
        // gleiche Folgen: jede Base ein Punkt
        let s = [0u8, 1, 2, 3, 0, 1, 2, 3, 2, 1];
        let x = extend(&s, &s, &sc, 5, 19);
        assert_eq!((x.score, x.query_len, x.target_len), (29, 10, 10));
        assert_eq!((x.global_score, x.global_target_len), (29, 10));
        // ein Fehltreffer am Ende lohnt nicht: Maximum davor, Ende trotzdem erreichbar
        let mut q = s;
        q[9] = 0;
        let x = extend(&q, &s, &sc, 5, 19);
        assert_eq!((x.score, x.query_len, x.target_len), (28, 9, 9));
        assert_eq!((x.global_score, x.global_target_len), (24, 10));
        // leere Referenz
        let x = extend(&s, &[], &sc, 5, 19);
        assert_eq!((x.score, x.query_len, x.target_len), (19, 0, 0));
    }
}
