//! Verlängern eines Seeds mit Smith-Waterman und affinen Lücken (BWA: `ksw_extend2`).
//!
//! Das Alignment ist am Seed verankert: Es beginnt mit dem Score `h0` des Seeds und läuft von
//! dort in eine Richtung. Gesucht ist die Stelle mit dem höchsten Score. Ein Pfad endet, sobald
//! sein Score auf 0 fällt. Lücken öffnen nur aus einem Treffer oder Fehltreffer heraus, nie
//! direkt aus einer anderen Lücke; ein Alignment wie `100M3I3D20M` entsteht so nicht.

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
pub fn extend(query: &[u8], target: &[u8], sc: &Scoring, end_bonus: i32, h0: i32) -> Extension {
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
