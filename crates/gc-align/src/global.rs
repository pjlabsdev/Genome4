//! Globales Alignment im Band mit affinen Lücken und Rückverfolgung (BWA: `ksw_global2`), und
//! daraus das CIGAR eines Bereichs (`bwa_gen_cigar2`).
//!
//! Anders als beim Verlängern steht hier fest, welches Stück Read auf welches Stück Referenz
//! gehört. Beide Enden sind festgelegt; gesucht ist der beste Weg dazwischen und seine Schritte.

use crate::extend::Scoring;
use gc_index::FmIndex;

/// Ein Schritt im CIGAR.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Base gegen Base, gleich oder ungleich (`M`).
    Match,
    /// Zusätzliche Basen im Read (`I`).
    Ins,
    /// Basen der Referenz fehlen im Read (`D`).
    Del,
    /// Read-Ende, das nicht zum Alignment gehört (`S`).
    SoftClip,
}

impl Op {
    pub fn letter(self) -> char {
        match self {
            Op::Match => 'M',
            Op::Ins => 'I',
            Op::Del => 'D',
            Op::SoftClip => 'S',
        }
    }
}

/// CIGAR als Folge (Schritt, Länge).
pub type Cigar = Vec<(Op, u32)>;

/// Schreibt ein CIGAR in der üblichen Form, etwa `57M36D93M`.
pub fn cigar_string(cigar: &[(Op, u32)]) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    for &(op, n) in cigar {
        let _ = write!(s, "{n}{}", op.letter());
    }
    s
}

const MINUS_INF: i32 = -0x4000_0000;

fn push(cigar: &mut Cigar, op: Op, n: u32) {
    match cigar.last_mut() {
        Some((last, len)) if *last == op => *len += n,
        _ => cigar.push((op, n)),
    }
}

/// Globales Alignment von `query` gegen `target` innerhalb des Bands `w` um die Hauptdiagonale.
/// Liefert den Score und das CIGAR. Bei gleich guten Wegen gilt die Reihenfolge von bwa:
/// Base gegen Base vor Deletion vor Insertion. Eine Lücke öffnet nur aus einem Basenpaar heraus.
pub fn global(query: &[u8], target: &[u8], sc: &Scoring, w: i32) -> (i32, Cigar) {
    let qlen = query.len();
    let tlen = target.len();
    let w = w.max(0) as usize;
    let oe_del = sc.o_del + sc.e_del;
    let oe_ins = sc.o_ins + sc.e_ins;
    let n_col = qlen.min(2 * w + 1);
    // Je Zelle: Richtung für H (Bits 0..1), für E der nächsten Zeile (2..3), für F der
    // nächsten Spalte (4..5).
    let mut z = vec![0u8; n_col * tlen];
    let mut h_row = vec![MINUS_INF; qlen + 1];
    let mut e_row = vec![MINUS_INF; qlen + 1];
    h_row[0] = 0;
    for (j, h) in h_row.iter_mut().enumerate().skip(1).take(w.min(qlen)) {
        *h = -(sc.o_ins + sc.e_ins * j as i32);
    }
    for i in 0..tlen {
        let beg = i.saturating_sub(w);
        let end = (i + w + 1).min(qlen);
        let mut f = MINUS_INF;
        let mut h1 = if beg == 0 {
            -(sc.o_del + sc.e_del * (i as i32 + 1))
        } else {
            MINUS_INF
        };
        let zi = &mut z[i * n_col..];
        for j in beg..end {
            let m = h_row[j] + sc.score(target[i], query[j]);
            let mut e = e_row[j];
            h_row[j] = h1;
            let mut d: u8 = if m >= e { 0 } else { 1 };
            let mut h = m.max(e);
            if h < f {
                d = 2;
                h = f;
            }
            h1 = h;
            let t = m - oe_del;
            e -= sc.e_del;
            if e > t {
                d |= 1 << 2;
            } else {
                e = t;
            }
            e_row[j] = e;
            let t = m - oe_ins;
            f -= sc.e_ins;
            if f > t {
                d |= 2 << 4;
            } else {
                f = t;
            }
            zi[j - beg] = d;
        }
        h_row[end] = h1;
        e_row[end] = MINUS_INF;
    }
    let score = h_row[qlen];
    let mut cigar = Cigar::new();
    let mut i = tlen as i64 - 1;
    let mut k = (tlen + w).min(qlen) as i64 - 1;
    let mut which = 0u8;
    while i >= 0 && k >= 0 {
        let beg = (i as usize).saturating_sub(w);
        let col = k as usize - beg;
        which = (z[i as usize * n_col + col] >> (which << 1)) & 3;
        match which {
            0 => {
                push(&mut cigar, Op::Match, 1);
                i -= 1;
                k -= 1;
            }
            1 => {
                push(&mut cigar, Op::Del, 1);
                i -= 1;
            }
            _ => {
                push(&mut cigar, Op::Ins, 1);
                k -= 1;
            }
        }
    }
    if i >= 0 {
        push(&mut cigar, Op::Del, i as u32 + 1);
    }
    if k >= 0 {
        push(&mut cigar, Op::Ins, k as u32 + 1);
    }
    cigar.reverse();
    (score, cigar)
}

/// Ergebnis von [`gen_cigar`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GlobalAlignment {
    pub score: i32,
    pub cigar: Cigar,
    /// Editierabstand: ungleiche Basen plus Basen in Lücken (`NM`). Eine Deletion ganz am Anfang
    /// oder Ende zählt nicht; sie wird später aus dem CIGAR entfernt.
    pub nm: u32,
    /// Die Referenz an den abweichenden Stellen (`MD`), in Richtung des Vorwärtsstrangs.
    pub md: String,
}

/// Alignment des Read-Stücks `query` gegen die Textpositionen `rb..re` (Vorwärtsstrang, dann
/// reverses Komplement). `None`, wenn das Stück leer ist oder über die Strangmitte reicht.
///
/// Liegt das Stück auf dem Gegenstrang, werden beide Folgen vor dem Alignment umgedreht, damit
/// Lücken in Wiederholungen auf dem Vorwärtsstrang so weit links wie möglich stehen. CIGAR und
/// `MD` sind dann schon in der Richtung des Vorwärtsstrangs.
pub fn gen_cigar(
    index: &FmIndex,
    query: &[u8],
    rb: u64,
    re: u64,
    sc: &Scoring,
    band: i32,
) -> Option<GlobalAlignment> {
    use std::fmt::Write;
    let l = index.reference().len();
    if query.is_empty() || rb >= re || (rb < l && re > l) || re > 2 * l {
        return None;
    }
    let mut rseq: Vec<u8> = (rb..re).map(|p| index.text_base(p)).collect();
    let mut q = query.to_vec();
    let reverse = rb >= l;
    if reverse {
        q.reverse();
        rseq.reverse();
    }
    let (score, cigar) = if q.len() == rseq.len() && band == 0 {
        let score = q.iter().zip(&rseq).map(|(&a, &b)| sc.score(b, a)).sum();
        (score, vec![(Op::Match, q.len() as u32)])
    } else {
        let half = (q.len() as i32 + 1) >> 1;
        let max_ins = (f64::from(half * sc.matching - sc.o_ins) / f64::from(sc.e_ins) + 1.0) as i32;
        let max_del = (f64::from(half * sc.matching - sc.o_del) / f64::from(sc.e_del) + 1.0) as i32;
        let max_gap = max_ins.max(max_del).max(1);
        let diff = (rseq.len() as i32 - q.len() as i32).abs();
        let w = ((max_gap + diff + 1) >> 1).min(band).max(diff + 3);
        global(&q, &rseq, sc, w)
    };
    // Umgedreht auf dem Gegenstrang heißt: Komplement des Vorwärtsstrangs.
    let letters: &[u8; 5] = if reverse { b"TGCAN" } else { b"ACGTN" };
    let letter = |code: u8| letters[usize::from(code.min(4))] as char;
    let (mut x, mut y, mut run) = (0usize, 0usize, 0u32);
    let (mut mismatches, mut gaps) = (0u32, 0u32);
    let mut md = String::new();
    for (k, &(op, n)) in cigar.iter().enumerate() {
        let n_us = n as usize;
        match op {
            Op::Match => {
                for i in 0..n_us {
                    if q[x + i] != rseq[y + i] {
                        let _ = write!(md, "{run}{}", letter(rseq[y + i]));
                        mismatches += 1;
                        run = 0;
                    } else {
                        run += 1;
                    }
                }
                x += n_us;
                y += n_us;
            }
            Op::Del => {
                if k > 0 && k < cigar.len() - 1 {
                    let _ = write!(md, "{run}^");
                    md.extend(rseq[y..y + n_us].iter().map(|&b| letter(b)));
                    run = 0;
                    gaps += n;
                }
                y += n_us;
            }
            Op::Ins => {
                x += n_us;
                gaps += n;
            }
            Op::SoftClip => x += n_us,
        }
    }
    let _ = write!(md, "{run}");
    Some(GlobalAlignment {
        score,
        cigar,
        nm: mismatches + gaps,
        md,
    })
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

    /// Unabhängige Rechnung ohne Band: Gotoh mit drei Tabellen, Lücke nur aus einem Basenpaar.
    fn reference_score(q: &[u8], t: &[u8], sc: &Scoring) -> i32 {
        let (n, m) = (t.len(), q.len());
        let inf = MINUS_INF;
        let mut mm = vec![vec![inf; m + 1]; n + 1];
        let mut e = vec![vec![inf; m + 1]; n + 1];
        let mut f = vec![vec![inf; m + 1]; n + 1];
        mm[0][0] = 0;
        for i in 0..=n {
            for j in 0..=m {
                if i > 0 {
                    e[i][j] = (mm[i - 1][j] - sc.o_del - sc.e_del).max(e[i - 1][j] - sc.e_del);
                }
                if j > 0 {
                    f[i][j] = (mm[i][j - 1] - sc.o_ins - sc.e_ins).max(f[i][j - 1] - sc.e_ins);
                }
                if i > 0 && j > 0 {
                    let best = mm[i - 1][j - 1].max(e[i - 1][j - 1]).max(f[i - 1][j - 1]);
                    mm[i][j] = best + sc.score(t[i - 1], q[j - 1]);
                }
            }
        }
        mm[n][m].max(e[n][m]).max(f[n][m])
    }

    fn cigar_score(q: &[u8], t: &[u8], cigar: &[(Op, u32)], sc: &Scoring) -> (i32, usize, usize) {
        let (mut x, mut y, mut s) = (0usize, 0usize, 0i32);
        for &(op, n) in cigar {
            let n_us = n as usize;
            match op {
                Op::Match => {
                    s += (0..n_us).map(|k| sc.score(t[y + k], q[x + k])).sum::<i32>();
                    x += n_us;
                    y += n_us;
                }
                Op::Ins => {
                    s -= sc.o_ins + sc.e_ins * n as i32;
                    x += n_us;
                }
                Op::Del => {
                    s -= sc.o_del + sc.e_del * n as i32;
                    y += n_us;
                }
                Op::SoftClip => x += n_us,
            }
        }
        (s, x, y)
    }

    #[test]
    fn known_cases_give_the_expected_cigar() {
        let sc = Scoring::default();
        let t: Vec<u8> = (0..60u32).map(|i| ((i * 7 + i / 3) % 4) as u8).collect();
        let (s, c) = global(&t, &t, &sc, 10);
        assert_eq!((s, cigar_string(&c).as_str()), (60, "60M"));
        let mut q = t.clone();
        q.drain(30..33);
        let (s, c) = global(&q, &t, &sc, 10);
        assert_eq!(s, 57 - 6 - 3);
        assert_eq!(
            c.iter().map(|x| x.0).collect::<Vec<_>>(),
            [Op::Match, Op::Del, Op::Match]
        );
        let mut q = t.clone();
        q.splice(20..20, [3, 3, 0, 3]);
        let (s, c) = global(&q, &t, &sc, 10);
        assert_eq!(s, 60 - 6 - 4);
        assert_eq!(
            c.iter().map(|x| x.0).collect::<Vec<_>>(),
            [Op::Match, Op::Ins, Op::Match]
        );
    }

    #[test]
    fn score_equals_an_unbanded_computation_and_the_cigar_adds_up_to_it() {
        let sc = Scoring::default();
        let mut state = 0x2545_f491_4f6c_dd1du64;
        for round in 0..400 {
            let tlen = 1 + (xorshift(&mut state) % 80) as usize;
            let t: Vec<u8> = (0..tlen)
                .map(|_| (xorshift(&mut state) % 4) as u8)
                .collect();
            let mut q = t.clone();
            for _ in 0..(xorshift(&mut state) % 5) {
                if q.is_empty() {
                    break;
                }
                let p = (xorshift(&mut state) as usize) % q.len();
                match xorshift(&mut state) % 3 {
                    0 => q[p] = (xorshift(&mut state) % 5) as u8,
                    1 => {
                        let n = 1 + (xorshift(&mut state) as usize % 4);
                        let end = (p + n).min(q.len());
                        q.drain(p..end);
                    }
                    _ => {
                        let n = 1 + (xorshift(&mut state) % 4);
                        for _ in 0..n {
                            q.insert(p, (xorshift(&mut state) % 4) as u8);
                        }
                    }
                }
            }
            if q.is_empty() {
                continue;
            }
            // Band so weit, dass es nichts abschneidet.
            let w = (tlen + q.len()) as i32;
            let (score, cigar) = global(&q, &t, &sc, w);
            assert_eq!(score, reference_score(&q, &t, &sc), "Runde {round}");
            let (s, x, y) = cigar_score(&q, &t, &cigar, &sc);
            assert_eq!((s, x, y), (score, q.len(), t.len()), "Runde {round}");
            // Ein enges Band kann nur schlechter sein, und sein CIGAR muss weiter aufgehen.
            let narrow = (tlen as i32 - q.len() as i32).abs() + 3;
            let (score_n, cigar_n) = global(&q, &t, &sc, narrow);
            assert!(score_n <= score);
            let (s, x, y) = cigar_score(&q, &t, &cigar_n, &sc);
            assert_eq!(
                (s, x, y),
                (score_n, q.len(), t.len()),
                "Runde {round}, enges Band"
            );
        }
    }
}
