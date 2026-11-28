//! Versuch Strobe-Paarweg (`GC_STROBE_PAARWEG=1`): Paare, bei denen der Strobeweg mindestens
//! einen Read nicht allein verorten kann, werden aus den Strobe-Ankern beider Reads gemeinsam
//! verortet, statt über den vollen Weg und die Partnersuche.
//!
//! Ablauf je Paar: Treffer beider Reads nach Diagonale bündeln (Anker), Ankerpaare in
//! Paarrichtung und erlaubtem Abstand nach der Zahl ihrer Treffer ordnen, die Anker der besten
//! Ankerpaare und die besten Einzelanker jedes Reads ausrichten (zuerst ohne Lücke, lokal nur bei
//! Verdacht). Die Bereiche gehen dann in die übliche Paarbewertung; eine Partnersuche gibt es für
//! diese Paare nicht. Hat ein Read keinen Anker oder findet sich kein Ankerpaar, geht das Paar den
//! bisherigen Weg.

use crate::region::Region;
use crate::sam::AlignOptions;
use gc_index::{FmIndex, StrobeIndex};

/// Schalter: `GC_STROBE_PAARWEG=1`.
pub fn mode() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_PAARWEG").is_ok_and(|v| v == "1"))
}

/// Ein schon verorteter Read bekommt die Orte seiner Anker dazu (Schalter zurück:
/// `GC_STROBE_PAAR_DAZU=0`).
pub fn merge_mode() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_PAAR_DAZU").map_or(true, |v| v != "0"))
}

/// Fast beste Orte ohne Partner schicken das Paar den bisherigen Weg (Schalter zurück:
/// `GC_STROBE_PAAR_OHNE_PARTNER=1`).
fn unpaired_mode() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_PAAR_OHNE_PARTNER").map_or(true, |v| v != "1"))
}

fn rescue_mode() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_PAAR_SUCHE").is_ok_and(|v| v == "1"))
}

/// Wie ein Paar aus den Ankern genommen wird.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// Jeder fast beste Ort hat seinen Partner: ohne Partnersuche.
    Sure,
    /// Mit Partnersuche, ohne vollen Weg.
    Rescue,
    /// Lage unsicher: bestes Ankerpaar, MAPQ höchstens [`unsure_cap`].
    Unsure,
}

/// Unsichere Paare selbst nehmen, mit kleiner MAPQ (`GC_STROBE_PAAR_UNSICHER=0..60`), statt
/// des bisherigen Wegs.
fn unsure_mode() -> bool {
    unsure_cap().is_some()
}

pub fn unsure_cap() -> Option<u8> {
    static V: std::sync::OnceLock<Option<u8>> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("GC_STROBE_PAAR_UNSICHER")
            .ok()
            .and_then(|v| v.parse().ok())
    })
}

pub static UNSURE: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

fn frac_mode() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_PAAR_WDH").map_or(true, |v| v != "0"))
}

fn rescue_needs_clean() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_PAAR_SUCHE_SAUBER").map_or(true, |v| v != "0"))
}

fn unique_mode() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_PAAR_EINDEUTIG").map_or(true, |v| v != "0"))
}

/// Unsichere Paare gehen den vollen Weg, aber ohne Partnersuche
/// (`GC_STROBE_PAAR_UNSICHER_VOLL=1`, dazu `GC_STROBE_PAAR_UNSICHER` setzen).
pub fn unsure_full() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("GC_STROBE_PAAR_UNSICHER_VOLL").is_ok_and(|v| v == "1"))
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

/// Größter Abstand eines Ankerpaares (Anfang des vorderen bis Ende des hinteren Reads).
fn pair_max() -> i64 {
    static V: std::sync::OnceLock<i64> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_usize("GC_STROBE_PAAR_MAX", 1500) as i64)
}

pub fn pair_max_pub() -> i64 {
    pair_max()
}

/// So viele der besten Ankerpaare werden ausgerichtet.
fn top_pairs() -> usize {
    static V: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_usize("GC_STROBE_PAAR_TOP", 16))
}

/// Dazu je Read so viele der besten Einzelanker.
fn top_single() -> usize {
    static V: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_usize("GC_STROBE_PAAR_EINZELN", 3))
}

/// Höchstens so viele Fehlpaarungen gelten ohne Lücke als sicher.
fn max_mismatch() -> usize {
    static V: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_usize("GC_STROBE_PAAR_FEHLER", 2))
}

/// Zähler für das Protokoll.
pub static TRIED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static TAKEN: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// (Strang, Diagonale) eines Ankers.
fn split(g: u64) -> (u64, i64) {
    (g >> 63, (g & !(1 << 63)) as i64 - (1 << 32))
}

/// Die Bereiche beider Reads aus ihren Ankern; `None`: bisheriger Weg.
pub fn pair_regions(
    index: &FmIndex,
    strobes: &StrobeIndex,
    reads: [&[u8]; 2],
    unique: [bool; 2],
    opt: &AlignOptions,
) -> Option<([Vec<Region>; 2], Kind, [f32; 2])> {
    use std::sync::atomic::Ordering::Relaxed;
    TRIED.fetch_add(1, Relaxed);
    if reads
        .iter()
        .any(|r| r.len() < 60 || r.iter().any(|&c| c > 3))
    {
        return None;
    }
    let full = reads.map(|r| crate::strobe::strobe_groups(strobes, r, &opt.strobe));
    // Anteil gesperrter Randstrobes als Wiederholungsanteil (wie bwa `frac_rep`), Schalter
    // `GC_STROBE_PAAR_WDH=0` setzt ihn auf 0.
    let frac = [0, 1].map(|r| {
        let (_, b, u) = full[r];
        if b + u == 0 || !frac_mode() {
            0.0
        } else {
            b as f32 / (b + u) as f32
        }
    });
    let groups = full.map(|g| g.0);
    if groups.iter().any(Vec::is_empty) {
        return None;
    }
    // Anker des zweiten Reads je Strang nach Diagonale.
    let mut by_diag: [Vec<(i64, usize)>; 2] = [Vec::new(), Vec::new()];
    for (k, &(_, g)) in groups[1].iter().enumerate() {
        let (s, d) = split(g);
        by_diag[s as usize].push((d, k));
    }
    for v in &mut by_diag {
        v.sort_unstable();
    }
    let (l0, l1) = (reads[0].len() as i64, reads[1].len() as i64);
    let max = pair_max();
    // (Treffer beider Anker, Anker Read 1, Anker Read 2)
    let mut cands: Vec<(u32, usize, usize)> = Vec::new();
    for (i, &(c0, g)) in groups[0].iter().enumerate().take(500) {
        let (s, p) = split(g);
        // Vorwärts liegender Read vorn, der andere auf dem Gegenstrang dahinter.
        let (lo, hi) = if s == 0 {
            (p - l1, p + max - l1)
        } else {
            (p + l0 - max, p + l0)
        };
        let v = &by_diag[1 - s as usize];
        let from = v.partition_point(|&(d, _)| d < lo);
        for &(d, k) in &v[from..] {
            if d > hi {
                break;
            }
            cands.push((c0 + groups[1][k].0, i, k));
        }
    }
    if cands.is_empty() {
        return None;
    }
    cands.sort_unstable_by(|a, b| b.0.cmp(&a.0).then((a.1, a.2).cmp(&(b.1, b.2))));
    // Alle Ankerpaare ab der Hälfte der besten Trefferzahl werden ausgerichtet; sind es zu
    // viele, könnte die wahre Lage verborgen bleiben: dann der bisherige Weg.
    let floor = (cands[0].0 / 2).max(2);
    let n_good = cands.iter().take_while(|c| c.0 >= floor).count();
    let crowded = n_good > top_pairs();
    if crowded && !unsure_mode() {
        return None;
    }
    let mut want: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
    for &(_, i, k) in cands.iter().take(n_good.clamp(1, top_pairs())) {
        if !want[0].contains(&i) {
            want[0].push(i);
        }
        if !want[1].contains(&k) {
            want[1].push(k);
        }
    }
    for (r, w) in want.iter_mut().enumerate() {
        for i in 0..groups[r].len().min(top_single()) {
            if !w.contains(&i) {
                w.push(i);
            }
        }
    }
    let mut out: [Vec<Region>; 2] = [Vec::new(), Vec::new()];
    for r in 0..2 {
        for &i in &want[r] {
            if let Some(reg) = align_anchor(index, reads[r], groups[r][i].1, opt) {
                // dieselbe Stelle nur einmal, die bessere bleibt
                // fast deckungsgleich (wie bwa: 95 % des kürzeren)
                let same = |o: &Region| {
                    let ov = o.re.min(reg.re) as i64 - o.rb.max(reg.rb) as i64;
                    let short = (o.re - o.rb).min(reg.re - reg.rb) as i64;
                    ov * 100 >= short * 95
                };
                match out[r].iter_mut().find(|o| same(o)) {
                    Some(o) if o.score >= reg.score => {}
                    Some(o) => *o = reg,
                    None => out[r].push(reg),
                }
            }
        }
        if out[r].is_empty() {
            return None;
        }
        out[r].sort_unstable_by(|a, b| b.score.cmp(&a.score).then(a.rb.cmp(&b.rb)));
    }
    // Jeder fast beste Ort eines Reads braucht einen Partner in Paarlage; sonst könnte dort ein
    // gleich gutes Paar liegen, das nur die Partnersuche fände: dann der bisherige Weg.
    if unpaired_mode() {
        let l = index.reference().len();
        // (Contig, Anfang vorwärts, Ende vorwärts, Gegenstrang)
        let at = |x: &Region| {
            if x.rb < l {
                (x.contig, x.rb as i64, x.re as i64, false)
            } else {
                (x.contig, (2 * l - x.re) as i64, (2 * l - x.rb) as i64, true)
            }
        };
        for r in 0..2 {
            // Ist der Partner für sich eindeutig, kann er bei den anderen Orten dieses Reads
            // höchstens mit seinem Zweitbesten liegen; das deckt die Paarbewertung über den
            // Abzug für ungepaarte Reads ab (Schalter zurück: `GC_STROBE_PAAR_EINDEUTIG=0`).
            if unique[1 - r] && unique_mode() {
                continue;
            }
            let best = out[r][0].score;
            for a in out[r]
                .iter()
                .filter(|a| a.score >= best - opt.pair.pen_unpaired)
            {
                let (ca, ab, ae, ar) = at(a);
                let paired = out[1 - r].iter().any(|b| {
                    let (cb, bb, be, br) = at(b);
                    if cb != ca || br == ar {
                        return false;
                    }
                    let (f, rv) = if ar {
                        ((bb, be), (ab, ae))
                    } else {
                        ((ab, ae), (bb, be))
                    };
                    let span = rv.1 - f.0;
                    (0..=max).contains(&span)
                });
                if !paired {
                    if unsure_mode() {
                        UNSURE.fetch_add(1, Relaxed);
                        return Some((out, Kind::Unsure, frac));
                    }
                    // Mit Partnersuche, aber ohne vollen Weg (`GC_STROBE_PAAR_SUCHE=1`).
                    // Mit gesperrten Randstrobes könnten Kopien verborgen sein: dann voller Weg.
                    if rescue_mode() && (frac[0] == 0.0 && frac[1] == 0.0 || !rescue_needs_clean())
                    {
                        return Some((out, Kind::Rescue, frac));
                    }
                    return None;
                }
            }
        }
    }
    if crowded {
        UNSURE.fetch_add(1, Relaxed);
        return Some((out, Kind::Unsure, frac));
    }
    TAKEN.fetch_add(1, Relaxed);
    Some((out, Kind::Sure, frac))
}

/// Richtet den Read am Anker `g` aus: ohne Lücke, wenn höchstens [`max_mismatch`] Basen
/// abweichen, sonst lokal; abgeschnittene Enden werden ohne Lücke angehängt, wenn das je Seite
/// höchstens den Abzug fürs Abschneiden kostet.
pub(crate) fn align_anchor(
    index: &FmIndex,
    read: &[u8],
    g: u64,
    opt: &AlignOptions,
) -> Option<Region> {
    let sc = &opt.region.scoring;
    let so = &opt.strobe;
    let len = read.len();
    let l = index.reference().len();
    let (strand, diag) = split(g);
    let fb = (diag - so.margin as i64).max(0) as u64;
    let fe = ((diag + len as i64 + so.margin as i64).max(0) as u64).min(l);
    if fe <= fb + len as u64 / 2 {
        return None;
    }
    let (tb, te) = if strand == 0 {
        (fb, fe)
    } else {
        (2 * l - fe, 2 * l - fb)
    };
    let target = index.text_range(tb, te);
    if target.iter().any(|&c| c > 3) {
        return None;
    }
    let off = if strand == 0 {
        diag - fb as i64
    } else {
        fe as i64 - diag - len as i64
    };
    let full_len = |rb: u64, score: i32| -> Option<Region> {
        let hit = index.hit_at(rb, len as u64)?;
        Some(region(
            rb,
            rb + len as u64,
            0,
            len as i32,
            hit.contig,
            score,
            score,
            sc,
        ))
    };
    if off >= 0 && off as usize + len <= target.len() {
        let t = &target[off as usize..off as usize + len];
        let mm = read.iter().zip(t).filter(|(a, b)| a != b).count();
        if mm <= max_mismatch() {
            let score = (len - mm) as i32 * sc.matching - mm as i32 * sc.mismatch;
            return full_len(tb + off as u64, score);
        }
    }
    let byte = len as i32 * sc.matching < 250;
    let r = crate::local::local_align(read, &target, sc, byte, 0, true);
    if r.score < opt.sam.min_score {
        return None;
    }
    let (mut rb, mut re, mut qb, mut qe) = (tb + r.tb as u64, tb + r.te as u64 + 1, r.qb, r.qe + 1);
    let mut truesc = r.score;
    let side = |q: &[u8], t: &[u8]| -> i32 {
        q.iter()
            .zip(t)
            .map(|(a, b)| if a == b { sc.matching } else { -sc.mismatch })
            .sum()
    };
    if qb > 0 {
        let n = qb as u64;
        if rb >= n && (rb < l) == (rb - n < l) {
            let d = side(&read[..qb as usize], &index.text_range(rb - n, rb));
            if d >= -opt.region.pen_clip5 {
                truesc += d;
                rb -= n;
                qb = 0;
            }
        }
    }
    if qe < len as i32 {
        let n = (len as i32 - qe) as u64;
        if re + n <= 2 * l && (re <= l) == (re + n <= l) {
            let d = side(&read[qe as usize..], &index.text_range(re, re + n));
            if d >= -opt.region.pen_clip3 {
                truesc += d;
                re += n;
                qe = len as i32;
            }
        }
    }
    let hit = index.hit_at(rb, re - rb)?;
    Some(region(rb, re, qb, qe, hit.contig, r.score, truesc, sc))
}

#[allow(clippy::too_many_arguments)]
fn region(
    rb: u64,
    re: u64,
    qb: i32,
    qe: i32,
    contig: u32,
    score: i32,
    truesc: i32,
    sc: &crate::extend::Scoring,
) -> Region {
    Region {
        rb,
        re,
        qb,
        qe,
        contig,
        score,
        truesc,
        w: sc.band_width,
        seedcov: qe - qb,
        seedlen0: qe - qb,
        sub: 0,
        sub_n: 0,
        secondary: -1,
        csub: 0,
        rescued: false,
    }
}
