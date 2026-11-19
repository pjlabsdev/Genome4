//! Seeds eines Reads sammeln und im Genom orten (BWA-MEM: `mem_collect_intv`, `mem_chain`).

use gc_index::{FmIndex, Smem};

/// Einstellungen der Seed-Suche. Die Vorgaben sind die von BWA-MEM.
#[derive(Debug, Clone, PartialEq)]
pub struct SeedOptions {
    /// Kürzester Seed (`-k`).
    pub min_seed_len: u32,
    /// Treffer ab `min_seed_len * split_factor` Basen werden in der Mitte neu angesetzt (`-r`).
    pub split_factor: f32,
    /// Nachsuche nur bei Treffern mit höchstens so vielen Fundstellen (`-y` in bwa: split_width).
    pub split_width: u64,
    /// Je Treffer höchstens so viele Fundstellen orten (`-c`); häufigere werden ausgedünnt.
    pub max_occ: u64,
    /// Dritte Runde: Seeds, die seltener als so oft vorkommen, von links nach rechts ohne
    /// Linksverlängerung suchen (`-y`, bwa: max_mem_intv). 0 schaltet die Runde ab.
    pub max_mem_intv: u64,
}

impl Default for SeedOptions {
    fn default() -> Self {
        Self {
            min_seed_len: 19,
            split_factor: 1.5,
            split_width: 10,
            max_occ: 500,
            max_mem_intv: 20,
        }
    }
}

/// Ein georteter Seed: `read[qbeg..qbeg + len]` gleicht dem Indextext ab `rbeg`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seed {
    /// Position im Indextext: unter L auf dem Vorwärtsstrang, ab L im reversen Komplement.
    pub rbeg: u64,
    pub qbeg: u32,
    pub len: u32,
    /// Contig, in dem der Seed liegt.
    pub contig: u32,
    /// Rohe Fundstellenzahl des Treffers, aus dem der Seed stammt.
    pub occurrences: u64,
}

/// Treffer aller drei Runden (SMEM, Nachsuche in langen Treffern, seltene Seeds von links nach
/// rechts), nach Anfang und Ende geordnet. Doppelte bleiben erhalten.
pub fn collect_intervals(index: &FmIndex, read: &[u8], opt: &SeedOptions) -> Vec<Smem> {
    let mut all = index.smems(read, opt.min_seed_len);
    // Nachsuche: Ein langer, seltener Treffer kann kürzere, häufigere Treffer verdecken, die zum
    // richtigen Ort führen. Deshalb in seiner Mitte neu ansetzen und nur Treffer zulassen, die
    // häufiger vorkommen als er selbst.
    let split_len = (opt.min_seed_len as f32 * opt.split_factor + 0.499) as u32;
    let first_round = all.len();
    let mut extra = Vec::new();
    for i in 0..first_round {
        let m = all[i];
        if m.len() < split_len || m.interval.size > opt.split_width {
            continue;
        }
        let mid = (m.start + m.end) as usize / 2;
        extra.clear();
        index.smems_at(read, mid, m.interval.size + 1, &mut extra);
        all.extend(extra.iter().filter(|e| e.len() >= opt.min_seed_len));
    }
    // Dritte Runde: seltene Seeds, die ein längerer, anderswo passender Treffer verdeckt.
    if opt.max_mem_intv > 0 {
        let mut x = 0;
        while x < read.len() {
            let (next, found) = index.seed_forward(read, x, opt.min_seed_len, opt.max_mem_intv);
            all.extend(found);
            x = next;
        }
    }
    // Doppelte bleiben stehen, wie in bwa: Derselbe Treffer kann aus mehreren Runden kommen, und
    // jeder ergibt eigene Seeds. In Wiederholungen (etwa einer langen Folge gleicher Basen)
    // entscheidet das mit, welche Ketten entstehen.
    all.sort_unstable_by_key(|m| (m.start, m.end));
    all
}

/// Anteil des Reads, der von Treffern mit mehr als `max_occ` Fundstellen überdeckt ist
/// (bwa: `frac_rep`). Die Treffer müssen nach Anfang geordnet sein.
pub fn repeat_fraction(intervals: &[Smem], max_occ: u64, read_len: usize) -> f32 {
    let (mut b, mut e, mut l_rep) = (0u32, 0u32, 0u32);
    for m in intervals.iter().filter(|m| m.interval.size > max_occ) {
        if m.start > e {
            l_rep += e - b;
            b = m.start;
            e = m.end;
        } else {
            e = e.max(m.end);
        }
    }
    l_rep += e - b;
    if read_len == 0 {
        0.0
    } else {
        l_rep as f32 / read_len as f32
    }
}

/// Ortet die Fundstellen der Treffer, siehe [`collect_seeds`]. Erst werden die Zeilen aller
/// Treffer des Reads gesammelt, dann gemeinsam geortet (`FmIndex::locate_many`): Die Wege durch
/// den Index laufen so im Gleichschritt statt nacheinander.
pub fn locate_seeds(index: &FmIndex, intervals: &[Smem], opt: &SeedOptions) -> Vec<Seed> {
    let (mut rows, mut owner) = (Vec::new(), Vec::new());
    for (i, m) in intervals.iter().enumerate() {
        let iv = m.interval;
        let step = if iv.size > opt.max_occ {
            iv.size / opt.max_occ
        } else {
            1
        };
        let mut k = 0;
        let mut count = 0;
        while k < iv.size && count < opt.max_occ {
            rows.push(iv.k + k);
            owner.push(i);
            k += step;
            count += 1;
        }
    }
    let mut positions = Vec::with_capacity(rows.len());
    index.locate_many(&rows, &mut positions);
    let mut seeds = Vec::with_capacity(rows.len());
    for (&rbeg, &i) in positions.iter().zip(&owner) {
        let m = &intervals[i];
        if let Some(hit) = index.hit_at(rbeg, u64::from(m.len())) {
            seeds.push(Seed {
                rbeg,
                qbeg: m.start,
                len: m.len(),
                contig: hit.contig,
                occurrences: m.interval.size,
            });
        }
    }
    seeds
}

/// Sammelt die Seeds eines Reads (Codes 0..3, sonst N) und ortet sie. Die Reihenfolge ist die
/// der Treffer (Anfang, Ende im Read), darin nach Zeile im Index. Fundstellen, die über eine
/// Contig- oder Stranggrenze reichen oder eine N-Strecke berühren, fallen weg.
pub fn collect_seeds(index: &FmIndex, read: &[u8], opt: &SeedOptions) -> Vec<Seed> {
    locate_seeds(index, &collect_intervals(index, read, opt), opt)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gc_index::{Reference, Variant};

    fn xorshift(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn seeds_match_the_text_and_reseeding_finds_repeats() {
        let mut rng = 0x9E37_79B9_7F4A_7C15u64;
        let mut s: Vec<u8> = (0..3000)
            .map(|_| b"ACGT"[(xorshift(&mut rng) & 3) as usize])
            .collect();
        // Drei Kopien eines 40-Basen-Stücks; die erste steht in einem längeren eindeutigen Umfeld.
        let copy = s[500..540].to_vec();
        s[1200..1240].copy_from_slice(&copy);
        s[2500..2540].copy_from_slice(&copy);
        // Die Kopien sollen sich nicht zufällig über ihre Ränder hinaus fortsetzen.
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        for at in [1199, 2499] {
            s[at] = other(s[499]);
        }
        for at in [1240, 2540] {
            s[at] = other(s[540]);
        }
        let mut r = Reference::new();
        r.push_contig("a", &s[..1500]);
        r.push_contig("b", &s[1500..]);
        let index = FmIndex::build(r, Variant::Kompakt, 0).unwrap().0;
        let text = index.reference().index_text();
        let code = |b: u8| gc_index::reference::encode(b).unwrap();

        // Read über die erste Kopie samt eindeutigem Umfeld.
        let read: Vec<u8> = s[470..570].iter().map(|&b| code(b)).collect();
        // Zuerst ohne die dritte Runde.
        let opt = SeedOptions {
            max_mem_intv: 0,
            ..SeedOptions::default()
        };
        let intervals = collect_intervals(&index, &read, &opt);
        // Erste Suche: ein Treffer über den ganzen Read. Nachsuche: die Kopie mit drei Fundstellen.
        assert!(
            intervals
                .iter()
                .any(|m| m.len() == 100 && m.interval.size == 1)
        );
        assert!(
            intervals
                .iter()
                .any(|m| m.start == 30 && m.end == 70 && m.interval.size == 3)
        );
        let seeds = collect_seeds(&index, &read, &opt);
        assert_eq!(seeds.len(), 4);
        for sd in &seeds {
            for i in 0..sd.len as usize {
                assert_eq!(
                    text[sd.rbeg as usize + i] - 1,
                    read[sd.qbeg as usize + i],
                    "{sd:?}"
                );
            }
        }
        assert_eq!(seeds.iter().filter(|sd| sd.contig == 1).count(), 1);

        // Mit der dritten Runde kommen kurze Seeds dazu; auch sie müssen zum Text passen.
        let with_third = collect_seeds(&index, &read, &SeedOptions::default());
        assert!(with_third.len() > seeds.len());
        for sd in &with_third {
            for i in 0..sd.len as usize {
                assert_eq!(text[sd.rbeg as usize + i] - 1, read[sd.qbeg as usize + i]);
            }
        }

        // Ohne Nachsuche (Schwelle unerreichbar) bleibt nur der lange Treffer.
        let no_split = SeedOptions {
            split_factor: 100.0,
            ..opt.clone()
        };
        assert_eq!(collect_seeds(&index, &read, &no_split).len(), 1);

        // Ausdünnen häufiger Treffer.
        let thin = SeedOptions {
            max_occ: 2,
            ..opt.clone()
        };
        let short: Vec<u8> = copy.iter().map(|&b| code(b)).collect();
        assert_eq!(collect_seeds(&index, &short, &thin).len(), 2);
    }

    /// In einer Folge gleicher Basen liefern mehrere Runden denselben Treffer. bwa behält die
    /// Doppelten, und sie entscheiden mit, welche Ketten entstehen. Deshalb bleiben sie stehen.
    #[test]
    fn repeated_intervals_are_kept() {
        let mut rng = 0x1234_5678_9ABC_DEF1u64;
        let mut s: Vec<u8> = (0..4000)
            .map(|_| b"CGT"[(xorshift(&mut rng) % 3) as usize])
            .collect();
        for at in [700, 1900, 3100] {
            s[at..at + 60].fill(b'A');
        }
        let mut r = Reference::new();
        r.push_contig("a", &s);
        let index = FmIndex::build(r, Variant::Kompakt, 0).unwrap().0;
        let code = |b: u8| gc_index::reference::encode(b).unwrap();
        let mut read: Vec<u8> = s[660..700].iter().map(|&b| code(b)).collect();
        read.extend(std::iter::repeat_n(code(b'A'), 75));
        read.extend(s[760..790].iter().map(|&b| code(b)));
        let intervals = collect_intervals(&index, &read, &SeedOptions::default());
        assert!(intervals.is_sorted_by_key(|m| (m.start, m.end)));
        let mut distinct: Vec<_> = intervals.iter().map(|m| (m.start, m.end)).collect();
        distinct.dedup();
        assert!(
            distinct.len() < intervals.len(),
            "kein doppelter Treffer unter {} Treffern",
            intervals.len()
        );
    }
}
