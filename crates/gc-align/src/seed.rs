//! Seeds eines Reads sammeln und im Genom orten (BWA-MEM: `mem_collect_intv`, `mem_chain`).

use gc_index::{FmIndex, ForwardSeed, Smem, SmemSearch};

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
    /// Die Seed-Suchen mehrerer Reads verschränkt rechnen ([`collect_intervals_many`]). Gleiche
    /// Treffer; aus heißt Read für Read wie bisher.
    pub interleave: bool,
    /// Erst alle Reads eines Stapels über den Schnellweg, dann die übrigen zu vollen Gruppen
    /// für die verschränkte Suche zusammenfassen. Aus: Gruppen nach der Reihenfolge der Reads,
    /// in denen nach dem Schnellweg oft nur wenige Reads übrig sind. Gleiche Treffer.
    pub regroup: bool,
    /// Der Schnellweg sucht seine Stücke gemeinsam ([`gc_index::FmIndex::intervals_many`]):
    /// Tabelle der kurzen Muster und überlappende Speicherzugriffe. Gleiche Bereiche; aus sucht
    /// jedes Stück einzeln von vorn.
    pub fast_pieces_together: bool,
    /// Beim Orten jede Zeile des Index je Read nur einmal nachschlagen. Gleiche Orte; aus ortet
    /// jede Zeile so oft, wie sie in den Treffern vorkommt.
    pub locate_unique: bool,
}

impl Default for SeedOptions {
    fn default() -> Self {
        Self {
            min_seed_len: 19,
            split_factor: 1.5,
            split_width: 10,
            max_occ: 500,
            max_mem_intv: 20,
            interleave: true,
            regroup: true,
            fast_pieces_together: true,
            locate_unique: true,
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

/// [`collect_intervals`] für mehrere Reads zugleich, mit genau denselben Ergebnissen. Die
/// Suchen der Reads laufen abwechselnd Schritt für Schritt: Erst fordert jede die Speicherzeilen
/// ihres nächsten Schritts an, dann rechnet jede ihren Schritt. Der Prozessor wartet so auf die
/// Zugriffe vieler Reads zugleich statt auf einen nach dem anderen (wie BWA-MEM2).
pub fn collect_intervals_many(
    index: &FmIndex,
    reads: &[&[u8]],
    opt: &SeedOptions,
) -> Vec<Vec<Smem>> {
    let mut drivers: Vec<Driver<'_>> = reads.iter().map(|&r| Driver::new(r)).collect();
    let mut active: Vec<usize> = (0..drivers.len()).collect();
    while !active.is_empty() {
        for &d in &active {
            drivers[d].prefetch(index);
        }
        active.retain(|&d| !drivers[d].advance(index, opt));
    }
    drivers.into_iter().map(|d| d.all).collect()
}

/// Ablauf der drei Runden von [`collect_intervals`] für einen Read, in Schritten.
struct Driver<'r> {
    read: &'r [u8],
    stage: Stage,
    job: Job,
    all: Vec<Smem>,
    extra: Vec<Smem>,
}

#[derive(Clone, Copy)]
enum Stage {
    /// Erste Runde, nächste Startposition.
    First(usize),
    /// Nachsuche: nächster Treffer der ersten Runde, Zahl der Treffer der ersten Runde.
    Split(usize, usize),
    /// Dritte Runde, nächste Startposition.
    Rare(usize),
    Done,
}

enum Job {
    Idle,
    Smem(SmemSearch),
    Seed(ForwardSeed),
}

impl<'r> Driver<'r> {
    fn new(read: &'r [u8]) -> Self {
        Driver {
            read,
            stage: Stage::First(0),
            job: Job::Idle,
            all: Vec::new(),
            extra: Vec::new(),
        }
    }

    fn prefetch(&self, index: &FmIndex) {
        match &self.job {
            Job::Smem(s) => s.prefetch(index, self.read),
            Job::Seed(s) => s.prefetch(index, self.read),
            Job::Idle => {}
        }
    }

    /// Ein Schritt. `true`, wenn der Read fertig ist.
    fn advance(&mut self, index: &FmIndex, opt: &SeedOptions) -> bool {
        let read = self.read;
        match &mut self.job {
            Job::Smem(search) => {
                let out = match self.stage {
                    Stage::First(_) => &mut self.all,
                    _ => &mut self.extra,
                };
                let Some(reach) = search.step(index, read, out) else {
                    return false;
                };
                self.job = Job::Idle;
                match self.stage {
                    Stage::First(_) => {
                        self.stage = Stage::First(reach);
                    }
                    Stage::Split(i, first_round) => {
                        let min = opt.min_seed_len;
                        self.all
                            .extend(self.extra.iter().filter(|e| e.len() >= min));
                        self.stage = Stage::Split(i + 1, first_round);
                    }
                    _ => unreachable!("SMEM-Suche nur in den ersten beiden Runden"),
                }
            }
            Job::Seed(search) => {
                let Some((next, found)) = search.step(index, read) else {
                    return false;
                };
                self.job = Job::Idle;
                self.all.extend(found);
                self.stage = Stage::Rare(next);
            }
            Job::Idle => {}
        }
        self.start_next(index, opt)
    }

    /// Legt die nächste Suche an, ohne sie zu rechnen. `true`, wenn keine mehr bleibt.
    fn start_next(&mut self, index: &FmIndex, opt: &SeedOptions) -> bool {
        let read = self.read;
        loop {
            match self.stage {
                Stage::First(x) if x < read.len() => {
                    self.job = Job::Smem(SmemSearch::new(index, read, x, 1, &self.all));
                    return false;
                }
                Stage::First(_) => {
                    self.all.retain(|m| m.len() >= opt.min_seed_len);
                    self.stage = Stage::Split(0, self.all.len());
                }
                Stage::Split(i, first_round) if i < first_round => {
                    let split_len = (opt.min_seed_len as f32 * opt.split_factor + 0.499) as u32;
                    let m = self.all[i];
                    if m.len() < split_len || m.interval.size > opt.split_width {
                        self.stage = Stage::Split(i + 1, first_round);
                        continue;
                    }
                    let mid = (m.start + m.end) as usize / 2;
                    self.extra.clear();
                    self.job = Job::Smem(SmemSearch::new(
                        index,
                        read,
                        mid,
                        m.interval.size + 1,
                        &self.extra,
                    ));
                    return false;
                }
                Stage::Split(..) => {
                    self.stage = if opt.max_mem_intv > 0 {
                        Stage::Rare(0)
                    } else {
                        Stage::Done
                    };
                }
                Stage::Rare(x) if x < read.len() => {
                    self.job = Job::Seed(ForwardSeed::new(
                        index,
                        read,
                        x,
                        opt.min_seed_len,
                        opt.max_mem_intv,
                    ));
                    return false;
                }
                Stage::Rare(_) => self.stage = Stage::Done,
                Stage::Done => {
                    self.all.sort_unstable_by_key(|m| (m.start, m.end));
                    return true;
                }
            }
        }
    }
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
    if opt.locate_unique {
        // Dieselbe Zeile kommt oft aus mehreren Treffern (derselbe Treffer aus mehreren Runden,
        // ein Treffer als Anfang eines längeren). Jede Zeile nur einmal orten.
        let mut unique = rows.clone();
        unique.sort_unstable();
        unique.dedup();
        let mut found = Vec::with_capacity(unique.len());
        index.locate_many(&unique, &mut found);
        positions.extend(rows.iter().map(|r| found[unique.binary_search(r).unwrap()]));
    } else {
        index.locate_many(&rows, &mut positions);
    }
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
    fn many_reads_at_once_give_the_same_intervals() {
        let mut rng = 0x2545_F491_4F6C_DD1Du64;
        let mut s: Vec<u8> = (0..40_000)
            .map(|_| b"ACGT"[(xorshift(&mut rng) & 3) as usize])
            .collect();
        // Wiederholungen und arme Folgen, damit alle drei Runden etwas zu tun haben.
        for k in 0..30 {
            let from = (xorshift(&mut rng) % 39_000) as usize;
            let to = (xorshift(&mut rng) % 39_000) as usize;
            let len = 20 + (xorshift(&mut rng) % 200) as usize;
            let len = len.min(40_000 - from.max(to));
            let piece = s[from..from + len].to_vec();
            s[to..to + len].copy_from_slice(&piece);
            if k % 5 == 0 {
                for b in &mut s[to..to + len / 3] {
                    *b = b'A';
                }
            }
        }
        let mut r = Reference::new();
        r.push_contig("a", &s[..25_000]);
        r.push_contig("b", &s[25_000..]);
        let code = |b: u8| gc_index::reference::encode(b).unwrap();
        for variant in [Variant::Standard, Variant::Kompakt] {
            let index = FmIndex::build(r.clone(), variant, 0).unwrap().0;
            let reads: Vec<Vec<u8>> = (0..400)
                .map(|n| {
                    let len = 30 + (xorshift(&mut rng) % 221) as usize;
                    let at = (xorshift(&mut rng) % (40_000 - len as u64)) as usize;
                    let mut read: Vec<u8> = s[at..at + len].iter().map(|&b| code(b)).collect();
                    if n % 3 == 1 {
                        read.reverse();
                        for b in &mut read {
                            *b = 3 - *b;
                        }
                    }
                    for b in &mut read {
                        match xorshift(&mut rng) % 60 {
                            0 => *b = (*b + 1) % 4,
                            1 => *b = 4,
                            _ => {}
                        }
                    }
                    if n % 7 == 0 {
                        read = (0..len).map(|_| (xorshift(&mut rng) & 3) as u8).collect();
                    }
                    read
                })
                .collect();
            let views: Vec<&[u8]> = reads.iter().map(|r| &r[..]).collect();
            for opt in [
                SeedOptions::default(),
                SeedOptions {
                    max_mem_intv: 0,
                    ..SeedOptions::default()
                },
            ] {
                let many = collect_intervals_many(&index, &views, &opt);
                let plain = SeedOptions {
                    locate_unique: false,
                    ..opt.clone()
                };
                for (read, got) in reads.iter().zip(&many) {
                    assert_eq!(got, &collect_intervals(&index, read, &opt));
                    // Jede Zeile einmal orten gibt dieselben Seeds.
                    assert_eq!(
                        locate_seeds(&index, got, &opt),
                        locate_seeds(&index, got, &plain)
                    );
                }
            }
        }
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
