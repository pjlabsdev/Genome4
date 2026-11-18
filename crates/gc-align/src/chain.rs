//! Seeds zu Ketten verbinden (BWA-MEM: `mem_chain`, `test_and_merge`, `mem_chain_weight`).
//!
//! Eine Kette ist eine Folge von Seeds, die auf demselben Strang desselben Contigs in gleicher
//! Reihenfolge und ungefähr auf derselben Diagonale liegen. Sie beschreibt einen Kandidaten für
//! den Ort des Reads; das Alignment füllt später die Lücken zwischen den Seeds.

use crate::kbtree::KbTree;
use crate::seed::Seed;

/// Einstellungen des Verkettens. Die Vorgaben sind die von BWA-MEM.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainOptions {
    /// Bandbreite: so weit dürfen zwei Seeds einer Kette von der gemeinsamen Diagonale abweichen (`-w`).
    pub band_width: i64,
    /// Größte Lücke zwischen zwei Seeds einer Kette (`max_chain_gap`).
    pub max_chain_gap: i64,
}

impl Default for ChainOptions {
    fn default() -> Self {
        Self {
            band_width: 100,
            max_chain_gap: 10_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chain {
    pub contig: u32,
    /// Seeds in der Reihenfolge, in der sie angehängt wurden; Positionen im Read und in der
    /// Referenz steigen.
    pub seeds: Vec<Seed>,
}

impl Chain {
    /// Position des ersten Seeds im Indextext.
    pub fn pos(&self) -> u64 {
        self.seeds[0].rbeg
    }

    /// Die Kette liegt im reversen Komplement (`l` = Länge der Referenz).
    pub fn is_reverse(&self, l: u64) -> bool {
        self.seeds[0].rbeg >= l
    }

    /// Gewicht: Zahl der Basen, die die Seeds überdecken, ohne Überlappungen doppelt zu zählen,
    /// im Read und in der Referenz getrennt gerechnet; das kleinere von beiden.
    pub fn weight(&self) -> u32 {
        let cover = |range: &dyn Fn(&Seed) -> (u64, u64)| -> u64 {
            let (mut total, mut end) = (0, 0);
            for s in &self.seeds {
                let (b, e) = range(s);
                if b >= end {
                    total += e - b;
                } else if e > end {
                    total += e - end;
                }
                end = end.max(e);
            }
            total
        };
        let q = cover(&|s| (u64::from(s.qbeg), u64::from(s.qbeg + s.len)));
        let r = cover(&|s| (s.rbeg, s.rbeg + u64::from(s.len)));
        q.min(r) as u32
    }

    /// Versucht, den Seed anzuhängen. `true`, wenn er zur Kette gehört (angehängt oder bereits
    /// in ihr enthalten).
    fn test_and_merge(&mut self, p: &Seed, l: u64, opt: &ChainOptions) -> bool {
        let first = self.seeds[0];
        let last = *self.seeds.last().expect("Kette hat mindestens einen Seed");
        if p.contig != self.contig {
            return false;
        }
        let qend = last.qbeg + last.len;
        let rend = last.rbeg + u64::from(last.len);
        if p.qbeg >= first.qbeg
            && p.qbeg + p.len <= qend
            && p.rbeg >= first.rbeg
            && p.rbeg + u64::from(p.len) <= rend
        {
            return true;
        }
        // Nie über die Grenze zwischen den Strängen verketten.
        if (last.rbeg < l || first.rbeg < l) && p.rbeg >= l {
            return false;
        }
        let x = i64::from(p.qbeg) - i64::from(last.qbeg);
        let y = p.rbeg as i64 - last.rbeg as i64;
        if y >= 0
            && x - y <= opt.band_width
            && y - x <= opt.band_width
            && x - i64::from(last.len) < opt.max_chain_gap
            && y - i64::from(last.len) < opt.max_chain_gap
        {
            self.seeds.push(*p);
            return true;
        }
        false
    }
}

/// Verbindet Seeds zu Ketten. `l` ist die Länge der Referenz (ein Strang). Jeder Seed wird mit
/// der Kette verglichen, deren Anfang in der Referenz am nächsten vor ihm liegt; passt er nicht,
/// beginnt er eine neue Kette. Ergebnis nach Position im Indextext geordnet.
pub fn chain_seeds(seeds: &[Seed], l: u64, opt: &ChainOptions) -> Vec<Chain> {
    // Schlüssel: Position des ersten Seeds. Gleiche Schlüssel sind häufig; der Baum verhält
    // sich dabei wie der von bwa (siehe `kbtree`).
    let mut tree: KbTree<Chain> = KbTree::new();
    for p in seeds {
        let merged = tree
            .lower_mut(p.rbeg)
            .is_some_and(|c| c.test_and_merge(p, l, opt));
        if !merged {
            tree.insert(
                p.rbeg,
                Chain {
                    contig: p.contig,
                    seeds: vec![*p],
                },
            );
        }
    }
    tree.into_values()
}

/// Einstellungen des Kettenfilters. Die Vorgaben sind die von BWA-MEM.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterOptions {
    /// Ketten unter diesem Gewicht entfallen (`-W`).
    pub min_chain_weight: u32,
    /// Zwei Ketten überlappen stark, wenn sie diesen Anteil der kürzeren im Read teilen.
    pub mask_level: f32,
    /// Eine stark überlappende Kette entfällt, wenn ihr Gewicht unter diesem Anteil der
    /// schwereren liegt (`-D`) ...
    pub drop_ratio: f32,
    /// ... und der Abstand der Gewichte mindestens das Doppelte hiervon beträgt (`-k`).
    pub min_seed_len: u32,
    pub max_chain_gap: i64,
}

impl Default for FilterOptions {
    fn default() -> Self {
        Self {
            min_chain_weight: 0,
            mask_level: 0.5,
            drop_ratio: 0.5,
            min_seed_len: 19,
            max_chain_gap: 10_000,
        }
    }
}

/// Filtert Ketten (BWA-MEM: `mem_chain_flt`). Ergebnis nach fallendem Gewicht geordnet.
/// Eine Kette entfällt, wenn sie im Read stark mit einer schwereren überlappt und deutlich
/// leichter ist. Von den Ketten, die eine behaltene stark überlappen, bleibt je behaltener Kette
/// die erste trotzdem erhalten: Sie wird gebraucht, um später die Eindeutigkeit zu beurteilen.
pub fn filter_chains(chains: Vec<Chain>, opt: &FilterOptions) -> Vec<Chain> {
    let mut a: Vec<(u32, Chain)> = chains
        .into_iter()
        .map(|c| (c.weight(), c))
        .filter(|(w, _)| *w >= opt.min_chain_weight)
        .collect();
    if a.is_empty() {
        return Vec::new();
    }
    // Reihenfolge bei gleichem Gewicht wie in bwa, sonst bleibt eine andere Kette "die erste".
    crate::ksort::introsort(&mut a, |x, y| x.0 > y.0);
    let beg = |c: &Chain| i64::from(c.seeds[0].qbeg);
    let end = |c: &Chain| {
        let s = c.seeds[c.seeds.len() - 1];
        i64::from(s.qbeg + s.len)
    };
    // 0 = entfällt, 1 = nur als erster Überlapper behalten, 2 = behalten mit Überlappung, 3 = frei
    let mut kept = vec![0u8; a.len()];
    let mut first: Vec<Option<usize>> = vec![None; a.len()];
    let mut chosen = vec![0usize];
    kept[0] = 3;
    for i in 1..a.len() {
        let mut large_overlap = false;
        let mut dropped = false;
        for &j in &chosen {
            let b_max = beg(&a[j].1).max(beg(&a[i].1));
            let e_min = end(&a[j].1).min(end(&a[i].1));
            if e_min > b_max {
                let li = end(&a[i].1) - beg(&a[i].1);
                let lj = end(&a[j].1) - beg(&a[j].1);
                let min_l = li.min(lj);
                if (e_min - b_max) as f32 >= min_l as f32 * opt.mask_level
                    && min_l < opt.max_chain_gap
                {
                    large_overlap = true;
                    if first[j].is_none() {
                        first[j] = Some(i);
                    }
                    if (a[i].0 as f32) < a[j].0 as f32 * opt.drop_ratio
                        && a[j].0 - a[i].0 >= opt.min_seed_len << 1
                    {
                        dropped = true;
                        break;
                    }
                }
            }
        }
        if !dropped {
            chosen.push(i);
            kept[i] = if large_overlap { 2 } else { 3 };
        }
    }
    for &j in &chosen {
        if let Some(i) = first[j]
            && kept[i] == 0
        {
            kept[i] = 1;
        }
    }
    a.into_iter()
        .zip(kept)
        .filter(|(_, k)| *k != 0)
        .map(|((_, c), _)| c)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seed(rbeg: u64, qbeg: u32, len: u32) -> Seed {
        Seed {
            rbeg,
            qbeg,
            len,
            contig: 0,
            occurrences: 1,
        }
    }

    #[test]
    fn seeds_on_one_diagonal_form_one_chain() {
        let seeds = [seed(1000, 0, 30), seed(1031, 31, 40), seed(1080, 78, 25)];
        let chains = chain_seeds(&seeds, 1_000_000, &ChainOptions::default());
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].seeds.len(), 3);
        // Read: 0..30, 31..71, 78..103 = 95; Referenz: 1000..1030, 1031..1071, 1080..1105 = 95
        assert_eq!(chains[0].weight(), 95);
    }

    #[test]
    fn filter_drops_light_overlapping_chains_but_keeps_the_first() {
        let chain = |rbeg: u64, qbeg: u32, len: u32| Chain {
            contig: 0,
            seeds: vec![seed(rbeg, qbeg, len)],
        };
        let chains = vec![
            chain(5000, 10, 20), // leicht, überlappt die schwere ganz: erster Überlapper, bleibt
            chain(1000, 0, 100), // schwer
            chain(7000, 20, 25), // leicht, überlappt: entfällt
            chain(9000, 120, 22), // leicht, daneben: bleibt
            chain(3000, 0, 80),  // fast so schwer, überlappt: bleibt
        ];
        let kept = filter_chains(chains, &FilterOptions::default());
        let pos: Vec<u64> = kept.iter().map(Chain::pos).collect();
        assert_eq!(pos, vec![1000, 3000, 9000]);
        // Ohne die fast so schwere Kette: Die Kette mit Gewicht 25 ist der erste Überlapper der
        // schweren und bleibt deshalb, obwohl sie deutlich leichter ist; die mit 20 entfällt.
        let chains = vec![
            chain(5000, 10, 20),
            chain(1000, 0, 100),
            chain(7000, 20, 25),
        ];
        let kept = filter_chains(chains, &FilterOptions::default());
        let pos: Vec<u64> = kept.iter().map(Chain::pos).collect();
        assert_eq!(pos, vec![1000, 7000]);
        assert!(filter_chains(Vec::new(), &FilterOptions::default()).is_empty());
    }

    #[test]
    fn overlapping_seeds_are_counted_once() {
        let c = Chain {
            contig: 0,
            seeds: vec![seed(100, 0, 30), seed(120, 20, 30)],
        };
        assert_eq!(c.weight(), 50);
    }

    #[test]
    fn far_seeds_strands_and_contigs_stay_apart() {
        let l = 10_000;
        let opt = ChainOptions::default();
        // weit weg von der Diagonale
        let chains = chain_seeds(&[seed(1000, 0, 30), seed(1500, 40, 30)], l, &opt);
        assert_eq!(chains.len(), 2);
        // anderer Strang
        let chains = chain_seeds(&[seed(9990, 0, 20), seed(10_020, 30, 20)], l, &opt);
        assert_eq!(chains.len(), 2);
        // anderer Contig
        let mut other = seed(1040, 40, 30);
        other.contig = 1;
        let chains = chain_seeds(&[seed(1000, 0, 30), other], l, &opt);
        assert_eq!(chains.len(), 2);
        // enthaltener Seed wird nicht doppelt aufgenommen
        let chains = chain_seeds(&[seed(1000, 0, 50), seed(1010, 10, 20)], l, &opt);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].seeds.len(), 1);
        // Referenz läuft rückwärts: neue Kette
        let chains = chain_seeds(&[seed(1050, 0, 30), seed(1000, 40, 30)], l, &opt);
        assert_eq!(chains.len(), 2);
    }
}
