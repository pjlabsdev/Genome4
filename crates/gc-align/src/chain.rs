//! Seeds zu Ketten verbinden (BWA-MEM: `mem_chain`, `test_and_merge`, `mem_chain_weight`).
//!
//! Eine Kette ist eine Folge von Seeds, die auf demselben Strang desselben Contigs in gleicher
//! Reihenfolge und ungefähr auf derselben Diagonale liegen. Sie beschreibt einen Kandidaten für
//! den Ort des Reads; das Alignment füllt später die Lücken zwischen den Seeds.

use crate::seed::Seed;
use std::collections::BTreeMap;

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
    // Schlüssel: Position des ersten Seeds, dazu eine laufende Nummer für gleiche Positionen.
    let mut tree: BTreeMap<(u64, u32), Chain> = BTreeMap::new();
    for (n, p) in seeds.iter().enumerate() {
        let merged = tree
            .range_mut(..=(p.rbeg, u32::MAX))
            .next_back()
            .is_some_and(|(_, c)| c.test_and_merge(p, l, opt));
        if !merged {
            tree.insert(
                (p.rbeg, n as u32),
                Chain {
                    contig: p.contig,
                    seeds: vec![*p],
                },
            );
        }
    }
    tree.into_values().collect()
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
