//! Bereiche eines Reads bereinigen, ordnen und bewerten (BWA-MEM: `mem_sort_dedup_patch`,
//! `mem_mark_primary_se`, `mem_approx_mapq_se`).

use crate::extend::Scoring;
use crate::region::Region;

/// Einstellungen. Die Vorgaben sind die von BWA-MEM.
#[derive(Debug, Clone, PartialEq)]
pub struct FinalizeOptions {
    /// Zwei Bereiche gelten als derselbe, wenn sie sich in Read und Referenz zu mehr als diesem
    /// Anteil des kürzeren überlappen.
    pub mask_level_redun: f32,
    /// Ein Bereich ist einem besseren untergeordnet, wenn sie sich im Read zu mindestens diesem
    /// Anteil des kürzeren überlappen.
    pub mask_level: f32,
    pub max_chain_gap: u64,
    pub min_seed_len: i32,
    /// MAPQ: unterhalb dieser Länge zählt die Länge nicht, darüber dämpft sie mit `fac / ln(l)`.
    pub mapq_coef_len: i32,
}

impl Default for FinalizeOptions {
    fn default() -> Self {
        Self {
            mask_level_redun: 0.95,
            mask_level: 0.5,
            max_chain_gap: 10_000,
            min_seed_len: 19,
            mapq_coef_len: 50,
        }
    }
}

/// Entfernt doppelte Bereiche und ordnet nach fallendem Score, dann Lage. Von zwei Bereichen,
/// die sich in Read und Referenz fast ganz überlappen, bleibt der mit dem höheren Score.
///
/// Noch nicht gebaut ist das Zusammenfügen zweier hintereinander liegender Bereiche über eine
/// große Lücke (`mem_patch_reg`).
pub fn sort_dedup(regions: &mut Vec<Region>, opt: &FinalizeOptions) {
    if regions.len() > 1 {
        regions.sort_by_key(|r| r.re);
        for i in 1..regions.len() {
            let p = regions[i];
            let prev = regions[i - 1];
            if p.contig != prev.contig || p.rb >= prev.re + opt.max_chain_gap {
                continue;
            }
            let mut j = i;
            while j > 0 {
                j -= 1;
                let q = regions[j];
                if p.contig != q.contig || p.rb >= q.re + opt.max_chain_gap {
                    break;
                }
                if q.qe == q.qb {
                    continue;
                }
                let or = q.re as i64 - p.rb as i64;
                let oq = if q.qb < p.qb {
                    q.qe - p.qb
                } else {
                    p.qe - q.qb
                };
                let mr = (q.re - q.rb).min(p.re - p.rb) as f32;
                let mq = (q.qe - q.qb).min(p.qe - p.qb) as f32;
                if or as f32 > opt.mask_level_redun * mr && oq as f32 > opt.mask_level_redun * mq {
                    if p.score < q.score {
                        regions[i].qe = regions[i].qb;
                        break;
                    }
                    regions[j].qe = regions[j].qb;
                }
            }
        }
        regions.retain(|r| r.qe > r.qb);
        regions.sort_by_key(|r| (std::cmp::Reverse(r.score), r.rb, r.qb));
        regions.dedup_by_key(|r| (r.score, r.rb, r.qb));
    }
}

/// Streuwert für die Reihenfolge gleich guter Bereiche (Thomas Wangs 64-Bit-Hash, wie in bwa).
fn hash_64(mut key: u64) -> u64 {
    key = key.wrapping_add(!(key << 32));
    key ^= key >> 22;
    key = key.wrapping_add(!(key << 13));
    key ^= key >> 8;
    key = key.wrapping_add(key << 3);
    key ^= key >> 15;
    key = key.wrapping_add(!(key << 27));
    key ^= key >> 31;
    key
}

/// Ordnet die Bereiche nach Score (bei Gleichstand nach einem Streuwert aus `read_id` und
/// Position, damit die Wahl nicht von der Lage im Genom abhängt) und bestimmt, welche primär
/// sind: Ein Bereich ist einem besseren untergeordnet, wenn er ihn im Read stark überlappt.
/// Der bessere merkt sich den Score des ersten solchen Bereichs (`sub`) und zählt die, die fast
/// so gut sind (`sub_n`).
pub fn mark_primary(regions: &mut [Region], read_id: u64, sc: &Scoring, opt: &FinalizeOptions) {
    if regions.is_empty() {
        return;
    }
    let mut keyed: Vec<(u64, Region)> = regions
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let mut r = *r;
            r.sub = 0;
            r.sub_n = 0;
            r.secondary = -1;
            (hash_64(read_id.wrapping_add(i as u64)), r)
        })
        .collect();
    keyed.sort_by_key(|(h, r)| (std::cmp::Reverse(r.score), *h));
    for (slot, (_, r)) in regions.iter_mut().zip(keyed) {
        *slot = r;
    }
    let near = (sc.matching + sc.mismatch)
        .max(sc.o_del + sc.e_del)
        .max(sc.o_ins + sc.e_ins);
    let mut primaries = vec![0usize];
    for i in 1..regions.len() {
        let mut parent = None;
        for &j in &primaries {
            let b_max = regions[j].qb.max(regions[i].qb);
            let e_min = regions[j].qe.min(regions[i].qe);
            if e_min > b_max {
                let min_l = (regions[i].qe - regions[i].qb).min(regions[j].qe - regions[j].qb);
                if (e_min - b_max) as f32 >= min_l as f32 * opt.mask_level {
                    if regions[j].sub == 0 {
                        regions[j].sub = regions[i].score;
                    }
                    if regions[j].score - regions[i].score <= near {
                        regions[j].sub_n += 1;
                    }
                    parent = Some(j);
                    break;
                }
            }
        }
        match parent {
            Some(j) => regions[i].secondary = j as i32,
            None => primaries.push(i),
        }
    }
}

/// Mapping-Qualität eines Bereichs, 0 bis 60 (bwa: `mem_approx_mapq_se`). Sie wächst mit dem
/// Abstand zum zweitbesten Score, sinkt mit der Zahl fast gleich guter Bereiche und mit dem
/// Anteil des Reads in sehr häufigen Seeds (`frac_rep`).
pub fn mapq(r: &Region, frac_rep: f32, sc: &Scoring, opt: &FinalizeOptions) -> u8 {
    let sub = if r.sub != 0 {
        r.sub
    } else {
        opt.min_seed_len * sc.matching
    };
    if sub >= r.score || r.score == 0 {
        return 0;
    }
    let l = (r.qe - r.qb).max((r.re - r.rb) as i32);
    let identity = 1.0
        - f64::from(l * sc.matching - r.score)
            / f64::from(sc.matching + sc.mismatch)
            / f64::from(l);
    let mut tmp = if l < opt.mapq_coef_len {
        1.0
    } else {
        // bwa hält den Faktor ln(50) = 3,912 in einer Ganzzahl und rechnet deshalb mit 3. Für
        // gleiche MAPQ wie bwa wird das hier bewusst nachgebildet.
        f64::from(opt.mapq_coef_len).ln().trunc() / f64::from(l).ln()
    };
    tmp *= identity * identity;
    let mut q =
        (6.02 * f64::from(r.score - sub) / f64::from(sc.matching) * tmp * tmp + 0.499) as i32;
    if r.sub_n > 0 {
        q -= (4.343 * f64::from(r.sub_n + 1).ln() + 0.499) as i32;
    }
    let q = q.clamp(0, 60);
    (f64::from(q) * (1.0 - f64::from(frac_rep)) + 0.499) as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    fn region(rb: u64, qb: i32, len: i32, score: i32) -> Region {
        Region {
            rb,
            re: rb + len as u64,
            qb,
            qe: qb + len,
            contig: 0,
            score,
            truesc: score,
            w: 100,
            seedcov: len,
            seedlen0: len,
            sub: 0,
            sub_n: 0,
            secondary: -1,
        }
    }

    #[test]
    fn duplicates_are_removed_and_the_better_one_stays() {
        let opt = FinalizeOptions::default();
        let mut v = vec![
            region(1000, 0, 100, 90),
            region(1000, 0, 100, 90), // genau gleich
            region(1001, 1, 99, 85),  // fast gleich, schlechter
            region(5000, 0, 100, 70), // anderer Ort
            region(1000, 0, 40, 40), // kürzer, im ersten enthalten: bleibt (nur 40 % Überlappung des ersten ...)
        ];
        sort_dedup(&mut v, &opt);
        let got: Vec<(u64, i32)> = v.iter().map(|r| (r.rb, r.score)).collect();
        assert_eq!(got, vec![(1000, 90), (5000, 70)]);
    }

    #[test]
    fn primary_secondary_and_mapq() {
        let (sc, opt) = (Scoring::default(), FinalizeOptions::default());
        // eindeutig: ein Bereich, ganzer Read von 150 Basen ohne Fehler
        let mut v = vec![region(1000, 0, 150, 150)];
        mark_primary(&mut v, 7, &sc, &opt);
        assert_eq!((v[0].secondary, v[0].sub), (-1, 0));
        assert_eq!(mapq(&v[0], 0.0, &sc, &opt), 60);
        // zwei gleich gute Orte: MAPQ 0, einer primär, einer untergeordnet
        let mut v = vec![region(1000, 0, 150, 150), region(9000, 0, 150, 150)];
        mark_primary(&mut v, 7, &sc, &opt);
        assert_eq!((v[0].secondary, v[1].secondary), (-1, 0));
        assert_eq!((v[0].sub, v[0].sub_n), (150, 1));
        assert_eq!(mapq(&v[0], 0.0, &sc, &opt), 0);
        // deutlich schlechterer zweiter Ort: MAPQ dazwischen
        let mut v = vec![region(9000, 0, 150, 120), region(1000, 0, 150, 150)];
        mark_primary(&mut v, 7, &sc, &opt);
        assert_eq!((v[0].score, v[0].sub, v[0].sub_n), (150, 120, 0));
        // 6,02 * 30 * (3 / ln 150)^2 = 64,7; gedeckelt auf 60
        assert_eq!(mapq(&v[0], 0.0, &sc, &opt), 60);
        // halber Read in sehr häufigen Seeds: halbe MAPQ
        assert_eq!(mapq(&v[0], 0.5, &sc, &opt), 30);
        // Bereiche in verschiedenen Teilen des Reads sind beide primär
        let mut v = vec![region(1000, 0, 70, 70), region(9000, 80, 70, 60)];
        mark_primary(&mut v, 7, &sc, &opt);
        assert_eq!((v[0].secondary, v[1].secondary), (-1, -1));
        // knapper Abstand: Score 150 gegen 145, ein fast gleich guter Bereich
        let mut v = vec![region(1000, 0, 150, 150), region(9000, 0, 150, 145)];
        mark_primary(&mut v, 7, &sc, &opt);
        // 6,02 * 5 * (3 / ln 150)^2 = 10,79 -> 11; minus (4,343 * ln 2 + 0,499) = 3 -> 8
        assert_eq!(mapq(&v[0], 0.0, &sc, &opt), 8);
        // Von bwa 0.7.19 für chr20 ausgegebene Werte bei 150 Basen: Score 145 gegen 130 ergibt 32,
        // Score 130 gegen 120 ergibt 19.
        let mut v = vec![region(1000, 0, 150, 145), region(9000, 0, 150, 130)];
        mark_primary(&mut v, 7, &sc, &opt);
        assert_eq!(mapq(&v[0], 0.0, &sc, &opt), 32);
        let mut v = vec![region(1000, 0, 150, 130), region(9000, 0, 150, 120)];
        mark_primary(&mut v, 7, &sc, &opt);
        assert_eq!(mapq(&v[0], 0.0, &sc, &opt), 19);
    }

    #[test]
    fn tie_order_depends_on_the_read_id_not_on_the_position() {
        let (sc, opt) = (Scoring::default(), FinalizeOptions::default());
        let mut first = 0;
        for id in 0..200u64 {
            let mut v = vec![region(1000, 0, 150, 150), region(9000, 0, 150, 150)];
            mark_primary(&mut v, id, &sc, &opt);
            first += usize::from(v[0].rb == 1000);
        }
        assert!((60..140).contains(&first), "{first}");
    }
}
