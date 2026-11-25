//! Ausrichtung gegen die Wahrheit im Read-Namen prüfen.
//!
//! Gezählt wird je primärem Datensatz (keine sekundären, keine supplementären): richtig, wenn
//! Contig gleich und das linke Ende samt weich und hart abgeschnittener Basen höchstens
//! `toleranz` Basen neben dem wahren liegt. Dazu je MAPQ die Zahl der Reads und der falschen
//! und die Zahl, die die MAPQ verspricht (Summe von 10^(-MAPQ/10)). Mehrere BAMs mit denselben
//! Reads werden Read für Read verglichen.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::bufread::MultiGzDecoder;
use serde::Serialize;

/// Wahrer Ort beider Reads eines Paars: (Contig-Nummer, linkes Ende 1-basiert).
pub type PairTruth = [(u32, u64); 2];

/// Liest `w.<nr>.<c1>.<p1>.<c2>.<p2>`.
pub fn parse_name(name: &[u8]) -> Option<(u64, PairTruth)> {
    let s = std::str::from_utf8(name).ok()?;
    let mut it = s.split('.');
    if it.next()? != "w" {
        return None;
    }
    let nr = it.next()?.parse().ok()?;
    let c1 = it.next()?.parse().ok()?;
    let p1 = it.next()?.parse().ok()?;
    let c2 = it.next()?.parse().ok()?;
    let p2 = it.next()?.parse().ok()?;
    if it.next().is_some() {
        return None;
    }
    Some((nr, [(c1, p1), (c2, p2)]))
}

fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut f = BufReader::new(File::open(path).with_context(|| path.display().to_string())?);
    let head = f.fill_buf()?;
    let gz = head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b;
    Ok(if gz {
        Box::new(BufReader::new(MultiGzDecoder::new(f)))
    } else {
        Box::new(f)
    })
}

/// Contig-Namen in der Reihenfolge der Simulation (`<prefix>.contigs.tsv`).
pub fn read_contigs(path: &Path) -> Result<Vec<String>> {
    let mut out = Vec::new();
    for line in open_text(path)?.lines() {
        let line = line?;
        let mut f = line.split('\t');
        let (Some(i), Some(name)) = (f.next(), f.next()) else {
            bail!("{}: Zeile ohne Nummer und Namen", path.display());
        };
        if i.parse::<usize>().ok() != Some(out.len()) {
            bail!("{}: Nummern nicht fortlaufend", path.display());
        }
        out.push(name.to_string());
    }
    Ok(out)
}

/// Bereiche aus einer BED-Datei, je Contig sortiert und verschmolzen.
pub struct Regions {
    pub name: String,
    by_contig: HashMap<String, Vec<(u64, u64)>>,
    pub bases: u64,
}

impl Regions {
    pub fn load(name: &str, path: &Path) -> Result<Self> {
        let mut by_contig: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        for line in open_text(path)?.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
                continue;
            }
            let mut f = line.split('\t');
            let (Some(c), Some(a), Some(b)) = (f.next(), f.next(), f.next()) else {
                bail!("{}: Zeile mit weniger als drei Spalten", path.display());
            };
            let (a, b): (u64, u64) = (a.parse()?, b.parse()?);
            by_contig.entry(c.to_string()).or_default().push((a, b));
        }
        let mut bases = 0;
        for v in by_contig.values_mut() {
            v.sort_unstable();
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(v.len());
            for &(a, b) in v.iter() {
                match merged.last_mut() {
                    Some(last) if a <= last.1 => last.1 = last.1.max(b),
                    _ => merged.push((a, b)),
                }
            }
            bases += merged.iter().map(|(a, b)| b - a).sum::<u64>();
            *v = merged;
        }
        Ok(Self {
            name: name.to_string(),
            by_contig,
            bases,
        })
    }

    /// Enthält der Bereich den 0-basierten Ort?
    pub fn contains(&self, contig: &str, pos0: u64) -> bool {
        let Some(v) = self.by_contig.get(contig) else {
            return false;
        };
        let i = v.partition_point(|&(a, _)| a <= pos0);
        i > 0 && pos0 < v[i - 1].1
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Outcome {
    Missing = 0,
    Unmapped = 1,
    Right = 2,
    Wrong = 3,
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct Counts {
    pub reads: u64,
    pub unmapped: u64,
    pub wrong: u64,
    pub mapq20: u64,
    pub wrong_mapq20: u64,
}

impl Counts {
    fn add(&mut self, o: Outcome, mapq: u8) {
        self.reads += 1;
        match o {
            Outcome::Unmapped => self.unmapped += 1,
            Outcome::Wrong => self.wrong += 1,
            _ => {}
        }
        if o != Outcome::Unmapped && mapq >= 20 && mapq != 255 {
            self.mapq20 += 1;
            if o == Outcome::Wrong {
                self.wrong_mapq20 += 1;
            }
        }
    }
}

/// Ergebnis für eine BAM.
pub struct Evaluation {
    pub label: String,
    /// Je MAPQ 0..=255: (Reads mit Ort, davon falsch).
    pub by_mapq: Vec<(u64, u64)>,
    pub all: Counts,
    pub regions: Vec<Counts>,
    pub records_skipped: u64,
    pub names_unknown: u64,
    /// Je Read (Paar-Nr. mal 2 plus Read 2): Ausgang und MAPQ.
    pub per_read: Vec<(Outcome, u8)>,
}

fn unclipped_start(rec: &[u8]) -> Option<i64> {
    let int = |at: usize| i32::from_le_bytes([rec[at], rec[at + 1], rec[at + 2], rec[at + 3]]);
    let pos = int(4) as i64;
    let l_name = rec[8] as usize;
    let n_cigar = u16::from_le_bytes([rec[12], rec[13]]) as usize;
    let mut clip = 0i64;
    for k in 0..n_cigar {
        let at = 32 + l_name + 4 * k;
        if at + 4 > rec.len() {
            return None;
        }
        let v = u32::from_le_bytes([rec[at], rec[at + 1], rec[at + 2], rec[at + 3]]);
        match v & 0xf {
            4 | 5 => clip += (v >> 4) as i64,
            _ => break,
        }
    }
    Some(pos - clip)
}

pub fn evaluate(
    label: &str,
    bam: &Path,
    contigs: &[String],
    regions: &[Regions],
    tolerance: u64,
    threads: usize,
    mut wrong_list: Option<&mut dyn std::io::Write>,
) -> Result<Evaluation> {
    let (mut reader, header) =
        gc_io::bam::open(bam, threads).with_context(|| bam.display().to_string())?;
    let bam_contigs: Vec<String> = gc_io::bam::contigs(&header)
        .into_iter()
        .map(|(n, _)| n)
        .collect();
    // Contig der Simulation -> Contig-Nummer in der BAM.
    let index: HashMap<&str, i32> = bam_contigs
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i as i32))
        .collect();
    let sim_to_bam: Vec<Option<i32>> = contigs
        .iter()
        .map(|n| index.get(n.as_str()).copied())
        .collect();

    let mut ev = Evaluation {
        label: label.to_string(),
        by_mapq: vec![(0, 0); 256],
        all: Counts::default(),
        regions: vec![Counts::default(); regions.len()],
        records_skipped: 0,
        names_unknown: 0,
        per_read: Vec::new(),
    };
    let mut buf = Vec::new();
    loop {
        buf.clear();
        if !gc_io::bam::read_encoded(&mut reader, &mut buf)? {
            break;
        }
        let rec = &buf[4..];
        if rec.len() < 32 {
            bail!("{}: Datensatz zu kurz", bam.display());
        }
        let flag = u16::from_le_bytes([rec[14], rec[15]]);
        if flag & (0x100 | 0x800) != 0 {
            ev.records_skipped += 1;
            continue;
        }
        let l_name = rec[8] as usize;
        let name = &rec[32..32 + l_name.saturating_sub(1)];
        let Some((nr, truth)) = parse_name(name) else {
            ev.names_unknown += 1;
            continue;
        };
        let second = flag & 0x80 != 0;
        let (tc, tp) = truth[second as usize];
        let mapq = rec[9];
        let outcome = if flag & 0x4 != 0 {
            Outcome::Unmapped
        } else {
            let tid = i32::from_le_bytes([rec[0], rec[1], rec[2], rec[3]]);
            let start = unclipped_start(rec).context("CIGAR abgeschnitten")? + 1;
            let same_contig = sim_to_bam.get(tc as usize).copied().flatten() == Some(tid);
            if same_contig && start.abs_diff(tp as i64) <= tolerance {
                Outcome::Right
            } else {
                if let Some(w) = wrong_list.as_deref_mut()
                    && mapq >= 20
                {
                    let found = usize::try_from(tid)
                        .ok()
                        .and_then(|t| bam_contigs.get(t))
                        .map_or("*", |s| s.as_str());
                    writeln!(
                        w,
                        "{label}\t{}\t{}\t{}:{}\t{found}:{start}\t{mapq}",
                        String::from_utf8_lossy(name),
                        1 + second as u8,
                        contigs.get(tc as usize).map_or("?", |s| s.as_str()),
                        tp
                    )?;
                }
                Outcome::Wrong
            }
        };
        if outcome != Outcome::Unmapped {
            ev.by_mapq[mapq as usize].0 += 1;
            if outcome == Outcome::Wrong {
                ev.by_mapq[mapq as usize].1 += 1;
            }
        }
        ev.all.add(outcome, mapq);
        if let Some(cname) = contigs.get(tc as usize) {
            for (k, r) in regions.iter().enumerate() {
                if r.contains(cname, tp - 1) {
                    ev.regions[k].add(outcome, mapq);
                }
            }
        }
        let slot = ((nr - 1) * 2 + second as u64) as usize;
        if ev.per_read.len() <= slot {
            ev.per_read.resize(slot + 1, (Outcome::Missing, 0));
        }
        ev.per_read[slot] = (outcome, mapq);
    }
    Ok(ev)
}

/// Stufen für Tabellen.
pub const BINS: &[(u8, u8)] = &[
    (0, 0),
    (1, 9),
    (10, 19),
    (20, 29),
    (30, 39),
    (40, 49),
    (50, 59),
    (60, 60),
];

#[derive(Debug, Serialize)]
pub struct Bin {
    pub mapq_from: u8,
    pub mapq_to: u8,
    pub reads: u64,
    pub wrong: u64,
    /// Was die MAPQ verspricht: Summe von 10^(-MAPQ/10).
    pub expected_wrong: f64,
}

pub fn bins(ev: &Evaluation) -> Vec<Bin> {
    BINS.iter()
        .map(|&(a, b)| {
            let mut bin = Bin {
                mapq_from: a,
                mapq_to: b,
                reads: 0,
                wrong: 0,
                expected_wrong: 0.0,
            };
            for q in a..=b {
                let (n, w) = ev.by_mapq[q as usize];
                bin.reads += n;
                bin.wrong += w;
                bin.expected_wrong += n as f64 * 10f64.powf(-(q as f64) / 10.0);
            }
            bin
        })
        .collect()
}

#[derive(Debug, Default, Serialize)]
pub struct Pairwise {
    pub a: String,
    pub b: String,
    pub both_right: u64,
    pub only_a_right: u64,
    pub only_b_right: u64,
    pub both_not_right: u64,
    /// Sicher (MAPQ ab 20) falsch in der einen, richtig in der anderen.
    pub a_confident_wrong_b_right: u64,
    pub b_confident_wrong_a_right: u64,
}

pub fn compare(a: &Evaluation, b: &Evaluation) -> Pairwise {
    let mut p = Pairwise {
        a: a.label.clone(),
        b: b.label.clone(),
        ..Default::default()
    };
    let n = a.per_read.len().max(b.per_read.len());
    let get = |e: &Evaluation, i: usize| e.per_read.get(i).copied();
    for i in 0..n {
        let (Some((oa, qa)), Some((ob, qb))) = (get(a, i), get(b, i)) else {
            continue;
        };
        if oa == Outcome::Missing && ob == Outcome::Missing {
            continue;
        }
        let ra = oa == Outcome::Right;
        let rb = ob == Outcome::Right;
        match (ra, rb) {
            (true, true) => p.both_right += 1,
            (true, false) => p.only_a_right += 1,
            (false, true) => p.only_b_right += 1,
            (false, false) => p.both_not_right += 1,
        }
        if oa == Outcome::Wrong && qa >= 20 && rb {
            p.a_confident_wrong_b_right += 1;
        }
        if ob == Outcome::Wrong && qb >= 20 && ra {
            p.b_confident_wrong_a_right += 1;
        }
    }
    p
}

pub fn parse_named(s: &str) -> std::result::Result<(String, PathBuf), String> {
    match s.split_once('=') {
        Some((n, p)) if !n.is_empty() && !p.is_empty() => Ok((n.to_string(), PathBuf::from(p))),
        _ => Err(format!("erwartet name=pfad, bekommen: {s}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_carry_the_truth() {
        assert_eq!(
            parse_name(b"w.12.3.1000.3.1250"),
            Some((12, [(3, 1000), (3, 1250)]))
        );
        assert_eq!(parse_name(b"x.1.0.1.0.1"), None);
        assert_eq!(parse_name(b"w.1.0.1.0"), None);
        assert_eq!(parse_name(b"w.1.0.1.0.1.9"), None);
    }

    fn record(pos: i32, cigar: &[(u32, u32)]) -> Vec<u8> {
        let mut r = vec![0u8; 32];
        r[4..8].copy_from_slice(&pos.to_le_bytes());
        r[8] = 2; // Name "a" plus NUL
        r[12..14].copy_from_slice(&(cigar.len() as u16).to_le_bytes());
        r.extend_from_slice(b"a\0");
        for &(len, op) in cigar {
            r.extend_from_slice(&(len << 4 | op).to_le_bytes());
        }
        r
    }

    #[test]
    fn clipped_bases_count_back_to_the_start() {
        assert_eq!(unclipped_start(&record(100, &[(150, 0)])), Some(100));
        assert_eq!(unclipped_start(&record(100, &[(5, 4), (145, 0)])), Some(95));
        assert_eq!(
            unclipped_start(&record(100, &[(2, 5), (3, 4), (145, 0)])),
            Some(95)
        );
    }

    #[test]
    fn regions_merge_and_contain() {
        let dir = std::env::temp_dir().join(format!("gv-bed-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bed = dir.join("r.bed");
        std::fs::write(&bed, "#x\nchr1\t10\t20\nchr1\t15\t30\nchr2\t0\t5\n").unwrap();
        let r = Regions::load("r", &bed).unwrap();
        assert_eq!(r.bases, 25);
        assert!(r.contains("chr1", 10));
        assert!(r.contains("chr1", 29));
        assert!(!r.contains("chr1", 30));
        assert!(!r.contains("chr1", 9));
        assert!(!r.contains("chr3", 1));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn pairwise_counts_each_read_once() {
        let mk = |v: Vec<(Outcome, u8)>| Evaluation {
            label: "x".into(),
            by_mapq: vec![(0, 0); 256],
            all: Counts::default(),
            regions: vec![],
            records_skipped: 0,
            names_unknown: 0,
            per_read: v,
        };
        let a = mk(vec![
            (Outcome::Right, 60),
            (Outcome::Wrong, 60),
            (Outcome::Right, 0),
            (Outcome::Missing, 0),
        ]);
        let b = mk(vec![
            (Outcome::Right, 60),
            (Outcome::Right, 60),
            (Outcome::Wrong, 30),
            (Outcome::Missing, 0),
        ]);
        let p = compare(&a, &b);
        assert_eq!(
            (
                p.both_right,
                p.only_a_right,
                p.only_b_right,
                p.both_not_right
            ),
            (1, 1, 1, 0)
        );
        assert_eq!(p.a_confident_wrong_b_right, 1);
        assert_eq!(p.b_confident_wrong_a_right, 1);
    }
}
