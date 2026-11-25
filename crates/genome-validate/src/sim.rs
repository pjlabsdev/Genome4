//! Read-Paare mit bekannter Herkunft erzeugen.
//!
//! Ablauf je Paar:
//! 1. Eine Probe mit zwei Haplotypen: auf die Referenz werden SNVs und kleine Indels gesetzt,
//!    jede Variante auf Haplotyp 1, 2 oder beiden. Die Varianten stehen in
//!    `<prefix>.varianten.tsv`, damit dieselbe Probe später auch Variantenfinder prüfen kann.
//! 2. Ein Stück eines Haplotyps wird gezogen, Länge normalverteilt, von beiden Enden gelesen.
//!    Zu jeder Base des Stücks ist der Ort auf der Referenz bekannt; daraus folgt das wahre
//!    linke Ende beider Reads.
//! 3. Lesefehler: je Base wird zuerst eine Qualität gezogen (zum Read-Ende hin öfter niedrig),
//!    dann tritt ein Fehler mit genau der Wahrscheinlichkeit auf, die die Qualität verspricht.
//!    Dazu seltene Einfügungen und Lücken von einer Base.
//!
//! Der Name trägt die Wahrheit: `w.<nr>.<contig1>.<pos1>.<contig2>.<pos2>`, Contig als Nummer
//! in der Reihenfolge der FASTA (Namen dazu in `<prefix>.contigs.tsv`), Positionen 1-basiert,
//! jeweils das linke Ende des Reads auf dem Vorwärtsstrang.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;

use crate::rng::Rng;

/// Fehlerprofil eines Geräts: die vier Qualitätsstufen und wie oft die niedrigen vorkommen.
#[derive(Debug, Clone, Copy, serde::Serialize)]
pub struct Profile {
    pub name: &'static str,
    /// Qualitäten der Stufen, von niedrig nach hoch.
    pub levels: [u8; 4],
    /// Anteil der Basen unter der höchsten Stufe am Read-Anfang und am Read-Ende.
    pub low_start: f64,
    pub low_end: f64,
    /// Einfügung oder Lücke von einer Base, je Base.
    pub indel_rate: f64,
}

/// NovaSeq meldet nur vier Qualitätsstufen (2, 12, 23, 37).
pub const ILLUMINA: Profile = Profile {
    name: "illumina",
    levels: [2, 12, 23, 37],
    low_start: 0.03,
    low_end: 0.15,
    indel_rate: 0.000_05,
};

/// DNBSEQ (MGI) mit seinen typischen Stufen.
pub const MGI: Profile = Profile {
    name: "mgi",
    levels: [5, 15, 28, 40],
    low_start: 0.02,
    low_end: 0.12,
    indel_rate: 0.000_02,
};

#[derive(Debug, Clone, serde::Serialize)]
pub struct Settings {
    pub pairs: u64,
    pub length: usize,
    pub insert_mean: f64,
    pub insert_sd: f64,
    pub seed: u64,
    pub snv_rate: f64,
    pub indel_rate: f64,
    pub max_indel: u64,
    pub profile: Profile,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Snv(u8),
    /// Fehlende Referenzbasen ab `pos`.
    Del(u32),
    /// Eingefügte Basen vor der Referenzbase `pos`; die Basen folgen aus der Saat.
    Ins(u32, u64),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Variant {
    pub pos: u64,
    /// Bit 0: Haplotyp 1, Bit 1: Haplotyp 2.
    pub haps: u8,
    pub kind: Kind,
}

const BASES: &[u8; 4] = b"ACGT";

fn inserted(seed: u64, len: u32) -> Vec<u8> {
    let mut r = Rng::new(seed);
    (0..len).map(|_| BASES[r.below(4) as usize]).collect()
}

/// Setzt Varianten auf einen Contig. Orte mit N bekommen keine.
pub fn place_variants(seq: &[u8], s: &Settings, rng: &mut Rng) -> Vec<Variant> {
    let mut out = Vec::new();
    let rate = s.snv_rate + s.indel_rate;
    if rate <= 0.0 {
        return out;
    }
    let mut pos = rng.gap(rate) - 1;
    let mut blocked_until = 0u64;
    while (pos as usize) < seq.len() {
        let base = seq[pos as usize].to_ascii_uppercase();
        if pos >= blocked_until && b"ACGT".contains(&base) {
            // 60 % mischerbig (je zur Hälfte auf einem der beiden), 40 % reinerbig.
            let u = rng.unit();
            let haps = if u < 0.3 {
                1
            } else if u < 0.6 {
                2
            } else {
                3
            };
            let kind = if rng.unit() < s.snv_rate / rate {
                let mut alt = BASES[rng.below(4) as usize];
                while alt == base {
                    alt = BASES[rng.below(4) as usize];
                }
                Kind::Snv(alt)
            } else {
                // Länge 1 bis max, kurze viel häufiger (geometrisch mit Hälfte je Stufe).
                let mut len = 1;
                while len < s.max_indel && rng.unit() < 0.5 {
                    len += 1;
                }
                if rng.unit() < 0.5 {
                    Kind::Del(len as u32)
                } else {
                    Kind::Ins(len as u32, rng.next_u64())
                }
            };
            if let Kind::Del(len) = kind {
                blocked_until = pos + len as u64 + 1;
            } else {
                blocked_until = pos + 1;
            }
            out.push(Variant { pos, haps, kind });
        }
        pos += rng.gap(rate);
    }
    out
}

/// Baut ab Referenzort `start` ein Stück von `len` Basen des Haplotyps `hap` (0 oder 1).
/// Gibt die Basen und zu jeder Base den Referenzort zurück, oder `None`, wenn das Stück über
/// das Contig-Ende läuft oder N enthält.
pub fn haplotype_piece(
    seq: &[u8],
    variants: &[Variant],
    hap: u8,
    start: u64,
    len: usize,
) -> Option<(Vec<u8>, Vec<u64>)> {
    let bit = 1u8 << hap;
    let mut bases = Vec::with_capacity(len);
    let mut coords = Vec::with_capacity(len);
    let mut v = variants.partition_point(|x| x.pos < start);
    // Beginnt das Stück in einer Lücke dieses Haplotyps, wird sie übersprungen.
    let mut p = start;
    if let Some(prev) = variants[..v].iter().rev().find(|x| x.haps & bit != 0)
        && let Kind::Del(l) = prev.kind
        && prev.pos + l as u64 > p
    {
        p = prev.pos + l as u64;
        v = variants.partition_point(|x| x.pos < p);
    }
    while bases.len() < len {
        if p as usize >= seq.len() {
            return None;
        }
        while v < variants.len() && variants[v].pos < p {
            v += 1;
        }
        let here = variants
            .get(v)
            .filter(|x| x.pos == p && x.haps & bit != 0)
            .map(|x| x.kind);
        match here {
            Some(Kind::Del(l)) => {
                p += l as u64;
                continue;
            }
            Some(Kind::Ins(l, seed)) => {
                for b in inserted(seed, l) {
                    bases.push(b);
                    coords.push(p);
                }
            }
            _ => {}
        }
        let base = match here {
            Some(Kind::Snv(alt)) => alt,
            _ => seq[p as usize].to_ascii_uppercase(),
        };
        if !b"ACGT".contains(&base) {
            return None;
        }
        bases.push(base);
        coords.push(p);
        p += 1;
    }
    bases.truncate(len);
    coords.truncate(len);
    Some((bases, coords))
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        x => x,
    }
}

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b)).collect()
}

/// Lesefehler und Qualitäten auf einen Read in Leserichtung. Verbraucht eine Lücke die letzte
/// Base des Stücks, liest das Gerät in den Adapter: dort steht eine zufällige Base.
pub fn sequence(source: &[u8], length: usize, p: &Profile, rng: &mut Rng) -> (Vec<u8>, Vec<u8>) {
    let mut seq = Vec::with_capacity(length);
    let mut qual = Vec::with_capacity(length);
    let mut i = 0;
    while seq.len() < length {
        let cycle = seq.len() as f64 / length as f64;
        let low = p.low_start + (p.low_end - p.low_start) * cycle * cycle;
        let level = if rng.unit() < low {
            // Unter den niedrigen Stufen ist die mittlere am häufigsten.
            let u = rng.unit();
            if u < 0.15 {
                0
            } else if u < 0.45 {
                1
            } else {
                2
            }
        } else {
            3
        };
        let q = p.levels[level];
        let err = 10f64.powf(-(q as f64) / 10.0);
        if rng.unit() < p.indel_rate {
            if rng.unit() < 0.5 {
                // Lücke: Quellbase fällt aus.
                i += 1;
                continue;
            }
            // Einfügung: zufällige Base, Quelle bleibt stehen.
            seq.push(BASES[rng.below(4) as usize]);
            qual.push(q + 33);
            continue;
        }
        let mut b = source
            .get(i)
            .copied()
            .unwrap_or_else(|| BASES[rng.below(4) as usize]);
        // Bei Qualität q ist die Base mit Wahrscheinlichkeit 10^(-q/10) falsch (bei Q2 also
        // fast immer zufällig).
        if rng.unit() < err {
            let mut x = BASES[rng.below(4) as usize];
            while x == b {
                x = BASES[rng.below(4) as usize];
            }
            b = x;
        }
        seq.push(b);
        qual.push(q + 33);
        i += 1;
    }
    (seq, qual)
}

pub struct Summary {
    pub contigs: usize,
    pub variants: u64,
    pub snvs: u64,
    pub indels: u64,
    pub pairs: u64,
    pub rejected: u64,
}

pub fn run(reference: &Path, prefix: &Path, s: &Settings) -> Result<Summary> {
    if s.length == 0 || s.insert_mean < s.length as f64 {
        bail!("Stücke müssen mindestens so lang sein wie ein Read");
    }
    let mut reader = gc_io::fasta::open(reference)?;
    let mut rec = gc_io::fasta::Record::default();
    let mut names = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    while reader.read_record(&mut rec)? {
        names.push(rec.name.clone());
        seqs.push(std::mem::take(&mut rec.seq));
    }
    if seqs.is_empty() {
        bail!("{}: keine Einträge", reference.display());
    }
    let path = |ext: &str| -> PathBuf {
        let mut p = prefix.as_os_str().to_owned();
        p.push(ext);
        PathBuf::from(p)
    };

    let mut contigs = BufWriter::new(File::create(path(".contigs.tsv"))?);
    for (i, (n, sq)) in names.iter().zip(&seqs).enumerate() {
        writeln!(contigs, "{i}\t{n}\t{}", sq.len())?;
    }
    contigs.flush()?;

    // Varianten: eigene Saat je Contig, damit sie nicht von der Reihenfolge der Reads abhängen.
    let mut variants = Vec::with_capacity(seqs.len());
    let mut vt = BufWriter::new(File::create(path(".varianten.tsv"))?);
    writeln!(vt, "#contig\tpos\tref\talt\thaplotypen")?;
    let (mut n_snv, mut n_indel) = (0u64, 0u64);
    for (i, sq) in seqs.iter().enumerate() {
        let mut r = Rng::new(s.seed ^ (0x5eed_0000 + i as u64));
        let vs = place_variants(sq, s, &mut r);
        for v in &vs {
            let gt = match v.haps {
                1 => "1|0",
                2 => "0|1",
                _ => "1|1",
            };
            let p = v.pos as usize;
            // Schreibweise wie VCF: Indels mit der Base davor als Anker.
            let (pos1, refa, alt) = match v.kind {
                Kind::Snv(a) => {
                    n_snv += 1;
                    (p + 1, vec![sq[p]], vec![a])
                }
                Kind::Del(l) if p > 0 => {
                    n_indel += 1;
                    (p, sq[p - 1..p + l as usize].to_vec(), vec![sq[p - 1]])
                }
                Kind::Ins(l, seed) if p > 0 => {
                    n_indel += 1;
                    let mut a = vec![sq[p - 1]];
                    a.extend(inserted(seed, l));
                    (p, vec![sq[p - 1]], a)
                }
                _ => continue,
            };
            writeln!(
                vt,
                "{}\t{}\t{}\t{}\t{}",
                names[i],
                pos1,
                String::from_utf8_lossy(&refa).to_ascii_uppercase(),
                String::from_utf8_lossy(&alt).to_ascii_uppercase(),
                gt
            )?;
        }
        variants.push(vs);
    }
    vt.flush()?;

    let ends: Vec<u64> = seqs
        .iter()
        .scan(0u64, |acc, sq| {
            *acc += sq.len() as u64;
            Some(*acc)
        })
        .collect();
    let total = *ends.last().unwrap();
    let open = |p: PathBuf| -> Result<GzEncoder<BufWriter<File>>> {
        let f = File::create(&p).with_context(|| p.display().to_string())?;
        Ok(GzEncoder::new(
            BufWriter::with_capacity(1 << 20, f),
            Compression::fast(),
        ))
    };
    let mut f1 = open(path("_1.fq.gz"))?;
    let mut f2 = open(path("_2.fq.gz"))?;
    let mut rng = Rng::new(s.seed);
    let (mut done, mut rejected) = (0u64, 0u64);
    while done < s.pairs {
        let size = (s.insert_mean + s.insert_sd * rng.gauss(0.0, 1.0)).round();
        let size = (size.max(s.length as f64)) as usize;
        let g = rng.below(total);
        let c = ends.partition_point(|&e| e <= g);
        let start = g - if c == 0 { 0 } else { ends[c - 1] };
        let hap = rng.below(2) as u8;
        let swap = rng.below(2) == 1;
        let Some((frag, coords)) = haplotype_piece(&seqs[c], &variants[c], hap, start, size) else {
            rejected += 1;
            continue;
        };
        let pos_left = coords[0] + 1;
        let pos_right = coords[size - s.length] + 1;
        let left = frag.clone();
        let right = reverse_complement(&frag);
        let (src1, pos1, src2, pos2) = if swap {
            (right, pos_right, left, pos_left)
        } else {
            (left, pos_left, right, pos_right)
        };
        let (s1, q1) = sequence(&src1, s.length, &s.profile, &mut rng);
        let (s2, q2) = sequence(&src2, s.length, &s.profile, &mut rng);
        done += 1;
        let name = format!("w.{done}.{c}.{pos1}.{c}.{pos2}");
        for (f, sq, q) in [(&mut f1, &s1, &q1), (&mut f2, &s2, &q2)] {
            f.write_all(b"@")?;
            f.write_all(name.as_bytes())?;
            f.write_all(b"\n")?;
            f.write_all(sq)?;
            f.write_all(b"\n+\n")?;
            f.write_all(q)?;
            f.write_all(b"\n")?;
        }
    }
    f1.finish()?.flush()?;
    f2.finish()?.flush()?;
    Ok(Summary {
        contigs: seqs.len(),
        variants: n_snv + n_indel,
        snvs: n_snv,
        indels: n_indel,
        pairs: done,
        rejected,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(pos: u64, haps: u8, kind: Kind) -> Variant {
        Variant { pos, haps, kind }
    }

    #[test]
    fn piece_without_variants_is_the_reference() {
        let seq = b"ACGTACGTAC";
        let (b, c) = haplotype_piece(seq, &[], 0, 2, 5).unwrap();
        assert_eq!(b, b"GTACG");
        assert_eq!(c, vec![2, 3, 4, 5, 6]);
    }

    #[test]
    fn variants_change_only_their_haplotype() {
        let seq = b"AAAAAAAAAA";
        let vs = [
            v(2, 1, Kind::Snv(b'C')),
            v(4, 2, Kind::Del(2)),
            v(7, 3, Kind::Ins(2, 9)),
        ];
        let ins = inserted(9, 2);
        let (b0, c0) = haplotype_piece(seq, &vs, 0, 0, 10).unwrap();
        let mut want = b"AACAAAA".to_vec();
        want.extend(&ins);
        want.push(b'A');
        assert_eq!(b0, want);
        assert_eq!(c0, vec![0, 1, 2, 3, 4, 5, 6, 7, 7, 7]);
        let (b1, c1) = haplotype_piece(seq, &vs, 1, 0, 8).unwrap();
        let mut want = b"AAAAA".to_vec();
        want.extend(&ins);
        want.push(b'A');
        assert_eq!(b1, want);
        // Ort 4 und 5 fehlen, die Einfügung sitzt vor Ort 7.
        assert_eq!(c1, vec![0, 1, 2, 3, 6, 7, 7, 7]);
    }

    #[test]
    fn piece_starting_inside_a_deletion_skips_it() {
        let seq = b"ACGTACGTAC";
        let vs = [v(2, 3, Kind::Del(3))];
        let (_, c) = haplotype_piece(seq, &vs, 0, 3, 3).unwrap();
        assert_eq!(c, vec![5, 6, 7]);
    }

    #[test]
    fn pieces_with_n_or_past_the_end_are_refused() {
        assert!(haplotype_piece(b"ACGNACGT", &[], 0, 0, 5).is_none());
        assert!(haplotype_piece(b"ACGT", &[], 0, 2, 5).is_none());
    }

    #[test]
    fn quality_promises_the_error_rate() {
        // Ohne Indels: gemessene Fehler je Stufe passen zur gemeldeten Qualität.
        let p = Profile {
            indel_rate: 0.0,
            low_start: 0.3,
            low_end: 0.3,
            ..ILLUMINA
        };
        let src = vec![b'A'; 200];
        let mut rng = Rng::new(5);
        let mut seen = [(0u64, 0u64); 4];
        for _ in 0..20_000 {
            let (s, q) = sequence(&src, 150, &p, &mut rng);
            for (b, q) in s.iter().zip(&q) {
                let l = p.levels.iter().position(|&x| x + 33 == *q).unwrap();
                seen[l].0 += 1;
                seen[l].1 += u64::from(*b != b'A');
            }
        }
        for (l, (n, e)) in seen.iter().enumerate() {
            // Die neue Base ist immer eine andere, also ist die Rate genau 10^(-q/10).
            let want = 10f64.powf(-(p.levels[l] as f64) / 10.0);
            let got = *e as f64 / *n as f64;
            assert!(
                (got - want).abs() < 0.1 * want + 0.0005,
                "Stufe {l}: {got} statt {want}"
            );
        }
    }

    #[test]
    fn same_seed_same_reads() {
        let src = b"ACGTTGCAACGTTGCAACGTTGCA".repeat(20);
        let a = sequence(&src, 150, &ILLUMINA, &mut Rng::new(11));
        let b = sequence(&src, 150, &ILLUMINA, &mut Rng::new(11));
        assert_eq!(a, b);
    }
}
