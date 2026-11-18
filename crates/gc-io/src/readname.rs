//! Herkunft einer Lesung aus ihrem Namen: Illumina (CASAVA 1.8+) und MGI/DNBSEQ.
//! Gebraucht für Read Groups und für optische Duplikate (Kachel und Koordinaten).

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Platform {
    Illumina,
    Mgi,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadOrigin {
    pub platform: Platform,
    pub flowcell: String,
    pub lane: u32,
    /// Kachel (Illumina) bzw. Feld aus Spalte und Reihe (MGI, `C001R002` → 1002).
    pub tile: u32,
    /// Lage auf der Kachel; bei MGI gibt es keine echten Koordinaten, dann `None`.
    pub xy: Option<(u32, u32)>,
}

/// Erkennt die Herkunft; `None`, wenn das Namensschema unbekannt ist.
pub fn parse(name: &[u8]) -> Option<ReadOrigin> {
    let name = std::str::from_utf8(crate::fastq::strip_mate_suffix(name)).ok()?;
    parse_illumina(name).or_else(|| parse_mgi(name))
}

/// `INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y` (7 Felder, ggf. mit UMI als 8. Feld).
fn parse_illumina(name: &str) -> Option<ReadOrigin> {
    let f: Vec<&str> = name.split(':').collect();
    if f.len() != 7 && f.len() != 8 {
        return None;
    }
    Some(ReadOrigin {
        platform: Platform::Illumina,
        flowcell: f[2].to_string(),
        lane: f[3].parse().ok()?,
        tile: f[4].parse().ok()?,
        xy: Some((f[5].parse().ok()?, f[6].parse().ok()?)),
    })
}

/// `V300012345L1C001R0020000123`: Flowcell, `L<lane>`, `C<spalte>R<reihe>`, laufende Nummer.
fn parse_mgi(name: &str) -> Option<ReadOrigin> {
    let l = name.rfind('L')?;
    let flowcell = &name[..l];
    if flowcell.len() < 6 || !flowcell.chars().all(|c| c.is_ascii_alphanumeric()) {
        return None;
    }
    let rest = &name[l + 1..];
    let c = rest.find('C')?;
    let lane: u32 = rest[..c].parse().ok()?;
    let rest = &rest[c + 1..];
    let r = rest.find('R')?;
    let col: u32 = rest[..r].parse().ok()?;
    let rest = &rest[r + 1..];
    if rest.len() < 4 || !rest.chars().all(|ch| ch.is_ascii_digit()) {
        return None;
    }
    let row: u32 = rest[..3].parse().ok()?;
    Some(ReadOrigin {
        platform: Platform::Mgi,
        flowcell: flowcell.to_string(),
        lane,
        tile: col * 1000 + row,
        xy: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn illumina_casava() {
        let o = parse(b"A00123:45:HXXXXDSXX:2:1101:12345:1000/1").unwrap();
        assert_eq!(o.platform, Platform::Illumina);
        assert_eq!(
            (o.flowcell.as_str(), o.lane, o.tile, o.xy),
            ("HXXXXDSXX", 2, 1101, Some((12345, 1000)))
        );
    }

    #[test]
    fn mgi_dnbseq() {
        let o = parse(b"V300012345L1C001R0020000123/2").unwrap();
        assert_eq!(o.platform, Platform::Mgi);
        assert_eq!(
            (o.flowcell.as_str(), o.lane, o.tile, o.xy),
            ("V300012345", 1, 1002, None)
        );
        let o = parse(b"E100012345L3C012R0340012345").unwrap();
        assert_eq!((o.lane, o.tile), (3, 12034));
    }

    #[test]
    fn unknown_scheme() {
        assert!(parse(b"SRR123.1").is_none());
        assert!(parse(b"read_42").is_none());
    }
}
