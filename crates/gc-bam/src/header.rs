//! Kopftext einer BAM nach dem Sortieren.

/// Setzt im Kopftext die Sortierung auf "nach Ort" und hängt eine `@PG`-Zeile für diesen Lauf an.
///
/// - Gibt es eine `@HD`-Zeile, wird ihr Feld `SO` auf `coordinate` gesetzt (oder angehängt) und
///   ein Feld `SS` entfernt, weil es die alte Sortierung beschreibt. Sonst kommt
///   `@HD\tVN:1.6\tSO:coordinate` an den Anfang.
/// - Die neue `@PG`-Zeile verweist mit `PP` auf die letzte vorhandene, wie bei samtools. Ihre
///   Kennung ist `id`; ist die schon vergeben, wird `.1`, `.2` und so weiter angehängt.
pub fn sorted_header(text: &str, id: &str, program_version: &str, command_line: &str) -> String {
    let mut lines: Vec<String> = text
        .lines()
        .filter(|l| !l.is_empty())
        .map(str::to_string)
        .collect();
    match lines.first_mut() {
        Some(first) if first.starts_with("@HD") => {
            let mut fields: Vec<String> = first
                .split('\t')
                .filter(|f| !f.starts_with("SS:"))
                .map(str::to_string)
                .collect();
            match fields.iter_mut().find(|f| f.starts_with("SO:")) {
                Some(f) => *f = "SO:coordinate".to_string(),
                None => fields.push("SO:coordinate".to_string()),
            }
            *first = fields.join("\t");
        }
        _ => lines.insert(0, "@HD\tVN:1.6\tSO:coordinate".to_string()),
    }
    let mut out = lines.join("\n");
    out.push('\n');
    with_program_line(&out, id, program_version, command_line)
}

/// Hängt an den Kopftext eine `@PG`-Zeile für diesen Lauf an. Sie verweist mit `PP` auf die
/// letzte vorhandene, wie bei samtools. Ihre Kennung ist `id`; ist die schon vergeben, wird `.1`,
/// `.2` und so weiter angehängt. Alle anderen Zeilen bleiben unverändert.
pub fn with_program_line(
    text: &str,
    id: &str,
    program_version: &str,
    command_line: &str,
) -> String {
    let mut lines: Vec<String> = text
        .lines()
        .filter(|l| !l.is_empty())
        .map(str::to_string)
        .collect();
    let ids: Vec<&str> = lines
        .iter()
        .filter(|l| l.starts_with("@PG\t"))
        .filter_map(|l| l.split('\t').find_map(|f| f.strip_prefix("ID:")))
        .collect();
    let mut own = id.to_string();
    let mut n = 0;
    while ids.contains(&own.as_str()) {
        n += 1;
        own = format!("{id}.{n}");
    }
    let mut pg = format!("@PG\tID:{own}\tPN:genome");
    if let Some(last) = ids.last() {
        pg.push_str(&format!("\tPP:{last}"));
    }
    // Tabulatoren und Zeilenenden im Aufruf würden die Zeile zerlegen.
    let command_line: String = command_line
        .chars()
        .map(|c| {
            if c == '\t' || c == '\n' || c == '\r' {
                ' '
            } else {
                c
            }
        })
        .collect();
    pg.push_str(&format!("\tVN:{program_version}\tCL:{command_line}"));
    lines.push(pg);
    let mut out = lines.join("\n");
    out.push('\n');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_gets_sort_order_and_program_line() {
        // ohne @HD
        let h = sorted_header(
            "@SQ\tSN:c\tLN:10\n@PG\tID:genome\tPN:genome\n",
            "genome.sort",
            "1",
            "x y",
        );
        assert_eq!(
            h,
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:10\n@PG\tID:genome\tPN:genome\n\
             @PG\tID:genome.sort\tPN:genome\tPP:genome\tVN:1\tCL:x y\n"
        );
        // @HD mit anderer Sortierung und SS; Kennung schon vergeben; Tabulator im Aufruf
        let h = sorted_header(
            "@HD\tVN:1.5\tSO:queryname\tSS:queryname:natural\tGO:query\n@PG\tID:genome.sort\n",
            "genome.sort",
            "2",
            "a\tb",
        );
        assert_eq!(
            h,
            "@HD\tVN:1.5\tSO:coordinate\tGO:query\n@PG\tID:genome.sort\n\
             @PG\tID:genome.sort.1\tPN:genome\tPP:genome.sort\tVN:2\tCL:a b\n"
        );
        // @HD ohne SO; kein @PG vorher
        let h = sorted_header("@HD\tVN:1.6\n", "s", "1", "c");
        assert_eq!(
            h,
            "@HD\tVN:1.6\tSO:coordinate\n@PG\tID:s\tPN:genome\tVN:1\tCL:c\n"
        );
        // leerer Kopf
        let h = sorted_header("", "s", "1", "c");
        assert_eq!(
            h,
            "@HD\tVN:1.6\tSO:coordinate\n@PG\tID:s\tPN:genome\tVN:1\tCL:c\n"
        );
    }
}
