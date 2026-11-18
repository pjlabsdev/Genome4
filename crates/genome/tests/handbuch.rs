//! Hält das Handbuch (`docs/handbuch`) und das Programm zusammen:
//! Verweise zeigen auf vorhandene Dateien, jeder Unterbefehl hat eine Werkzeugseite,
//! jede Protokoll-Kennung aus dem Quelltext ist dokumentiert.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("..")
}

fn files_with_extension(dir: &Path, extension: &str, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).unwrap_or_else(|e| panic!("{}: {e}", dir.display())) {
        let path = entry.expect("Verzeichniseintrag").path();
        if path.is_dir() {
            files_with_extension(&path, extension, out);
        } else if path.extension().is_some_and(|e| e == extension) {
            out.push(path);
        }
    }
}

#[test]
fn verweise_im_handbuch_zeigen_auf_vorhandene_dateien() {
    let mut pages = Vec::new();
    files_with_extension(&root().join("docs").join("handbuch"), "md", &mut pages);
    assert!(
        pages.len() >= 20,
        "Handbuch unvollständig: {} Seiten",
        pages.len()
    );
    let mut broken = Vec::new();
    for page in &pages {
        let text = fs::read_to_string(page).expect("Seite lesbar");
        let mut rest = text.as_str();
        while let Some(at) = rest.find("](") {
            rest = &rest[at + 2..];
            let Some(end) = rest.find(')') else { break };
            let target = rest[..end].split('#').next().unwrap_or("");
            if target.is_empty() || target.contains("://") {
                continue;
            }
            let resolved = page.parent().expect("Ordner").join(target);
            if !resolved.exists() {
                broken.push(format!("{} → {target}", page.display()));
            }
        }
    }
    assert!(
        broken.is_empty(),
        "Verweise ins Leere:\n{}",
        broken.join("\n")
    );
}

#[test]
fn jeder_unterbefehl_hat_eine_werkzeugseite() {
    let output = Command::new(env!("CARGO_BIN_EXE_genome"))
        .arg("--help")
        .output()
        .expect("genome --help");
    let help = String::from_utf8(output.stdout).expect("UTF-8");
    let commands: Vec<&str> = help
        .lines()
        .skip_while(|l| !l.starts_with("Commands:"))
        .skip(1)
        .take_while(|l| !l.trim().is_empty())
        .filter_map(|l| l.split_whitespace().next())
        .filter(|c| *c != "help")
        .collect();
    assert!(commands.len() >= 3, "Unterbefehle nicht erkannt:\n{help}");
    for command in commands {
        let page = root()
            .join("docs")
            .join("handbuch")
            .join("werkzeuge")
            .join(format!("{command}.md"));
        assert!(page.exists(), "Werkzeugseite fehlt: {}", page.display());
    }
}

/// Zeichenketten der Form `WORT.WORT` in Großbuchstaben sind Protokoll-Kennungen.
fn is_code(s: &str) -> bool {
    let parts: Vec<&str> = s.split('.').collect();
    parts.len() >= 2
        && parts[0].chars().count() >= 3
        && parts
            .iter()
            .all(|p| !p.is_empty() && p.chars().all(|c| c.is_uppercase() || c == '_'))
}

#[test]
fn jede_protokoll_kennung_ist_dokumentiert() {
    let documented = fs::read_to_string(
        root()
            .join("docs")
            .join("handbuch")
            .join("ereignisse-und-fehlercodes.md"),
    )
    .expect("Seite der Kennungen");
    let mut sources = Vec::new();
    files_with_extension(&root().join("crates"), "rs", &mut sources);
    let mut missing = Vec::new();
    let mut seen = 0;
    for source in &sources {
        let text = fs::read_to_string(source).expect("Quelltext lesbar");
        for literal in text.split('"').skip(1).step_by(2) {
            if is_code(literal) {
                seen += 1;
                if !documented.contains(&format!("`{literal}`")) {
                    missing.push(format!("{literal} ({})", source.display()));
                }
            }
        }
    }
    assert!(seen >= 10, "Kennungen im Quelltext nicht erkannt");
    missing.sort();
    missing.dedup();
    assert!(
        missing.is_empty(),
        "nicht dokumentierte Kennungen:\n{}",
        missing.join("\n")
    );
}
