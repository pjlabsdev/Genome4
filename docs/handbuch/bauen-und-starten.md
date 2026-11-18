# Bauen und Starten

Stand: gebaut.

## Voraussetzung

Rust in der aktuellen stabilen Version (mindestens 1.88). Sonst nichts: keine Systembibliotheken,
kein C-Compiler für eigene Teile. Die CI baut mit der jeweils neuesten stabilen Version auf macOS
(arm64), Linux (x86-64 und arm64) und Windows (x86-64).

## Bauen und Testen

```bash
cargo build --release          # Programm: target/release/genome
cargo test --workspace         # alle Tests
cargo fmt --all --check        # Formatierung
cargo clippy --workspace --all-targets -- -D warnings
```

Die drei Prüfbefehle laufen in der CI bei jedem Pull Request auf allen vier Systemen.

## Aufruf

Ein Programm, Unterbefehle je Werkzeug:

```bash
genome --help
genome fastq-check R1.fq.gz R2.fq.gz
genome bam-roundtrip ein.bam aus.bam
genome index build referenz.fa referenz.gcfm
genome index locate referenz.gcfm ACGTTGCATTGACC
genome index count referenz.gcfm sequenzen.txt
```

Jedes Werkzeug ist unter [Werkzeuge](werkzeuge/README.md) beschrieben.

## Ausgaben

- Das Ergebnis kommt als JSON auf die Standardausgabe ([Ergebnis-Hülle](formate/ergebnis-huelle.md)).
- Lesbare Protokollzeilen kommen auf die Fehlerausgabe, dieselben Ereignisse als JSON-Zeilen in
  die Protokolldatei des Laufs ([Protokollierung](protokollierung.md)).
- Endet ein Lauf mit Fehler, ist der Rückgabewert des Programms ungleich 0.
