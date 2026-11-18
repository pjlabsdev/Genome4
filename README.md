# Genome 4.0

Genome 4.0 ist die Neufassung von Genome komplett in Rust: alle Analysewerkzeuge
(von FASTQ bis zum Bericht) und die App selbst. Läuft auf macOS, Linux und Windows.

Privates Repository. Einzelne Werkzeuge werden später als eigene öffentliche
Repositories mit DOI veröffentlicht (wie hla-wgs und brca12-wgs).

- Plan: [`docs/plan/Genome-4.0-Plan.md`](docs/plan/Genome-4.0-Plan.md)
- Erster Schritt: FASTQ → BAM (`genome fq2bam`), Abschnitt 0 des Plans

## Aufbau

| Crate | Aufgabe |
|---|---|
| `gc-index` | Suffix-Array (SA-IS) und FM-Index über das Referenzgenom |
| `gc-io` | FASTQ, FASTA, BAM, BGZF lesen und schreiben |
| `gc-log` | Ereignisstrom und Lauf-Protokoll für alle Werkzeuge |
| `gc-schema` | versionierte JSON-Ausgaben (`genome.<werkzeug>.vN`) |
| `genome` | ein Programm: Werkzeuge, Arbeitsprozesse, später die Oberfläche |

## Bauen und testen

```bash
cargo build --release
cargo test --workspace
```

Heute verfügbar (Meilenstein M1):

```bash
genome fastq-check R1.fq.gz R2.fq.gz   # prüft Paare vollständig, Kennzahlen als JSON
genome bam-roundtrip in.bam out.bam    # Rundreise-Test, parallele BGZF-Kompression
```

Meilenstein M2, Teil 1 (siehe [`docs/bench/M2-fm-index.md`](docs/bench/M2-fm-index.md)):

```bash
genome index build ref.fa ref.gcfm --contig chr20   # FM-Index mit Selbstprüfung und Lauf-Protokoll
genome index locate ref.gcfm ACGTTGCA...            # Treffer mit Contig, Position, Strang
genome index count ref.gcfm sequenzen.txt           # exakte Trefferzahl je Zeile
```

## Regeln

Siehe [`AGENTS.md`](AGENTS.md).
