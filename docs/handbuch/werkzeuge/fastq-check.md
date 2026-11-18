# `genome fastq-check`

Stand: gebaut. Tests im Quelltext; eine Messung gegen fastp oder FastQC steht aus.

## Zweck
Prüft FASTQ-Dateien vollständig, bevor gerechnet wird, und gibt Kennzahlen aus. Fehler brechen
sofort mit Datei und Zeile ab, nicht erst nach Stunden im Alignment.

## Aufruf
```bash
genome fastq-check R1.fq.gz            # einzelne Datei
genome fastq-check R1.fq.gz R2.fq.gz   # gepaart
```

## Eingaben
FASTQ, unkomprimiert oder gzip (auch mehrteilig), erkannt am Dateianfang, nicht an der Endung.

## Ausgaben
JSON auf der Standardausgabe, Schema `genome.fastq-check.v1`:

| Feld | Bedeutung |
|---|---|
| `paired` | gepaart geprüft |
| `pairs` | Zahl der Paare (nur gepaart) |
| `r1`, `r2` | je Datei: `reads`, `bases`, `min_len`, `max_len`, `q30_pct` (Anteil Basen mit Qualität ab 30), `gc_pct`, `n_pct` |
| `origin` | aus dem Namen der ersten Lesung: `platform` (Illumina oder Mgi), `flowcell`, `lane`; `null`, wenn das Namensschema unbekannt ist |
| `seconds` | Laufzeit |

## Verfahren
[FASTQ-Prüfung und Herkunft der Lesungen](../verfahren/fastq-pruefung.md).

## Validierung
Tests in `gc-io` (`fastq.rs`, `readname.rs`): gültige Dateien, jede Fehlerart, vertauschte und
ungleich lange Paare, Illumina- und MGI-Namen.

## Grenzen
- Kein Trimmen, keine Adaptererkennung (kommt mit `fq2bam`).
- Kennzahlen sind einfache Summen; keine Verteilung je Zyklus.
- Schreibt noch kein Lauf-Protokoll über `gc-log`.

## Fehler und was zu tun ist
Die Meldung nennt Datei und Zeile. Die Fehlerarten stehen unter
[FASTQ-Prüfung](../verfahren/fastq-pruefung.md). Häufigste Ursache ist eine unvollständig kopierte
Datei: neu kopieren und Prüfsumme vergleichen.
