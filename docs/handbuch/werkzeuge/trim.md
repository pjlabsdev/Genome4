# `genome trim`

Stand: gebaut; der Vergleich Read für Read mit fastp 1.4.0 läuft mit `bench/m9`.

## Zweck
Read-Paare so trimmen wie fastp 1.4.0 mit `--detect_adapter_for_pe` und sonst Vorgaben, und die
behaltenen Paare als FASTQ schreiben. Das Werkzeug dient dem Vergleich mit fastp; in der Kette
trimmt `genome fq2bam --trimmen` ohne Zwischendateien.

## Aufruf
`genome trim --fastq R1 --fastq2 R2 --aus A1.fq --aus2 A2.fq [--log datei]`

| Option | Bedeutung | Vorgabe |
|---|---|---|
| `--fastq`, `--fastq2` | Eingabe Read 1 und Read 2 (roh, gzip oder BGZF) | |
| `--aus`, `--aus2` | Ausgabe Read 1 und Read 2, unkomprimiert | |
| `--log` | Lauf-Protokoll als JSON-Zeilen | `<aus>.log.jsonl` |

## Eingaben
Zwei FASTQ-Dateien mit Paaren in gleicher Reihenfolge ([FASTQ-Prüfung](../verfahren/fastq-pruefung.md)).

## Ausgaben
Die behaltenen Paare, in der Reihenfolge der Eingabe, Kopfzeile `@name kommentar` (ein
Leerzeichen als Trenner). Auf der Standardausgabe `genome.trim.v1` mit `trimmen` (Zähler wie im
Bericht von `fq2bam`) und `sekunden`.

## Verfahren
[Trimmen](../verfahren/trimmen.md). Stufen `fq2bam.trimmen_vorbereiten` und `trim.schreiben`.

## Validierung
Regeltests in `gc-io`. Vergleich mit fastp: `bench/m9/trim-vergleich.py`, Ergebnis folgt.

## Grenzen
Nur Paare. Ein Kern beim Schreiben.

## Fehler und was zu tun ist
Fehlerhafte FASTQ brechen mit `LAUF.FEHLGESCHLAGEN` und der Ursache im Feld `fehler` ab;
`genome fastq-check` zeigt die Stelle.

## Bericht
Keiner; die Zähler stehen in der JSON-Ausgabe.

## Literatur
Chen et al. 2018 ([Literatur](../literatur.md)).
