# `genome fq2bam`

Stand: geplant. Diese Seite gibt nur wieder, was der Plan (Abschnitt 0) festlegt. Verfahren,
Zahlen und Grenzen kommen mit den Meilensteinen dazu.

## Zweck
Ein Programm vom gepaarten FASTQ bis zur sortierten, duplikatmarkierten, indizierten BAM (oder
CRAM) mit Qualitätsauswertung. Ersetzt in Genome 3.x fastp, `bwa mem` 0.7.19 und samtools
(fixmate, sort, markdup, index, flagstat).

## Meilensteine

| M | Inhalt | Stand |
|---|---|---|
| M1 | FASTQ und BAM lesen und schreiben, Basismessung der heutigen Kette | gebaut; Basismessung läuft |
| M2 | Suffix-Array und FM-Index, Protokoll-Grundlage | gebaut und geprüft |
| M3 | Seeds, Chaining, SIMD-Erweiterung, Einzel-Lesungen auf chr20 | offen |
| M4 | Paarung, MAPQ, supplementäre Alignments, ALT, ganzes Genom | offen |
| M5 | Sortieren, Duplikate, Index, QC, Trimmen; erster Werkzeugbericht; Schattenbetrieb in 3.x | Sortieren ([`genome sort`](sort.md)) und Duplikate ([`genome markdup`](markdup.md)) gebaut und geprüft, Rest offen |
| M6 | eingefrorene Version, volle Testrunde, Forschungsbericht | offen |

## Vergleichsmaßstab
`bwa mem` 0.7.19 mit samtools 1.23.1, wie Genome 3.x sie heute aufruft
([Bestand](../../bench/M1-Bestand-Genome3-FASTQ-BAM.md)). Messziele: Plan Abschnitt 0.4.

## Literatur
Li 2013; Li & Durbin 2009; Vasimuddin et al. 2019. Siehe [Literatur](../literatur.md).
