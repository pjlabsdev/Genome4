# Werkzeuge

Jedes Werkzeug hat eine Seite nach derselben [Vorlage](vorlage.md).

## Gebaut

| Befehl | Aufgabe | Stand | Seite |
|---|---|---|---|
| `genome fastq-check` | FASTQ vollständig prüfen, Kennzahlen | gebaut | [fastq-check](fastq-check.md) |
| `genome bam-roundtrip` | BAM lesen und unverändert neu schreiben | gebaut und geprüft | [bam-roundtrip](bam-roundtrip.md) |
| `genome sort` | BAM nach Ort sortieren, mit festem Speicher | gebaut und geprüft | [sort](sort.md) |
| `genome markdup` | Duplikate markieren, Verfahren wie samtools | gebaut und geprüft | [markdup](markdup.md) |
| `genome bam-index` | Index `.bai` zu einer sortierten BAM | gebaut und geprüft | [bam-index](bam-index.md) |
| `genome flagstat` | Datensätze nach Flags zählen | gebaut und geprüft | [flagstat](flagstat.md) |
| `genome fq2bam` | FASTQ zur sortierten, duplikatmarkierten, indizierten BAM in einem Lauf | erste Fassung, auf chr20 geprüft | [fq2bam](fq2bam.md) |
| `genome trim` | Paare trimmen wie fastp 1.4.0, als FASTQ (für den Vergleich) | gebaut und geprüft | [trim](trim.md) |
| `genome index` | FM-Index bauen, prüfen, abfragen | gebaut und geprüft | [index](index.md) |

## Geplant

Reihenfolge und Inhalt stehen im [Plan](../../plan/Genome-4.0-Plan.md), Abschnitte 0, 6 und 7.
Eine Seite entsteht jeweils mit dem Werkzeug.

| Werkzeug | Aufgabe | ersetzt | Welle |
|---|---|---|---|
| `gannot` | Annotation | Ensembl VEP | 2 |
| `gqc` | Kontamination, Identität, Geschlecht | VerifyBamID2, Somalier, mosdepth | 2 |
| `ghap`, `glift`, `gmicro`, `gtel`, `gpgs` | Haplogruppen, Liftover, Viren, Telomere, polygene Scores | Haplogrep, Kraken2, TelomereHunter2 u. a. | 2 |
| `ghla` | HLA (pjlabshla, fertig, Einbau offen) | HLA*LA, T1K | 2 |
| `gpgx`, `gkir`, `gstr`, `gloci` | Pharmakogenetik, KIR, Repeats, Spezialgene (Lp(a), SMN u. a.) | Aldy, T1K, ExpansionHunter, KILDA | 3 |
| `gcall`, `gcnv`, `gsv` | kleine Varianten, Kopienzahl, strukturelle Varianten | bcftools, GATK, CNVkit, Delly, Manta | 4 und 5 |
