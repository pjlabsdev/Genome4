# M5: Sortieren nach Ort

Stand: 19.11.2026. Verfahren: [Sortieren nach Ort mit festem Speicher](../handbuch/verfahren/sortieren.md).
Werkzeug: [`genome sort`](../handbuch/werkzeuge/sort.md).

## Fragen
1. Ist die Reihenfolge Datensatz für Datensatz die von `samtools sort`?
2. Geht dabei etwas verloren oder wird etwas verändert?
3. Hängt das Ergebnis vom Speicher ab (mit und ohne Teile auf der Platte)?

## Aufbau
- Prüfer: `samtools sort` 1.24 aus pysam 0.24.1. Verglichen wird der Text von `samtools view`,
  also alle Felder jedes Datensatzes, in der Reihenfolge der Datei.
- Eingaben:
  - `wie_ausgerichtet`: die 100.000 Paare aus [M3 Paare](M3-paare.md) auf chr20, als BAM aus
    `genome index sam --bam` mit Lesegruppe (200.000 Datensätze, 68 MB unkomprimiert, 972 ohne
    eigenen Ort).
  - `gemischt`: dieselben Datensätze in zufälliger Reihenfolge.
  - `kuenstlich_7_contigs`: 300.000 erzeugte Datensätze auf 7 Contigs, Positionen 0 bis 199 (also
    im Mittel über 200 Datensätze je Ort), beide Stränge, 5 % ohne Contig, 5 % ohne eigenen Ort
    mit Contig und Position.
- Je Eingabe drei Läufe mit `--speicher-mb` 2048, 16 und 2, `--kerne 1`.
- Skript: `bench/m5/compare_sort.py`, Ergebnis `bench/m5/sort_vergleich.json`.
- Rechner: Cloud-Maschine mit 4 virtuellen Kernen, geteilt mit anderen. Je Einstellung ein Lauf.

## Ergebnis

| Eingabe | Datensätze | Speicher (MB) | Teile | Folge wie samtools | Menge wie Eingabe | `quickcheck` |
|---|---|---|---|---|---|---|
| wie ausgerichtet | 200.000 | 2048 | 0 | ja | ja | bestanden |
| wie ausgerichtet | 200.000 | 16 | 4 | ja | ja | bestanden |
| wie ausgerichtet | 200.000 | 2 | 32 | ja | ja | bestanden |
| gemischt | 200.000 | 2048 | 0 | ja | ja | bestanden |
| gemischt | 200.000 | 16 | 4 | ja | ja | bestanden |
| gemischt | 200.000 | 2 | 32 | ja | ja | bestanden |
| künstlich, 7 Contigs | 300.000 | 2048 | 0 | ja | ja | bestanden |
| künstlich, 7 Contigs | 300.000 | 16 | 1 | ja | ja | bestanden |
| künstlich, 7 Contigs | 300.000 | 2 | 8 | ja | ja | bestanden |

In allen neun Läufen steht `SO:coordinate` im Kopf, und die übrigen Kopfzeilen sind unverändert.
Die Prüfsumme über alle Zeilen ist je Eingabe für alle drei Speichergrößen dieselbe.

## Zeiten
Nur zur Einordnung, nicht als Vergleich: Die Dateien sind klein, es gab je einen Lauf, und samtools
lief innerhalb von Python.

| Eingabe | `genome sort` 2048 / 16 / 2 MB | `samtools sort` |
|---|---|---|
| wie ausgerichtet | 0,86 / 1,15 / 1,15 s | 1,92 s |
| gemischt | 0,84 / 1,23 / 1,20 s | 1,97 s |
| künstlich | 0,19 / 0,58 / 0,54 s | 0,30 s |

## Grenzen
- Kleine Dateien, ein Chromosom bei den echten Daten. Ein ganzes Genom (Laufzeit, Speicherspitze,
  Platz für die Teile) ist nicht gemessen.
- Kein Tempo-Vergleich mit samtools unter gleichen Bedingungen.
- Zusätzliche Alignments (Flags 256 und 2048) kommen in den Daten nicht vor; für die Reihenfolge
  spielen sie keine Rolle, weil nur Contig, Position und Strang zählen.

## Wiederholen
```bash
genome index sam chr20.schnell.gcfm --fastq pairs_1.fq --fastq2 pairs_2.fq \
  --bam paare.bam --lesegruppe '@RG\tID:l1\tSM:p'
python bench/m5/compare_sort.py target/release/genome paare.bam arbeit bench/m5/sort_vergleich.json
```
