# M5: `genome fq2bam` als Ganzes

Stand: 19.11.2026. Werkzeug: [`genome fq2bam`](../handbuch/werkzeuge/fq2bam.md).

## Fragen
1. Gibt der eine Befehl dasselbe wie die einzelnen Schritte, unabhängig von Kernen und Speicher?
2. Stimmt das Ergebnis mit der Kette aus samtools überein?
3. Sind Index und Zahlen im Bericht richtig?

## Aufbau
- Referenz chr20 (GRCh38), Index `chr20.schnell.gcfm`. Reads: die 130.268 Paare mit 30.268
  erzeugten Kopien aus [M5 Duplikate](M5-duplikate.md).
- Prüfer: samtools 1.24 aus pysam 0.24.1.
- `bench/m5/compare_fq2bam.py`, Ergebnis `bench/m5/fq2bam_vergleich.json`.
- Rechner: Cloud-Maschine mit 4 virtuellen Kernen, geteilt mit anderen. Je Einstellung ein Lauf.

## Ergebnis

### 1. Gegen die einzelnen Schritte

| Lauf | Teile auf der Platte | Zeilen gleich wie `index sam`, `sort`, `markdup` | Zeit |
|---|---|---|---|
| 4 Kerne, 2048 MB | 0 | alle 260.536 | 24,4 s |
| 1 Kern, 2048 MB | 0 | alle 260.536 | 96,8 s |
| 4 Kerne, 8 MB | 10 | alle 260.536 | 25,2 s |

### 2. Gegen die Kette mit samtools
`genome index sam`, dann `samtools fixmate -m`, `sort`, `markdup`: alle 260.536 Zeilen gleich,
wenn TLEN und die Reihenfolge der Zusatzangaben nicht zählen. TLEN unterscheidet sich in 1.318
Zeilen (Paare, die voneinander weg zeigen oder auf demselben Strang liegen; dort bleibt der Wert
von bwa). Markiert sind in beiden dieselben 56.067 Datensätze.

### 3. Index und Zahlen
- `samtools quickcheck` bestanden.
- `samtools idxstats` mit dem Index von fq2bam gleich wie mit einem Index von samtools.
- Zahlen im Bericht gleich wie `samtools flagstat` (sieben verglichene Zeilen) und wie die
  Statistik von `samtools markdup` (sechs Zahlen, auch die geschätzte Bibliotheksgröße).

### Wo die Zeit bleibt (4 Kerne, 24,4 s)

| Schritt | Zeit |
|---|---|
| Index laden | 0,2 s |
| Ausrichten und Teile schreiben | 23,4 s |
| Zusammenführen, Duplikate markieren, BAM schreiben | 0,7 s |
| BAM-Index | 0,1 s |

Fast die ganze Zeit ist Ausrichten. Alles danach macht hier 3 % aus.

## Grenzen
- Ein Chromosom, 0,6-fache Abdeckung, erzeugte Reads. Die Zeiten sagen nichts über ein ganzes
  Genom: Dort ist der Index 50-mal größer, Reads aus Wiederholungen sind teurer, und das Sortieren
  schreibt Teile auf die Platte. Ein Vergleich mit den 11 h 40 min von Genome 3.x ist daraus nicht
  ableitbar.
- Kein Zeitvergleich mit bwa und samtools unter gleichen Bedingungen in dieser Messung.

## Wiederholen
```bash
python bench/m5/compare_fq2bam.py target/release/genome chr20.schnell.gcfm \
  d_1.fq d_2.fq arbeit bench/m5/fq2bam_vergleich.json
```
