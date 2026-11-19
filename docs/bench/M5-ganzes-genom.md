# M5: `genome fq2bam` auf dem ganzen Genom

Stand: 19.11.2026. Werkzeug: [`genome fq2bam`](../handbuch/werkzeuge/fq2bam.md).

## Fragen
1. Läuft der eine Befehl auf dem ganzen Referenzgenom durch, und landen die Reads am wahren Ort?
2. Gilt die Gleichheit mit den einzelnen Schritten und mit der Kette aus samtools auch hier?
3. Was braucht der Lauf an Zeit, Speicher und Platz für die Teile?

## Aufbau
- Referenz GRCh38 ohne ALT-Contigs (195 Contigs, 3.099.922.541 Basen), Index
  `grch38.kompakt.gcfm` (4,84 GB, blockweise gebaut in 16,6 Minuten).
- Reads: 2.000.000 erzeugte Paare zu 150 Basen über alle Contigs, gewichtet nach Länge, Stücke
  von 400 ± 50 Basen, 0,2 % Substitutionen, keine Indels, feste Qualität
  (`bench/m5/simulate_genome_pairs.py`, Saat 7). Das ist eine 0,2-fache Abdeckung. Der wahre Ort
  steht im Namen.
- Wahrheit: `bench/m5/check_truth.py` (Contig gleich und linkes Ende höchstens 5 Basen daneben,
  weiche Enden zurückgerechnet), Ergebnis `bench/m5/ganzes_genom_wahrheit.json`.
- Vergleich: `bench/m5/compare_fq2bam.py` wie in [M5 fq2bam](M5-fq2bam.md), Prüfer samtools 1.24
  aus pysam 0.24.1, Ergebnis `bench/m5/ganzes_genom_vergleich.json`.
- Rechner: Cloud-Maschine mit 4 virtuellen Kernen und 15 GB RAM, geteilt mit anderen. Je
  Einstellung ein Lauf. Zeit, CPU-Zeit und Speicherspitze des Hauptlaufs (512 MB Sortierspeicher)
  über `getrusage`, Ergebnis `bench/m5/ganzes_genom_lauf.json`.

## Ergebnis

### 1. Orte gegen die Wahrheit

| MAPQ | Reads | am wahren Ort | Anteil |
|---|---|---|---|
| 30 und mehr | 3.862.782 | 3.862.773 | 99,9998 % |
| 1 bis 29 | 42.290 | 41.772 | 98,78 % |
| 0 | 94.928 | 31.065 | 32,73 % |
| alle | 4.000.000 | 3.935.610 | 98,39 % |

Alle Reads haben einen Ort. Die falschen liegen fast alle bei MAPQ 0, also dort, wo das Genom
den Ort nicht eindeutig hergibt und das Werkzeug das auch sagt. Die Ausgabe ist sortiert,
`samtools quickcheck` nimmt sie an.

### 2. Gleichheit

| Lauf | Teile auf der Platte | Zeilen gleich wie `index sam`, `sort`, `markdup` | Zeit |
|---|---|---|---|
| 4 Kerne, 2048 MB | 0 | alle 4.000.000 | 347 s |
| 1 Kern, 2048 MB | 0 | alle 4.000.000 | 1.457 s |
| 4 Kerne, 8 MB | 164 | alle 4.000.000 | 362 s |

Gegen `samtools fixmate -m`, `sort`, `markdup`: alle 4.000.000 Zeilen gleich, wenn TLEN und die
Reihenfolge der Zusatzangaben nicht zählen. TLEN unterscheidet sich in 6 Zeilen. Markiert sind in
beiden dieselben 4 Datensätze. `idxstats`, `flagstat` und die Statistik der Duplikate stimmen.

### 3. Zeit, Speicher, Platz (4 Kerne, 512 MB Sortierspeicher)

| Größe | Wert |
|---|---|
| Laufzeit | 348 s, 11.500 Reads je Sekunde |
| CPU-Zeit | 1.322 s |
| Speicherspitze | 5,24 GB |
| Index laden | 3,7 s |
| Ausrichten und Teile schreiben | 336,7 s (97 %) |
| Zusammenführen, Duplikate markieren, BAM schreiben | 6,1 s |
| BAM-Index | 0,9 s |
| Datensätze roh | 1,38 GB |
| Teile auf der Platte | 2 Teile, 0,32 GB |
| BAM | 0,29 GB |

Mit 4 Kernen ist der Lauf 4,2-mal so schnell wie mit einem.

## Grenzen
- Erzeugte Reads ohne Indels und mit fester Qualität. Echte Reads sind schwerer auszurichten,
  und echte Qualitäten lassen sich schlechter packen: Teile und BAM werden je Read größer.
- 0,2-fache Abdeckung. Bei 4 gefundenen Duplikaten sagt die Messung nichts über das Markieren
  bei voller Abdeckung, dafür steht [M5 Duplikate](M5-duplikate.md) auf chr20.
- Kein Lauf mit 30-facher Abdeckung. Rechnet man die 11.500 Reads je Sekunde gleichmäßig hoch,
  bräuchte ein Genom mit 650 Millionen Reads auf diesen 4 Kernen rund 16 Stunden und etwa 70 GB
  für die Teile, mit echten Qualitäten eher mehr. Das ist eine Schätzung und keine Messung. Ein Vergleich mit den 11 h 40 min von
  Genome 3.x auf 10 Kernen ist daraus nicht ableitbar, dafür braucht es denselben Rechner und
  dieselben Reads.
- Referenz ohne ALT-Contigs. Die ALT-Behandlung von bwa (`.alt`-Datei) ist nicht geprüft.

## Wiederholen
```bash
genome index build grch38.fa.gz grch38.kompakt.gcfm --verfahren blockweise
python bench/m5/simulate_genome_pairs.py grch38.fa.gz 2000000 150 400 50 0.002 7 w_1.fq w_2.fq
genome fq2bam grch38.kompakt.gcfm --fastq w_1.fq --fastq2 w_2.fq --aus w.bam --speicher-mb 512
python bench/m5/check_truth.py w.bam
python bench/m5/compare_fq2bam.py target/release/genome grch38.kompakt.gcfm \
  w_1.fq w_2.fq arbeit bench/m5/ganzes_genom_vergleich.json
```
