# M5: `genome fq2bam` mit echten Reads

Stand: 19.11.2026. Werkzeug: [`genome fq2bam`](../handbuch/werkzeuge/fq2bam.md).

## Fragen
1. Schreibt Genome 4.0 für echte Reads dieselben Zeilen wie bwa mem?
2. Gibt `genome fq2bam` dasselbe wie die Kette mit samtools?
3. Wie schnell ist es mit echten Reads?

## Aufbau
- Reads: GIAB HG002 (Illumina HiSeq, 2 × 148 Basen), der Ausschnitt chr20:10.000.000 bis
  12.000.000 aus `HG002.GRCh38.60x.1.bam` (wie in [M1](M1-bam-rundreise.md)). Daraus alle Paare,
  deren beide Reads im Ausschnitt stehen, zurück in die Richtung des Geräts gedreht:
  421.564 Paare (`bench/m5/bam_to_pairs.py`).
- Referenz: chr20 aus GRCh38, Index `chr20.kompakt.gcfm`.
- bwa: das in pybwa 2.3.2 eingebaute bwa mem 0.7.18 mit Vorgaben, Vergleich Zeile für Zeile und
  Feld für Feld mit `bench/m3/compare_bwa_pairs.py`, Ergebnis `bench/m5/echte_reads_bwa.json`.
- samtools 1.24 aus pysam 0.24.1, Vergleich mit `bench/m5/compare_fq2bam.py`, Ergebnis
  `bench/m5/echte_reads_fq2bam.json`.
- Rechner: Cloud-Maschine mit 4 virtuellen Kernen, geteilt mit anderen. Je Einstellung ein Lauf.

## Ergebnis

### 1. Gegen bwa mem

| | Zahl |
|---|---:|
| Reads | 843.128 |
| Zeilen (mit 878 ergänzenden) | 844.006 |
| Reads, bei denen alle Zeilen in allen Feldern gleich sind | 843.109 (99,998 %) |
| Paare ganz gleich | 421.545 |
| Reads mit Abweichung | 19, alle nur in `XS` |

Die 19 Reads haben bei bwa `XS` 19 bis 22 und bei uns 0. Bei allen 19 passt ein Stück von
19 Basen des Reads zufällig auf Zufallsbasen, mit denen bwa eine N-Strecke der Referenz füllt
(`bench/m5/check_n_hits.py`: 19 von 19). Genome 4.0 legt dort bewusst nichts an
([SAM-Ausgabe](../handbuch/verfahren/sam-ausgabe.md), Feinheiten). Ort, Strang, MAPQ, CIGAR und
alle anderen Angaben dieser Reads sind gleich.

Zwei Dinge bleiben aus dem Vergleich heraus, weil pybwa seine Datensätze selbst baut und dabei
vom Programm bwa abweicht (geprüft am Quelltext von bwa, `mem_aln2sam`, und am Patch 0006 von
pybwa):
- Flag 32 bei einem Partner ohne Ort (1.127 Zeilen): bwa setzt es nach dem Strang des Reads,
  pybwa nicht. Genome 4.0 hält sich an bwa.
- Basen und Qualitäten hart gekürzter Zeilen auf dem Gegenstrang (474 Zeilen): pybwa schneidet
  am falschen Ende. Bei Genome 4.0 passen die Basen zur CIGAR und zu `MD`.

### 2. Zwei Fehler, die erst die echten Reads zeigten

**Doppelte Seeds.** Vor der Korrektur wichen 38 Reads ab, 19 davon in CIGAR, `SA`, `XA` oder
MAPQ. Sie liegen in Wiederholungen, die meisten an einer Folge von 60 A bei chr20:34.228.133.
Dort liefern mehrere Runden der Seed-Suche denselben Treffer. bwa behält solche Doppelten, und
sie entscheiden mit, welche Ketten entstehen und welche von gleich guten Lagen gewinnt. Genome
4.0 hatte sie entfernt. Jetzt bleiben sie stehen
([Seeds und Ketten](../handbuch/verfahren/seeds-und-ketten.md)); die Ketten dieses Beispiels
sind danach dieselben wie in bwa (geprüft gegen `mem_chain` und `mem_chain_flt` aus der
Bibliothek). Bei den erzeugten Reads aus [M3 Paare](M3-paare.md) bleibt die Übereinstimmung mit bwa dabei gleich.

**Partner-Angaben auf ergänzenden Zeilen.** Gegen die Kette mit samtools wichen 40 Zeilen ab,
alle ergänzend (Flag 2048). `samtools fixmate` gibt ihnen Ort, Strang, `MQ` und `MC` der
Hauptzeile des Partners. bwa schreibt dort `MC` mit hartem statt weichem Kürzen und, wenn der
Partner keinen Ort hat, den Ort der eigenen Zeile. `--partner-angaben` und `fq2bam` machen es
jetzt wie fixmate.

### 3. Gegen die Kette mit samtools (nach den Korrekturen)

| Lauf | Teile auf der Platte | Zeilen gleich wie `index sam`, `sort`, `markdup` | Zeit |
|---|---|---|---|
| 4 Kerne, 2048 MB | 0 | alle 844.006 | 20,0 s |
| 1 Kern, 2048 MB | 0 | alle 844.006 | 68,1 s |
| 4 Kerne, 8 MB | 34 | alle 844.006 | 20,7 s |

Gegen `samtools fixmate -m`, `sort`, `markdup`: alle 844.006 Zeilen gleich, wenn TLEN und die
Reihenfolge der Zusatzangaben nicht zählen. TLEN unterscheidet sich in 764 Zeilen (wie in
[M5 Duplikate](M5-duplikate.md) beschrieben). Markiert sind in beiden dieselben 7.967
Datensätze. `quickcheck`, `idxstats`, `flagstat` und die Statistik der Duplikate stimmen.

### 4. Tempo
- `fq2bam` auf einem Kern: 68,1 s für 843.128 Reads, also 12.400 Reads je Sekunde, alle Schritte
  eingeschlossen. Auf 4 Kernen 20,0 s, davon 17,9 s Ausrichten.
- bwa meldet für dieselben Reads 50,8 s CPU-Zeit auf einem Kern, nur für das Ausrichten
  (16.600 Reads je Sekunde). Mit dem kompakten Index brauchen wir damit je Kern rund 1,3-mal so
  lange wie bwa, wie bei den erzeugten Reads in [M3 Tempo](M3-tempo.md).
- Echte Reads laufen hier rund 5-mal so schnell wie die erzeugten aus
  [M5 fq2bam](M5-fq2bam.md) (2.700 Reads je Sekunde und Kern). Das gilt für bwa genauso.

## Grenzen
- Ein Ausschnitt von 2 Mb eines Chromosoms, Referenz nur chr20. Reads, die bei einer Suche im
  ganzen Genom anderswo landen würden, haben hier keine Wahl.
- Die Reads stammen aus einer schon ausgerichteten Datei: Paare, deren Partner außerhalb des
  Ausschnitts liegt, fehlen. Schwierige Paare sind damit seltener als in einer rohen FASTQ.
- Verglichen ist mit bwa 0.7.18 aus pybwa und samtools 1.24, nicht mit den Programmen bwa 0.7.19
  und samtools 1.23.1 von Genome 3.x. Ob samtools 1.23.1 die ergänzenden Zeilen schon so
  behandelt wie 1.24, ist nicht geprüft.
- Zeiten: je ein Lauf auf einer geteilten Maschine.

## Wiederholen
```bash
python bench/m5/bam_to_pairs.py hg002_chr20.bam e_1.fq e_2.fq
genome index sam chr20.kompakt.gcfm --fastq e_1.fq --fastq2 e_2.fq > wir.sam
python bench/m3/compare_bwa_pairs.py chr20.fa e_1.fq e_2.fq wir.sam bench/m5/echte_reads_bwa.json
python bench/m5/compare_fq2bam.py target/release/genome chr20.kompakt.gcfm \
  e_1.fq e_2.fq arbeit bench/m5/echte_reads_fq2bam.json
```
