# M4: Ausrichtung gegen die Wahrheit, mit bwa auf denselben Reads

Stand: 25.11.2026. Werkzeuge: [`genome fq2bam`](../handbuch/werkzeuge/fq2bam.md) und das
Prüfprogramm [`genome-validate`](../handbuch/werkzeuge/genome-validate.md). Verfahren:
[Wahrheit erzeugen und Ausrichtung prüfen](../handbuch/verfahren/wahrheit-simulieren.md).

## Fragen
1. Landen die Reads am wahren Ort, auch mit Varianten der Probe, Indels und echten
   Qualitätsstufen?
2. Ist Genome 4 dabei mindestens so gut wie bwa, Read für Read auf denselben Reads?
3. Stimmt die MAPQ mit der tatsächlichen Fehlerrate überein?

Das Messziel aus dem Plan (Abschnitt 0.4): falsch platzierte Reads mit MAPQ ab 20 höchstens so
viele wie bei bwa, MAPQ mindestens so gut kalibriert.

## 1. Probe in der Cloud: chr20

### Aufbau
- Referenz: chr20 aus GRCh38 (64.444.167 Basen), Index `chr20.kompakt.gcfm`.
- Reads: `genome-validate reads-erzeugen --referenz chr20.fa --ausgabe c20 --paare 200000
  --saat 3` (Profil Illumina, Stücke 400 ± 60, 150 Basen). Die Probe hat 71.620 Varianten
  (63.701 SNVs, 7.919 Indels bis 10 Basen); 1.638 Stücke mit N wurden verworfen und ersetzt.
  Einstellungen: `bench/m4/ergebnisse/chr20_reads.json`.
- Genome 4: `genome fq2bam chr20.kompakt.gcfm --fastq c20_1.fq.gz --fastq2 c20_2.fq.gz
  --aus c20.g4.bam --ohne-bericht`, Stand main 02fdd93 plus diesem Zweig.
- bwa: das in pybwa 2.3.2 eingebaute bwa mem mit Vorgaben und Blöcken von 10 Mio. Basen
  (`bench/m4/pybwa_bam.py`). Das ist nicht `bwa 0.7.19`; der Vergleich mit genau dieser Version
  folgt auf dem ganzen Genom am Mac (Abschnitt 2).
- Prüfung: `genome-validate ausrichtung-pruefen --contigs c20.contigs.tsv --bam genome4=c20.g4.bam
  --bam bwa=c20.bwa.bam --json chr20_wahrheit.json --falsche chr20_sicher_falsch.tsv`, Toleranz
  20 Basen. Ergebnis: `bench/m4/ergebnisse/chr20_wahrheit.json`.

### Ergebnis

| | Genome 4 | bwa |
|---|---|---|
| Reads | 400.000 | 400.000 |
| ohne Ort | 0 | 0 |
| nicht am wahren Ort | 8.545 (2,14 %) | 8.535 (2,13 %) |
| Reads mit MAPQ ab 20 | 386.865 | 386.865 |
| davon falsch | 13 (0,0034 %) | 13 (0,0034 %) |

Je MAPQ-Stufe (falsch / erwartet heißt: tatsächlich falsche Reads und die Zahl, die die MAPQ
verspricht):

| MAPQ | Reads Genome 4 | falsch / erwartet | Reads bwa | falsch / erwartet |
|---|---|---|---|---|
| 0 | 10.218 | 8.174 / 10.218 | 10.218 | 8.164 / 10.218 |
| 1 bis 9 | 1.508 | 310 / 495,2 | 1.506 | 310 / 494,5 |
| 10 bis 19 | 1.409 | 48 / 64,8 | 1.411 | 48 / 64,9 |
| 20 bis 29 | 1.251 | 7 / 5,5 | 1.249 | 7 / 5,5 |
| 30 bis 39 | 462 | 1 / 0,2 | 462 | 1 / 0,2 |
| 40 bis 49 | 3.892 | 4 / 0,3 | 3.894 | 4 / 0,3 |
| 50 bis 59 | 742 | 0 / 0,0 | 742 | 0 / 0,0 |
| 60 | 380.518 | 1 / 0,4 | 380.518 | 1 / 0,4 |

Read für Read: beide richtig 391.451, nur Genome 4 richtig 4, nur bwa richtig 14, beide nicht
8.531. Sicher (MAPQ ab 20) falsch nur in einem der beiden: 0. Die 13 sicher falschen Reads sind
in beiden dieselben (`chr20_sicher_falsch.tsv`).

### Was das heißt
- Auf chr20 ist Genome 4 so genau wie bwa: dieselben Reads sicher falsch, dieselben Stufen.
- Die 18 Reads, die nur einer der beiden richtig hat, haben bei beiden MAPQ 0: Der Ort ist
  gleich gut, die Wahl zwischen gleich guten Orten fällt anders. Das bleibt zu erklären (pybwa
  zählt die Reads für diese Wahl über alle Blöcke anders als unser Lauf).
- Die MAPQ ist bei beiden zwischen 20 und 49 zu zuversichtlich: 12 falsche statt 6 erwarteter.
  Hier kann Genome 4 besser werden als bwa. 10 der 13 sicher falschen Reads sind fünf ganze
  Paare, die gemeinsam an einer Kopie desselben Abschnitts 410 kb bis 920 kb entfernt liegen; 3
  liegen 26 bis 64 Basen neben dem wahren Ort.
- Bei MAPQ 0 liegen 80 % falsch. Das ist richtig so: MAPQ 0 sagt, dass der Ort nicht eindeutig
  ist.

## 2. Ganzes Genom am Mac gegen bwa 0.7.19

Steht aus. Geplant: dieselbe Messung auf hs38d1 (wie pj's Daten), je 1 Mio. Paare mit Profil
Illumina und MGI, Genome 4 gegen `bwa mem` 0.7.19 aus Genome 3, ausgewertet zusätzlich in den
schweren Bereichen der GIAB-Schichten (Segmentduplikationen, schwer abbildbare Bereiche).

## Grenzen
- Siehe [`genome-validate`](../handbuch/werkzeuge/genome-validate.md): Fehlermodell nicht aus
  echten Läufen gelernt, Probe ohne Strukturvarianten. Die Zahlen vergleichen die beiden
  Ausrichter auf denselben Reads; als Fehlerrate für echte Daten taugen sie nicht.
- chr20 allein: Reads aus anderen Chromosomen, die sich auf chr20 verirren könnten, fehlen. Das
  ganze Genom misst Abschnitt 2.
- Ein Lauf mit einer Saat. Bei 13 sicher falschen Reads sind Unterschiede von wenigen Reads
  Zufall.
