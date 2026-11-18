# M3: ganze SAM-Zeilen gegen bwa

Stand: 18.11.2026. Verfahren: [Ausgabe als SAM](../handbuch/verfahren/sam-ausgabe.md).

## Aufbau
- Referenz chr20 aus GRCh38; Satz A (200.000 Reads, 1 % Austauschfehler) und Satz B (dazu Indels
  bis 10 Basen), wie in [M3 CIGAR](M3-cigar.md). Einzelne Reads, keine Paare.
- Genome 4.0: `genome index sam chr20.kompakt.gcfm --fastq reads.fq > unser.sam`.
- bwa: bwa mem aus pybwa 2.3.2 mit Vorgaben. Vergleich mit `bench/m3/compare_bwa_sam.py`: je
  Read alle Zeilen, Feld für Feld (Flag, Contig, Position, MAPQ, CIGAR, Basen, Qualitäten, `NM`,
  `MD`, `AS`, `XS`, `SA`, `XA`), Zuordnung nach Reihenfolge.

## Ergebnis

| | Satz A | Satz B |
|---|---:|---:|
| Reads | 200.000 | 200.000 |
| alle Zeilen in allen Feldern gleich | 199.964 | 199.956 |
| andere Zahl von Zeilen | 0 | 0 |
| Reads mit `XA` bei bwa | 5.281 | 5.366 |
| Reads mit ergänzenden Zeilen bei bwa | 0 | 0 |
| Reads ohne Ort bei bwa | 0 | 0 |

Abweichende Reads nach den Feldern, die sich unterscheiden:

| Felder | Satz A | Satz B |
|---|---:|---:|
| XS | 25 | 27 |
| pos | 7 | 8 |
| XA | 1 | 4 |
| XA, XS, mapq | 2 | 1 |
| XS, mapq | 1 | 1 |
| mapq | 0 | 1 |
| AS, MD, NM, XS, cigar, pos | 0 | 1 |
| AS, MD, NM, cigar, pos | 0 | 1 |

## Einordnung
- **Nur `XS`** (25 und 27 Reads), auch **`XS` mit MAPQ**: bwa findet einen Seed von 19 oder 20
  Basen in einer mit Zufallsbasen gefüllten N-Strecke; Genome 4.0 verwirft solche Fundstellen
  bewusst. An einem Read nachvollzogen (`chr20:13337599:-:1`: bwa meldet einen Bereich von 20
  Basen bei chr20:64.334.606, dort stehen in der Referenz nur N), bei den übrigen nach dem
  gleichen Muster angenommen, nicht einzeln geprüft.
- **Position**, und die Zeilen mit CIGAR: die Reads am Zentromer aus
  [M3 CIGAR](M3-cigar.md), nicht geklärt.
- **Nur `XA`**: nicht geklärt.
- Im ersten Lauf wich `XA` bei 322 Reads (Satz A) ab. Ursache war die Schwelle von 0,8, die bwa
  mit einfacher Genauigkeit hält; siehe Verfahrensseite.

## Was der Vergleich nicht zeigt
- Kein Read in beiden Sätzen ist geteilt oder ohne Ort. Ergänzende Zeilen, `SA` und Reads ohne
  Ort sind nur durch eigene Tests belegt.
- Verglichen mit dem in pybwa eingebauten bwa, nicht mit dem Programm bwa 0.7.19.
- Tempo: 81 s (Satz A) und 86 s (Satz B) für 200.000 Reads auf einem Kern, nichts optimiert.

Rohdaten: `bench/m3/ergebnisse/sam_vergleich.json`, `sam_vergleich_indel.json`.
