# M3: erster Vergleich mit bwa (Ort, Strang, MAPQ, Score)

Stand: 18.11.2026. Verfahren: [Bereiche bewerten](../handbuch/verfahren/bereiche-bewerten.md) und
die Seiten davor. Verglichen wird je Read der primäre Bereich von Genome 4.0 mit dem primären
Alignment von bwa mem. CIGAR und Ausgabe als BAM gibt es noch nicht; sie sind nicht Teil dieses
Vergleichs.

## Aufbau
- Referenz: chr20 aus GRCh38. Reads: 200.000 erzeugte zu 150 Basen, 1 % Austauschfehler, beide
  Stränge (`bench/m3/simulate_reads.py`, Saat 1), einzeln ausgerichtet (kein Paar).
- Genome 4.0: `genome index ketten chr20.kompakt.gcfm --fastq sim150.fq --tabelle unser.tsv`.
- bwa: das in pybwa 2.3.2 eingebaute bwa mem mit Vorgaben, über `bench/m3/compare_bwa_align.py`.
  Das ist nicht das Programm `bwa 0.7.19` von pj's Mac; welche bwa-Version pybwa einbettet, ist
  hier nicht geprüft. Der Vergleich mit genau 0.7.19 steht aus.
- Verglichen: Contig, Position des linken Endes, Strang, MAPQ, und der Alignment-Score (bwa `AS`
  gegen unseren besten lokalen Score).

## Ergebnis

| | Reads | Anteil |
|---|---|---|
| Ort und Strang gleich | 199.992 | 99,996 % |
| davon auch Score gleich | 199.990 | |
| davon auch MAPQ gleich | 199.926 | 99,963 % aller Reads |
| Ort gleich, MAPQ verschieden | 66 | 0,033 % |
| anderer Ort, bei beiden MAPQ 0 | 7 | |
| anderer Ort, sonst | 1 | |
| nur bei einem der beiden ohne Alignment | 0 | |

Verteilung unserer MAPQ: 186.390 Reads mit 60, 3.027 mit 30 bis 59, 2.550 mit 1 bis 29, 8.033
mit 0.

## Die Abweichungen
- **7 Reads, anderer Ort, beide MAPQ 0:** gleich gute Orte; welcher gemeldet wird, ist bei beiden
  Programmen Zufall nach einem Streuwert. Bei den übrigen rund 6.900 Reads mit gleich guten Orten
  fällt die Wahl gleich aus, weil der Streuwert von bwa nachgebildet ist.
- **1 Read, anderer Ort:** `chr20:63035123:+:3`. Genome 4.0 meldet den wahren Ort mit 150 Basen
  ohne Lücke, Score 135, MAPQ 60. bwa meldet 36 Basen weiter links ein Alignment mit einer
  Deletion von 36 Basen (`57M36D93M`), Score 98, MAPQ 16. Ursache ist sehr wahrscheinlich das in
  Genome 4.0 noch fehlende Zusammenfügen zweier Bereiche (`mem_patch_reg`), das bwa hier anwendet;
  geprüft ist das nicht. Das ist ein Einzelfall und kein Beleg, dass Genome 4.0 genauer ist.
- **66 Reads, gleiche Lage, andere MAPQ:** Abweichungen in beide Richtungen, von -40 bis +45. Die
  Ursache ist nicht untersucht. Naheliegend sind Unterschiede beim zweitbesten Bereich (`sub`,
  `sub_n`), etwa durch das fehlende Zusammenfügen oder durch die Reihenfolge bei Gleichstand.
- **2 Reads, gleiche Lage, anderer Score:** nicht untersucht.

## Der Fund zur MAPQ
Im ersten Lauf wich die MAPQ bei 5.762 Reads ab, immer nach oben. Ursache: bwa hält einen Faktor
seiner MAPQ-Formel, gedacht als ln(50) = 3,912, in einer Ganzzahl und rechnet mit 3. Nach dem
Nachbilden sind es 66. Einzelheiten in
[Bereiche bewerten](../handbuch/verfahren/bereiche-bewerten.md).

## Gegen die Wahrheit
Aus demselben Lauf (siehe [M3 Bereiche](M3-bereiche.md)): Der beste Bereich liegt bei 99,67 % der
Reads am wahren Ort. Da Ort und Strang bei 199.992 von 200.000 Reads mit bwa übereinstimmen,
liegt bwa auf diesen Daten praktisch gleich. Eine Auswertung "falsch bei MAPQ ab 30" für beide
Programme fehlt noch.

## Was der Vergleich nicht zeigt
- Nur chr20, nur erzeugte Reads mit Austauschfehlern, keine Indels, keine echten Reads, keine
  Paare. Indels und geteilte Reads sind genau die Fälle, in denen das fehlende Zusammenfügen
  und der fehlende CIGAR zählen.
- Kein Tempovergleich: bwa brauchte hier 57 s, Genome 4.0 74 s, beide auf einem Kern; Genome 4.0
  rechnet aber noch keinen CIGAR und schreibt keine Ausgabe, und nichts ist optimiert.

Rohdaten: `bench/m3/ergebnisse/vergleich_bwa.json`, `mapq.json`.
