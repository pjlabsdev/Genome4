# M3: Read-Paare gegen bwa und gegen den wahren Ort

Stand: 18.11.2026. Verfahren: [Read-Paare](../handbuch/verfahren/read-paare.md).

## Aufbau
- Referenz chr20 aus GRCh38.
- 100.000 erzeugte Paare, je 150 Basen, Stücklänge normalverteilt (Mittel 400, Streuung 60),
  1 % Austauschfehler, Indels bis 10 Basen mit Rate 0,0005 je Base. Bei 1 % der Paare ist Read 2
  eine Zufallsfolge (findet keinen Ort), bei 1 % stammt Read 2 aus einem anderen Stück.
  `bench/m3/simulate_pairs.py ... 100000 150 400 60 0.01 0.0005 0.01 0.01 3`. Kein
  realistisches Fehlermodell.
- Genome 4.0: `genome index sam chr20.kompakt.gcfm --fastq r1.fq --fastq2 r2.fq > unser.sam`,
  Blöcke von 10 Millionen Basen (33.334 Paare), also drei Blöcke.
- bwa: bwa mem aus pybwa 2.3.2 (bwa 0.7.18) mit Vorgaben, Paare, dieselben Blöcke. Vergleich mit
  `bench/m3/compare_bwa_pairs.py`: je Read alle Zeilen, Feld für Feld (Flag, Contig, Position,
  MAPQ, CIGAR, Contig und Position des Partners, TLEN, Basen, Qualitäten, `NM`, `MD`, `MC`, `MQ`,
  `AS`, `XS`, `SA`, `XA`), Zuordnung nach Reihenfolge.
- Die Ausgabe von pybwa trägt die Flag-Bits 2 (richtiges Paar), 64 und 128 (Read 1, Read 2)
  nicht. Sie sind aus dem Vergleich herausgenommen.

## Insertgröße
Beide Seiten schätzen je Block dasselbe (Richtung FR, die anderen drei Richtungen haben keine
Paare):

| Block | Paare in der Schätzung | Perzentile 25/50/75 | Mittel | Streuung | richtiges Paar von bis |
|---|---:|---|---:|---:|---|
| 0 | 30.030 | 359 / 399 / 440 | 399,42 | 59,94 | 116 bis 683 |
| 1 | 30.022 | 358 / 398 / 439 | 398,41 | 59,81 | 115 bis 682 |
| 2 | 30.050 | 359 / 399 / 439 | 398,54 | 59,51 | 119 bis 679 |

Mittel, Streuung und Grenzen stimmen mit der Ausgabe von bwa in allen angezeigten Stellen
überein. Erzeugt wurde mit Mittel 400 und Streuung 60; der Abstand der Anfänge ist die
Stücklänge minus 1.

## Ganze Zeilen gegen bwa

| | |
|---|---:|
| Paare | 100.000 |
| Paare, bei denen alle Zeilen beider Reads in allen Feldern gleich sind | 99.942 |
| Reads | 200.000 |
| Reads mit allen Feldern gleich | 199.908 |
| andere Zahl von Zeilen | 0 |
| Reads ohne Ort bei bwa | 972 |
| Reads mit ergänzenden Zeilen bei bwa | 0 |

Die 92 abweichenden Reads nach den Feldern, die sich unterscheiden:

| Felder | Reads |
|---|---:|
| Position und Position des Partners (dazu teils TLEN, AS, MD, NM) | 62 |
| XS (dazu teils XA) | 19 |
| MAPQ und MQ | 6 |
| Flag | 3 |
| XA | 2 |

## Einordnung
- **Position (62 Reads).** Paare in langen Wiederholungen am Zentromer mit MAPQ 0: Es gibt
  viele gleich gute Orte, und bwa wählt einen anderen als wir. Um die Stufe zu finden, wurden
  im ersten Block die Paar-Funktionen von bwa einzeln aufgerufen
  (`bench/m3/trace_bwa_pairs.py`) und Stufe für Stufe mit unseren Zwischenständen verglichen:
  Bei 51 von 66.668 Reads unterscheiden sich schon die Bereiche des einzelnen Reads vor jeder
  Paarung. Das Ergebnis der Paarbewertung weicht bei 8 von 33.334 Paaren ab, bei 7 davon waren
  die Bereiche vorher schon verschieden. An einem Paar nachvollzogen (`chr20:27503588:...:9977`):
  Read 1 hat 459 Bereiche, einer davon ist ein anderer als bei bwa; Partnersuche, Rangfolge und
  Paarbewertung liefern mit gleichen Bereichen dasselbe. Vermutete Ursache, nicht bewiesen:
  Seeds mit mehr als 500 Fundstellen werden ausgedünnt, und welche Fundstellen bleiben, hängt
  von der Ordnung im Index ab, die sich an N-Strecken von bwa unterscheidet (andere Füllbasen).
  Das ist dieselbe Gruppe wie die bisher ungeklärten Reads am Zentromer bei einzelnen Reads.
- **XS, XA (21 Reads).** bwa findet einen Seed in einer mit Zufallsbasen gefüllten N-Strecke
  oder verlängert ein Alignment um eine Base in sie hinein. An zwei Reads nachvollzogen
  (bwa meldet `chr20,-28646295,21S129M`, die erste dieser 129 Basen ist in der Referenz ein N),
  bei den übrigen nach dem gleichen Muster angenommen.
- **MAPQ (6 Reads, 3 Paare).** Paare am Zentromer; die MAPQ des Paares hängt an der Zahl fast
  gleich guter Paare, die sich aus dem ersten Punkt ergibt. Nicht einzeln geprüft.
- **Flag (3 Reads).** Read mit Ort auf dem Gegenstrang, Partner ohne Ort: Wir setzen Bit 32
  (Partner auf dem Gegenstrang), weil der Partner an den Ort des Reads gestellt wird; so steht
  es im Quelltext von bwa. pybwa liefert das Bit nicht, wie auch die Bits 2, 64 und 128. Klärt
  der Vergleich mit dem Programm bwa.

## Gegen den wahren Ort
`bench/m3/evaluate_pairs.py`. Gezählt sind die 199.028 Reads mit bekanntem Ort (alle außer den
972 fremden). Richtig heißt: linkes Ende höchstens 20 Basen neben dem wahren. Dieselben Reads
einmal als Paare und einmal einzeln ausgerichtet:

| | als Paare | einzeln |
|---|---:|---:|
| richtig | 194.778 (97,86 %) | 192.677 (96,81 %) |
| falsch | 4.250 | 6.351 |
| ohne Ort | 0 | 0 |
| mit MAPQ ab 30 | 191.846 | 188.328 |
| falsch mit MAPQ ab 30 | 0 | 1 |
| richtig, aber MAPQ 0 | 966 | 1.666 |
| fremde Reads mit Ort | 0 von 972 | 0 von 972 |

Die Paarung macht 2.101 Reads zusätzlich richtig und hebt 3.518 Reads auf MAPQ ab 30, ohne
einen falschen Ort mit hoher MAPQ. Die falschen Orte liegen fast alle bei MAPQ 0 in
Wiederholungen von chr20 (Zentromer); dort ist der Ort aus 150 Basen nicht bestimmbar.
Bei den 1.031 Paaren mit Read 2 aus einem anderen Stück bleibt Read 1 an seinem Ort
(1.996 von 2.062 Reads richtig).

## Nachtrag 19.11.2026: N-Stellen als unbekannt
Seit diesem Stand sieht der Aligner an N-Stellen der Referenz "unbekannt" (Abzug 1 wie bei N im
Read) statt der Füllbase des Index ([Read-Paare](../handbuch/verfahren/read-paare.md), Grenzen).
Die Zahlen oben sind neu gemessen. Gegenüber vorher ändern sich 3 von 200.000 Reads (zweimal MQ
und MAPQ um einen Punkt, einmal `NM` in `XA`); die Übereinstimmung mit bwa geht von 199.907 auf
199.908 Reads, die Zahl ganz gleicher Paare bleibt 99.942. Gegen den wahren Ort ändert sich
nichts.

## Was die Messung nicht zeigt
- Nur chr20 und erzeugte Reads ohne realistisches Fehlermodell. Echte Daten und das ganze
  Genom kommen mit GIAB.
- Kein Paar im Satz hat einen geteilten Read; dieser Zweig ist nur durch den Quelltext belegt.
- Flag-Bits 2, 64, 128 nicht gegen bwa geprüft (siehe oben), nur durch eigene Tests.
- Verglichen mit dem in pybwa eingebauten bwa 0.7.18, nicht mit dem Programm bwa 0.7.19.
- Tempo: 155 s für 100.000 Paare auf einem Kern, davon etwa die Hälfte Partnersuche
  (mehrere Kerne: [M3 Kerne und BAM](M3-kerne-bam.md)). Das war der Stand ohne Vektorbefehle
  im lokalen Alignment; neuer Stand: [M3 Tempo](M3-tempo.md).

Rohdaten: `bench/m3/ergebnisse/paare_vergleich.json`, `paare_stufen.json`,
`paare_wahrheit.json`, `paare_wahrheit_einzeln.json`.

## Nachtrag vom 19.11.2026
Der Vergleich nimmt jetzt zwei Abweichungen von pybwa gegenüber dem Programm bwa heraus (Flag 32
bei einem Partner ohne Ort, Basen hart gekürzter Zeilen auf dem Gegenstrang, siehe
[M5 Echte Reads](M5-echte-reads.md)). Damit gelten 3 Reads mehr als gleich: 199.911 Reads und
99.945 Paare. `bench/m3/ergebnisse/paare_vergleich.json` ist neu geschrieben. Dass doppelte
Seeds jetzt erhalten bleiben, ändert an den übrigen Abweichungen nichts (gleiche Zahlen je Feld).
