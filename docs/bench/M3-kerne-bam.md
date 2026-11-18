# M3: Mehrere Kerne und BAM-Ausgabe

Stand: 18.11.2026. Verfahren: [Ausrichten auf mehreren Kernen](../handbuch/verfahren/mehrere-kerne.md),
[BAM lesen und schreiben](../handbuch/verfahren/bam-schreiben.md).

## Fragen
1. Ist die Ausgabe für jede Zahl von Kernen dieselbe?
2. Enthält die BAM dieselben Zeilen wie die SAM-Ausgabe, und liest samtools sie ohne Beanstandung?
3. Wie viel schneller wird es mit mehr Kernen?

## Aufbau
- Referenz chr20 aus GRCh38, Index `chr20.kompakt.gcfm`.
- Paare: die 100.000 Paare aus [M3 Paare](M3-paare.md) (200.000 Reads, drei Blöcke).
- Einzeln: die 200.000 Reads aus [M3 SAM](M3-sam.md) (drei Blöcke).
- `bench/m3/measure_cores.py`: `genome index sam` mit `--kerne 1`, `2` und `4`; verglichen wird
  die ganze Ausgabe ohne die `@PG`-Zeile (dort steht der Aufruf). Danach ein Lauf mit `--bam`;
  `samtools quickcheck` und `samtools view -h` (samtools aus pysam) wandeln die BAM zurück in
  SAM-Text, der Zeile für Zeile mit der SAM-Ausgabe verglichen wird.
- Rechner: Cloud-Maschine mit 4 virtuellen Kernen (Intel Xeon, 2,3 GHz), geteilt mit anderen.
  Je Einstellung ein Lauf.

## Ergebnis

### Gleichheit

| Satz | Zeilen | 1, 2 und 4 Kerne | BAM gegen SAM |
|---|---|---|---|
| 100.000 Paare | 200.000 | Byte für Byte gleich | alle Zeilen gleich, `quickcheck` bestanden |
| 200.000 einzelne Reads | 200.000 | Byte für Byte gleich | alle Zeilen gleich, `quickcheck` bestanden |

Die Zeilen sind dieselben wie in den Vergleichen mit bwa ([Paare](M3-paare.md), [SAM](M3-sam.md));
an den dortigen Zahlen ändert sich nichts.

### Laufzeit

| Satz | 1 Kern | 2 Kerne | 4 Kerne | 4 Kerne, BAM | Reads je Sekunde (4 Kerne) |
|---|---|---|---|---|---|
| 100.000 Paare | 161,8 s | 76,2 s | 38,4 s | 38,5 s | 5.200 |
| 200.000 einzelne Reads | 72,1 s | 36,4 s | 15,3 s | 15,1 s | 13.100 |

Größe der BAM: 17,0 MB (Paare), 15,2 MB (einzeln).

Mit 4 Kernen ist es 4,2-mal (Paare) und 4,7-mal (einzeln) so schnell wie mit einem. Mehr als
4-mal kann die Verteilung allein nicht bringen; der Überschuss ist Streuung zwischen einzelnen
Läufen auf einer geteilten Maschine. Belegt ist damit: Die Arbeit verteilt sich bis 4 Kerne
ohne erkennbaren Verlust. Mehr sagt die Messung nicht.

## Was die Messung nicht zeigt
- Mehr als 4 Kerne. Der Zielrechner hat 10 (4 schnelle, 6 sparsame); das Verhalten dort ist offen.
- Das Tempo an sich ist noch weit vom Ziel: Für ein Genom mit 30-facher Abdeckung (rund 600
  Millionen Reads) wären es bei 5.200 Reads je Sekunde etwa 32 Stunden. Die Rechnung selbst
  (lokales Alignment, Verlängerung, Seed-Suche) ist noch nicht auf Tempo gebaut; das ist der
  nächste Schritt und wird an denselben Prüffällen gemessen.
- Kein Vergleich der Laufzeit mit bwa auf demselben Rechner (das Programm bwa steht hier nicht
  zur Verfügung).
- Je Einstellung nur ein Lauf, keine Streuung angegeben.
- Arbeitsspeicher nicht gemessen.

Rohdaten: `bench/m3/ergebnisse/kerne_paare.json`, `kerne_einzeln.json`.
