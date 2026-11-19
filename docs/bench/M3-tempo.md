# M3: Tempo des Aligners, erster Schritt

Stand: 19.11.2026. Verfahren: [Read-Paare](../handbuch/verfahren/read-paare.md) (lokales
Alignment), [Seed verlängern](../handbuch/verfahren/seed-verlaengern.md).

## Fragen
1. Wo bleibt die Rechenzeit?
2. Was bringt das lokale Alignment mit Vektorbefehlen, und bleibt die Ausgabe gleich?
3. Wie weit sind wir je Kern von bwa entfernt?

## Aufbau
- Daten und Rechner wie in [M3 Kerne und BAM](M3-kerne-bam.md): chr20, 100.000 Paare und
  200.000 einzelne Reads, Cloud-Maschine mit 4 virtuellen Kernen, geteilt mit anderen, je
  Einstellung ein Lauf.
- Wo die Zeit bleibt: `valgrind --tool=callgrind` auf den ersten 2.000 Paaren (vorher) und 500
  Paaren (nachher), ein Kern. Gezählt werden ausgeführte Befehle, nicht Zeit; Wartezeit auf den
  Arbeitsspeicher (Index-Zugriffe) ist darin nicht enthalten.
- Unsere Laufzeit: `bench/m3/measure_cores.py` (prüft zugleich, dass die Ausgabe für jede
  Kernzahl gleich ist und die BAM stimmt).
- bwa: `bench/m3/time_bwa.py`, bwa mem aus pybwa 2.3.2 (bwa 0.7.18), ein Kern, dieselben Reads.
  Gemessen ist nur der Aufruf des Aligners. pybwa legt je Zeile ein Objekt an, das Programm bwa
  schreibt Text; dafür ist bei uns das Lesen der FASTQ und das Schreiben der Ausgabe in der Zeit
  enthalten. Der Vergleich gibt die Größenordnung, keine genaue Zahl.

## Ergebnis

### Wo die Zeit blieb (Paare, Anteil der ausgeführten Befehle)

| Teil | vorher | nachher |
|---|---|---|
| Lokales Alignment der Partnersuche | 75 % | 9 % |
| Verlängerung der Seeds | 15 % | 51 % |
| Zählen im Index (`occ`) | 5 % | 17 % |
| Rest | 5 % | 23 % |

Je Paar wurden vorher 46,8 Millionen Befehle ausgeführt, nachher 14,4 Millionen.

### Was geändert wurde
- **Lokales Alignment.** Die Spuren sind jetzt feste Felder (16 Bytes oder 8 Werte mit 16 Bit),
  die Schritte gesättigte Additionen, Subtraktionen und Maxima. Der Übersetzer macht daraus
  Vektorbefehle. Kein systemnaher Code, deshalb auf macOS, Linux und Windows dasselbe Programm.
- **Verlängerung der Seeds.** Punkte je Referenzbase vorab als Tabelle, innere Schleife ohne
  Grenzprüfungen. Kleiner Gewinn.

### Ausgabe
Unverändert. Die Prüfsumme der ganzen Ausgabe (ohne `@PG`) ist für Paare und einzelne Reads
dieselbe wie vor der Änderung, mit 1, 2 und 4 Kernen. Die schnelle Fassung des lokalen
Alignments liefert in 3.000 Zufallsfällen jeden Wert wie die erste Fassung und in den 1.000
Prüffällen jeden Wert wie `ksw_align2` aus bwa.

### Laufzeit

| Satz | Kerne | vorher | nachher | bwa (pybwa) |
|---|---|---|---|---|
| 100.000 Paare | 1 | 161,8 s | 89,1 s | 67,7 s |
| 100.000 Paare | 2 | 76,2 s | 44,7 s | |
| 100.000 Paare | 4 | 38,4 s | 22,1 s | |
| 200.000 einzelne Reads | 1 | 72,1 s | 68,7 s | 59,5 s |
| 200.000 einzelne Reads | 2 | 36,4 s | 34,4 s | |
| 200.000 einzelne Reads | 4 | 15,3 s | 17,0 s | |

- Paare sind 1,8-mal so schnell wie vorher. Bei einzelnen Reads liegt die Änderung im Bereich
  der Streuung zwischen Läufen (die Zeile mit 4 Kernen ist sogar langsamer gemessen).
- Je Kern brauchen wir rund 1,3-mal (Paare) und 1,15-mal (einzeln) so lange wie bwa.
- Die Paarung kostet bei uns rund 18 s auf einem Kern: Dieselben 200.000 Reads einzeln
  ausgerichtet brauchen 71 s.

## Zweiter Schritt: Zählen im Index für alle vier Basen

Die Seed-Suche braucht bei jeder Verlängerung die Zählung aller vier Basen an beiden Enden des
Bereichs. Bisher waren das acht einzelne Aufrufe von `occ`, jetzt zwei von `occ4`, die den
Zählblock je einmal durchgehen ([FM-Index](../handbuch/verfahren/fm-index.md)). Die Ausgabe ist
wieder unverändert (gleiche Prüfsummen, beide Index-Varianten, 1 und 4 Kerne).

Seed-Suche allein (`genome index smem`, 200.000 Reads, ein Kern): 8,5 s statt 12,1 s.

| Satz | Kerne | Index kompakt | Index schnell | bwa (pybwa) |
|---|---|---|---|---|
| 100.000 Paare | 1 | 81,8 s | 68,0 s | 67,7 s |
| 100.000 Paare | 4 | 21,1 s | 17,2 s | |
| 200.000 einzelne Reads | 1 | 61,8 s | 49,0 s | 59,5 s |
| 200.000 einzelne Reads | 4 | 15,4 s | 13,1 s | |

- Mit dem kompakten Index (12,5 Bit je Base) liegen wir je Kern bei rund 1,2-mal (Paare) und
  1,04-mal (einzeln) der Zeit von bwa.
- Mit dem schnellen Index (24 Bit je Base, Stichprobe des Suffix-Arrays alle 8 statt 32 Zeilen)
  sind es Gleichstand bei Paaren und 0,82-mal bei einzelnen Reads. Der Index von bwa hat
  Zählblöcke alle 128 und Stichproben alle 32 Zeilen, entspricht also dem kompakten; der
  schnelle tauscht Speicher gegen Zeit. Für das ganze Genom wäre er rund 9 GB groß statt 4,8 GB.
- Rohdaten: `bench/m3/ergebnisse/tempo2_paare_kompakt.json`, `tempo2_paare_schnell.json`,
  `tempo2_einzeln_kompakt.json`, `tempo2_einzeln_schnell.json`.

## Einordnung
Der Aligner ist ein Nachbau von bwa mem und rechnet dieselben Schritte. Mehr als etwa
Gleichstand mit bwa je Kern ist auf diesem Weg nicht zu erwarten. Was darüber hinausgeht,
braucht andere Mittel: Index-Zugriffe mit Vorabladen, Verlängerung vieler Seeds gleichzeitig mit Vektorbefehlen (der Weg von BWA-MEM2,
Vasimuddin et al. 2019). Das sind die nächsten Schritte; die Ausgabe muss dabei Byte für Byte
gleich bleiben.

## Was die Messung nicht zeigt
- Nur chr20. Beim ganzen Genom ist der Index 50-mal so groß und passt nicht mehr in den
  schnellen Zwischenspeicher; der Anteil der Index-Zugriffe wird dort deutlich höher sein.
- Der Vergleich mit bwa ist ungenau (siehe Aufbau) und läuft gegen 0.7.18 aus pybwa, nicht
  gegen das Programm bwa 0.7.19.
- Je Einstellung ein Lauf auf einer geteilten Maschine. Unterschiede unter etwa 10 % sind nicht
  belegt.
- Nur Intel mit SSE2. Auf ARM (NEON) laufen die Tests im CI, das Tempo ist dort nicht gemessen.
- Arbeitsspeicher nicht gemessen.

Rohdaten: `bench/m3/ergebnisse/tempo_paare.json`, `tempo_einzeln.json`, `tempo_bwa_paare.json`,
`tempo_bwa_einzeln.json`; Stand vorher: `kerne_paare.json`, `kerne_einzeln.json`.
