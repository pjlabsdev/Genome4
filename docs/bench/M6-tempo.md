# M6: Tempo des Aligners, zweiter Schritt (Index)

Stand: 19.11.2026. Verfahren: [FM-Index](../handbuch/verfahren/fm-index.md),
[Seeds und Ketten](../handbuch/verfahren/seeds-und-ketten.md). Erster Schritt:
[M3 Tempo](M3-tempo.md).

## Fragen
1. Wo bleibt die Zeit bei echten Reads und beim ganzen Genom?
2. Was bringen die Änderungen am Index, und bleibt die Ausgabe gleich?
3. Was wurde versucht und verworfen?

## Aufbau
- Reads: GIAB HG002, Paare zu 148 Basen aus chr20:10 bis 12 Mb wie in
  [M5 Echte Reads](M5-echte-reads.md); dazu die 100.000 erzeugten Paare aus [M3](M3-paare.md).
- Index: GRCh38 ohne ALT (`grch38.kompakt.gcfm`, 4,84 GB) und chr20 in beiden Varianten.
- Vorher: Stand von main (8714a75). Nachher: dieser Stand. Beide mit denselben Einstellungen
  gebaut, im Wechsel gemessen.
- Messung: `bench/m6/measure_tempo.py`, je Programm drei Läufe; Zeit ist die des Ausrichtens aus
  dem Protokoll (`ALIGN.FERTIG`), ohne Laden des Index. Die Tabelle nennt den mittleren der drei
  Läufe. Rohdaten: `bench/m6/tempo.json`.
- Wo die Zeit bleibt: Der laufende Prozess wurde 200- bis 300-mal angehalten und der Ort
  notiert (ein Kern, ganzes Genom). Das misst Zeit samt Warten auf den Arbeitsspeicher, anders
  als das Zählen der Befehle in M3. Ungenau auf einige Prozentpunkte.
- Rechner: Cloud-Maschine mit 4 virtuellen Kernen (Intel Xeon, 2,1 GHz), 15 GB RAM, geteilt mit
  anderen. Der Rechner, auf dem die Läufe am Ende stattfinden (Apple M5), wurde nicht gemessen.

## Ergebnis

### Wo die Zeit blieb (echte Paare, ganzes Genom, ein Kern)

| Teil | vorher | nachher |
|---|---|---|
| Seed-Suche (SMEM, Nachsuche, dritte Runde) | 28 % | 39 % |
| Orten der Seeds | 28 % | 24 % |
| Lokales Alignment der Partnersuche | 18 % | 15 % |
| Verlängern der Seeds | 17 % | 17 % |
| Rest | 9 % | 5 % |

Mit dem ganzen Genom verschiebt sich das Bild gegenüber chr20: Der Index ist größer als jeder
Zwischenspeicher des Prozessors, die Arbeit im Index macht mehr als die Hälfte der Zeit aus.
Ein Read hat im Mittel 12 Treffer und 27 Seeds, die Seed-Suche braucht je Read rund 570
(chr20) bis 730 (ganzes Genom) Verlängerungen.

### Was geändert wurde
- **Schieben statt Dividieren.** Der Block und die Stelle im Block wurden bei jedem Zugriff mit
  einer Division durch den Abstand der Zählstände berechnet. Der Abstand steht in der Datei, der
  Übersetzer konnte daraus nichts Billigeres machen. Beide Abstände sind Zweierpotenzen; jetzt
  wird geschoben und maskiert. Das war der größte einzelne Gewinn beim Orten.
- **Orten im Gleichschritt.** Alle Seeds eines Reads werden gemeinsam geortet
  (`locate_many`): Die Wege von bis zu 32 Zeilen laufen nebeneinander, der Prozessor wartet auf
  mehrere Stellen des Arbeitsspeichers zugleich.
- **Beide Enden in einem Zugriff.** Liegen beide Enden eines Bereichs im selben Wort der BWT,
  wird nur einmal gezählt und das Stück dazwischen dazu (`occ4_pair`).
- **Zählen mit gesammelten Teilsummen** statt einer vollen Zählung je Wort.

In den Teilen für sich gemessen (50.000 Reads, ganzes Genom, ein Kern, mit einem
Hilfsprogramm, das nur für diese Messung bestand): Orten 2,04 s auf rund 1,55 s, Seed-Suche
3,0 s auf rund 2,7 s. Ausgeführte Befehle der Seed-Suche (callgrind, 3.000 Reads): 11 % weniger.

### Ausgabe
Unverändert. Die Prüfsumme der ganzen SAM-Ausgabe ohne `@PG` ist in jedem der sechs Sätze der
Tabelle vorher und nachher dieselbe. Der Stand vorher ist der, der in
[M5 Echte Reads](M5-echte-reads.md) Zeile für Zeile mit bwa mem verglichen wurde.

### Laufzeit

| Satz | Index | Kerne | vorher | nachher | Zeit nachher zu vorher |
|---|---|---|---|---|---|
| 50.000 echte Paare | GRCh38 kompakt | 1 | 17,1 s | 14,4 s | 0,84 |
| 50.000 echte Paare | GRCh38 kompakt | 4 | 4,88 s | 3,87 s | 0,79 |
| 100.000 echte Paare | chr20 kompakt | 1 | 16,5 s | 14,7 s | 0,89 |
| 100.000 echte Paare | chr20 schnell | 1 | 12,5 s | 12,0 s | 0,96 |
| 421.564 echte Paare | chr20 kompakt | 4 | 18,8 s | 17,1 s | 0,91 |
| 100.000 erzeugte Paare | chr20 kompakt | 1 | 84,6 s | 77,7 s | 0,92 |

Am ganzen Genom, dem Fall, auf den es ankommt, braucht das Ausrichten rund ein Fünftel weniger
Zeit. Auf chr20 ist der Gewinn kleiner, weil der kleine Index zum guten Teil im
Zwischenspeicher liegt. Gegen bwa mem (pybwa, 16.600 Reads je Sekunde und Kern auf chr20,
[M5 Echte Reads](M5-echte-reads.md)): mit dem kompakten Index jetzt 13.600 Reads je Sekunde, also
rund 1,2-mal die Zeit von bwa statt 1,3-mal. Am ganzen Genom gibt es noch keine Messung von bwa.

### Nachtrag: Verlängern der Seeds
Zwei Änderungen an [Seed verlängern](../handbuch/verfahren/seed-verlaengern.md), beide mit
gleichem Ergebnis in 20.000 Zufallsfällen und gleicher Ausgabe in den Sätzen unten:
- **Frühes Ende**, sobald keine spätere Zeile das Ergebnis noch ändern kann.
- **Rechnung in Zeilen** mit 16-Bit-Zahlen in Blöcken von acht Zellen.

Ausgeführte Befehle in der Verlängerung (callgrind, 2.000 echte Paare, chr20): 566 Millionen
vorher, 478 Millionen mit frühem Ende, 408 Millionen mit beidem (28 % weniger, im ganzen Lauf
5 %). Bei echten Reads hat ein Aufruf im Mittel rund 160 Zeilen zu rund 10 Zellen; für so
schmale Zeilen bringen Vektorbefehle wenig.

Laufzeit, je ein bis drei Läufe im Wechsel, Rechner zur Zeit der Messung unruhig:

| Satz | Index | Kerne | Stand der Tabelle oben | mit Verlängerung neu |
|---|---|---|---|---|
| 50.000 echte Paare | GRCh38 kompakt | 1 | 14,5 s und 15,4 s | 14,6 s und 14,7 s |
| 100.000 echte Paare | chr20 kompakt | 1 | 15,3 s, 15,4 s, 16,9 s | 15,1 s, 15,5 s, 16,2 s |
| 100.000 erzeugte Paare | chr20 kompakt | 4 | 20,0 s | 19,0 s |

Mit echten Reads ist kein Unterschied messbar, mit erzeugten rund 5 %.

### Nachtrag: dichtere Tabellen im Index
Der Index wurde versuchsweise mit anderen Abständen gebaut (Zählstand alle `b` Basen,
Suffix-Array-Stichprobe jede `s`-te Zeile). Die Ausgabe ist in jedem Fall dieselbe.

chr20, 100.000 echte Paare, ein Kern, mittlerer von drei Läufen:

| `b` | `s` | Zeit | zu kompakt | Größe am ganzen Genom (gerechnet) |
|---|---|---|---|---|
| 128 | 32 (kompakt) | 15,1 s | 1,00 | 4,8 GB |
| 128 | 16 | 13,9 s | 0,92 | 5,8 GB |
| 64 | 32 | 13,8 s | 0,92 | 6,4 GB |
| 64 | 16 | 13,2 s | 0,87 | 7,4 GB |
| 128 | 8 | 12,9 s | 0,86 | 7,7 GB |
| 64 | 8 (schnell) | 12,5 s | 0,83 | 9,3 GB |

Ganzes Genom mit `b` = 128 und `s` = 8 (Datei 7,75 GB, Bau 18,3 min auf 4 Kernen, Stichprobe
von 100.000 Abfragen ohne Abweichung), 50.000 echte Paare, ein Kern, je zwei Läufe: 12,7 s und
13,0 s gegen 14,6 s und 14,7 s mit dem kompakten Index, also rund 0,87-mal die Zeit für 2,9 GB
mehr Arbeitsspeicher.

Entscheidung vom 19.11.2026: `b` = 128 und `s` = 8 ist als Variante `standard` die Vorgabe von
`genome index build`. `kompakt` und `schnell` bleiben wählbar, vorhandene Indexdateien bleiben
lesbar. Die Messungen weiter oben in diesem Bericht gelten für den Index, der dort genannt ist.
Der Spitzenbedarf von `genome fq2bam` am ganzen Genom (5,24 GB mit kompakt, siehe
[M5](M5-ganzes-genom.md)) steigt mit der Vorgabe rechnerisch um 2,9 GB auf rund 8,2 GB.

Gemessen nach dem Einbau (GRCh38 ohne ALT, 421.564 echte Paare, 4 Kerne, je ein Lauf von
`genome fq2bam`, Spitze als `VmHWM` des Prozesses): 8,07 GB mit `standard` gegen 5,24 GB mit
`kompakt`, also 2,83 GB mehr. Ausrichten 26,7 s gegen 36,7 s. Die beiden BAM-Dateien sind in
allen 844.609 Datensätzen und im Kopf bis auf die `@PG`-Zeile gleich. Der Spitzenbedarf beim Bau
des Index mit `standard` ist nicht gemessen.

### Nachtrag: Tabelle für kurze Muster in der Seed-Suche
Profil nach dem Einbau des Index `standard` (ganzes Genom, echte Reads, 161 Stichproben des
Aufrufstapels beim Ausrichten): Seed-Suche 50 %, Ortung 15 %, Verlängern 12 %. Gezählt an 50.000
Paaren laufen 42 % aller Schritte der Seed-Suche auf Bereichen mit mehr als 1.000 Zeilen, 30 %
auf 9 bis 1.000, 11 % auf 2 bis 8 und 17 % auf genau einer Zeile.

Die Bereiche aller Muster bis Länge 10 werden deshalb beim Laden einmal berechnet und in der
Suche nachgeschlagen (siehe [Seed-Suche](../handbuch/verfahren/seed-suche.md)). GRCh38 mit dem
Index `standard`, 50.000 echte Paare, Sekunden des Ausrichtens, zwei Messreihen zu je drei
Läufen:

| Kerne | vorher | nachher | Verhältnis der mittleren Läufe |
|---|---|---|---|
| 1 | 12,4 s, 12,5 s, 12,7 s | 11,1 s, 11,1 s, 11,4 s | 0,89 |
| 1 (Wiederholung) | 12,4 s, 12,9 s, 13,1 s | 11,3 s, 11,7 s, 11,7 s | 0,91 |
| 4 | 3,57 s, 3,65 s, 3,78 s | 3,04 s, 3,07 s, 3,08 s | 0,84 |

Die erste Messreihe lief zum Teil neben einem Übersetzungslauf, die Wiederholung ungestört. Die
Ausgabe ist in allen Läufen dieselbe (md5 ohne `@PG`). Die Tabelle braucht rund 34 MB, ihr Bau
beim Laden ist in der Gesamtlaufzeit nicht messbar (10,4 s gegen 10,4 s für 2.000 Paare samt
Laden des Index).

Eine Tabelle bis Länge 12 (rund 540 MB) war nicht schneller (11,2 s, 11,7 s, 12,0 s auf einem
Kern; 3,04 s, 3,08 s, 3,10 s auf vier) und verlängert das Laden um 2,5 s. Es bleibt bei 10.

### Nachtrag: Zählblöcke vorab anfordern
Profil über einen ganzen Lauf (GRCh38, Index `standard`, 421.564 echte Paare, ein Kern, 594
Stichproben des Aufrufstapels): Seed-Suche 37 %, lokales Ausrichten 18 %, Verlängern 16 %,
Ortung 8 %, Laden des Index 7 %, Ketten 4 %.

Beim Schritt nach links in der Seed-Suche werden jetzt die Zählblöcke aller lebenden Bereiche
erst angefordert und dann gezählt (siehe [Seed-Suche](../handbuch/verfahren/seed-suche.md)).
Sekunden des Ausrichtens, vorher und nachher abwechselnd gestartet:

| Satz | Kerne | vorher | nachher |
|---|---|---|---|
| 50.000 echte Paare | 1 | 11,4 s, 11,0 s, 11,4 s | 11,1 s, 11,2 s, 11,6 s |
| 50.000 echte Paare | 4 | 3,55 s, 3,20 s, 3,37 s | 3,15 s, 3,33 s, 3,36 s |
| 421.564 echte Paare | 4 | 32,9 s, 26,2 s, 25,4 s | 26,4 s, 25,6 s, 24,9 s |

Der erste Lauf vorher (32,9 s) lief mit kaltem Plattencache und zählt nicht. In jedem Paar
direkt hintereinander gestarteter Läufe auf dem großen Satz ist die neue Fassung rund 2 %
schneller; das liegt an der Grenze dessen, was auf dieser Maschine messbar ist. Ein erster
Versuch mit Umschalter im selben Programm zeigte rund 5 %. Ausgabe in allen Läufen gleich.

Die Seed-Suche bleibt damit der größte Teil. Was übrig ist, hängt Schritt für Schritt vom
vorigen ab (Schritt nach rechts) oder ist schon nachgeschlagen; ein weiterer großer Gewinn ist
dort nicht in Sicht.

### Versucht und verworfen
- **Zählen ohne Verzweigung** (immer alle Wörter des Blocks, nicht gefragte Basen ausgeblendet):
  24 % mehr Befehle in der Seed-Suche und langsamer, wieder ausgebaut.
- **Zählblöcke an Speicherzeilen ausrichten** (64 Bytes, ein Block der Variante kompakt in genau
  einer Zeile): kein messbarer Unterschied, wieder ausgebaut.
- **Seed-Suche mehrerer Reads abwechselnd führen** (je Schritt erst für 16 Reads die Blöcke
  anfassen, dann rechnen), geprüft gegen die einfache Suche und mit gleicher Ausgabe: 10 bis
  25 % langsamer als die einfache Suche, wieder ausgebaut. Die Suche hängt auf diesem Rechner am
  Rechnen (rund 470 Befehle je Verlängerung), und die beiden Zugriffe einer Verlängerung laufen
  schon von selbst nebeneinander. Die Verwaltung der unterbrechbaren Suche kostete mehr, als das
  Warten einbrachte.
- **Orten über mehrere Reads hinweg bündeln** statt je Read: höchstens 7 % des Ortens, im
  ganzen Lauf nicht messbar, weil ein Read mit im Mittel 27 Seeds die Wege schon füllt. Nicht
  übernommen.
- **Bau nur für den Prozessor des Messrechners** (unter anderem mit dessen Befehl zum Zählen von
  Bits): rund 4 % auf chr20. Nicht übernommen, das Programm soll auf jedem Rechner dasselbe sein.

## Grenzen
- Gemessen auf einer geteilten Cloud-Maschine mit kleinen Speicherseiten. Auf dem Mac kann das
  Verhältnis von Rechnen zu Warten auf den Speicher anders sein; dort steht eine Messung aus.
- Große Speicherseiten des Betriebssystems für den Index wurden nicht untersucht. Ein kurzer
  Versuch mit der Systemeinstellung dafür war nicht eindeutig.
- 50.000 Paare am ganzen Genom sind ein kleiner Satz aus einem Ausschnitt von chr20. Ein Lauf
  mit voller Abdeckung steht aus.
- Die Seed-Suche ist jetzt der größte Teil (40 % der ausgeführten Befehle bei echten Reads auf
  chr20). Ihr Schritt kostet rund 470 Befehle, davon rund 200 je Zählung im Index; billiger
  wird er nur mit einem anderen Aufbau der Zählblöcke.
