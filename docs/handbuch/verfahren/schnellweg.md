# Schnellweg für leichte Reads

Die meisten Reads eines menschlichen Genoms passen in voller Länge ohne Lücke an genau eine Stelle,
mit höchstens einer Abweichung. Für sie liefert das volle Verfahren (Seed-Suche mit Nachsuche,
Ketten, Verlängern) genau einen Bereich über den ganzen Read. Der Schnellweg findet diesen Bereich
direkt und überspringt die Seed-Suche. `fq2bam --schnellweg an` schaltet ihn ein; Vorgabe ist
`aus`. Quelltext: `crates/gc-align/src/fast.rs`.

## Ablauf
1. Den Read in Stücke von mindestens 21 Basen teilen (7 bei 148 Basen) und jedes Stück im
   FM-Index suchen, auf beiden Strängen.
2. Hat ein Stück mehr als 8 Fundstellen, geht der Read den vollen Weg.
3. Den Anfang wählen: Haben alle Stücke höchstens eine Fundstelle, gibt das erste Stück mit einer
   Fundstelle den Anfang (höchstens ein Stück darf keine haben). Sonst wird jede Fundstelle in einen
   Anfang des Reads umgerechnet; gewählt wird der Anfang, den alle Stücke bis auf höchstens eines
   teilen. Der ganze Read wird dort Base für Base mit der Referenz verglichen: höchstens eine
   Abweichung, kein N, innerhalb eines Contigs.
4. Jedes Stück ohne die Abweichung muss beim gewählten Anfang liegen, das Stück mit der Abweichung
   nicht. Sonst voller Weg.
5. Weitere Kopien: Alle anderen Anfänge sind weitere Stellen, an denen Teile des Reads genau
   passen. Liegt eine davon näher als Readlänge plus 256 Basen am gewählten Anfang (Wiederholung
   am Ort) oder gibt es mehr als 4 getrennte, geht der Read den vollen Weg. Sonst wird der ganze
   Read gegen jede Kopie lokal ausgerichtet, mit 128 Basen Rand auf jeder Seite. Ist der beste
   Score einer Kopie höher als die Grenze unten, geht der Read den vollen Weg.
6. Der Bereich bekommt die Punkte, die das volle Verfahren vergibt: über den ganzen Read, und als
   besten lokalen Wert den ohne das Ende hinter der Abweichung, wenn diese weniger als fünf Basen
   vor einem Ende liegt.

## Warum die MAPQ gleich bleibt
bwa senkt die MAPQ bei 148 Basen erst unter 60, wenn eine zweite Kopie höchstens etwa 16 Punkte
schlechter ist, also höchstens vier Ereignisse (Abweichung oder Lücke) vom Read entfernt. Eine
solche Kopie passt in mindestens drei der sieben Stücke genau. Höchstens eines davon ist das Stück
mit der Abweichung; ein anderes hat dann zwei Fundstellen, und der Read geht den vollen Weg.
Dasselbe gilt bis sechs Ereignisse, also für alle Kopien, die bwa als `XA` meldet (ab 0,8 des
besten Scores).

Eine erste Fassung mit drei Stücken war schneller, gab aber auf chr20 bei 16 Reads mit naher
zweiter Kopie MAPQ 60 statt 40 bis 57. Sie wurde verworfen.

## Grenze für weitere Kopien
Eine Kopie mit Score `x` neben dem besten Score `b` ändert die Ausgabe von bwa nur, wenn
- `x >= 0,8 b` (dann steht sie in `XA`), oder
- `b - x` unter etwa 17 liegt (dann sinkt die MAPQ unter 60, und die Kopie wird Ausgangspunkt der
  Partnersuche und kann ein Paar bilden, das mit dem besten Paar konkurriert).

Der Schnellweg verlangt `x < 0,8 b` und `x <= b - 31`, bei 148 Basen also `x <= 117`. Der Score,
den das volle Verfahren einer Kopie gibt, ist höchstens der beste lokale Score im Fenster, denn
auch das Verlängern vom Seed aus ergibt ein lokales Alignment in diesem Fenster. Die Kopie bleibt
im vollen Weg als Nebenbereich stehen und bestimmt dort `XS`; der Schnellweg schreibt `XS:i:0`.

`--schnellweg einfach` schaltet die Kopien ab: Dann schickt jedes Stück mit mehr als einer
Fundstelle den Read den vollen Weg (Stand vor dieser Erweiterung).

## Stücke gemeinsam suchen
Die Stücke eines Reads werden gemeinsam gesucht (`FmIndex::intervals_many`). Die letzten bis zu
zehn Basen jedes Stücks kommen fertig aus der Tabelle der kurzen Muster, die beim Laden des Index
mit denselben Schritten gerechnet wird. Die übrigen Schritte laufen für alle Stücke im
Gleichschritt, und nach jedem Schritt wird der Zählblock des nächsten schon angefordert, sodass
sich die Wartezeiten auf den Speicher überlappen. Die Bereiche sind dieselben wie bei der Suche
Stück für Stück (Test `intervals_many_equals_interval`). Der Schalter
`SeedOptions::fast_pieces_together = false` sucht jedes Stück wieder einzeln von vorn.

Seit 2026-11-27 läuft der Schnellweg in Gruppen von 32 Reads (`fast_regions_group`): Die Stücke
aller Reads der Gruppe werden zu je 32 verschränkt gesucht, und die Orte aller Reads gemeinsam
nachgeschlagen (`FmIndex::locate_many`), statt Read für Read einzeln durch die BWT zu gehen. Im
Profil auf echten Reads ging fast die Hälfte der Zeit des Schnellwegs in dieses Nachschlagen. Das
Ergebnis ist Read für Read dasselbe (Test `the_group_gives_the_same_regions_as_read_by_read`).
Schalter zurück: `GC_SCHNELL_GRUPPE=0`.

| Messung (2026-11-26) | einzeln | gemeinsam |
|---|---|---|
| Schnellweg allein, 200.000 HG002-MGI-Reads chr20, Index GRCh38, 1 Kern, A-B-B-A (`FAST=1 seedbench`) | 5,68 s / 4,98 s | 3,08 s / 2,97 s |
| `fq2bam` auf 200.000 dieser Paare | BAM byte-gleich | BAM byte-gleich |

## Kopien ohne Ausrichtung ausschließen
Vor der Ausrichtung einer weiteren Kopie wird eine Obergrenze ihres lokalen Scores bestimmt:
die Zahl der Read-Positionen, deren 5-mer im Referenzfenster vorkommt, plus 4 (eine Maske aus
1.024 Bits je Fenster). Ein lokales Alignment besteht aus exakten Läufen, getrennt durch
Fehltreffer (−4) oder Lücken (mindestens −7). Ein Lauf der Länge `L` bringt `L` Punkte und
enthält `max(L − 4, 0)` 5-mer-Starts, also höchstens 4 Punkte mehr; jede Trennung hebt diesen
Überschuss wieder auf. Liegt die Grenze bei oder unter der Grenze für Kopien, kann die Kopie den
Read nicht vom Schnellweg schicken, und ihre Ausrichtung entfällt.

Die Herleitung gilt nur ohne N: Gegen ein N kostet eine Base nur 1 Punkt. Enthält das Fenster
ein N (auch Lochstellen der Referenz), wird wie bisher ausgerichtet; ebenso bei einer anderen
Bewertung als der Vorgabe. Schalter: `SeedOptions::fast_copy_filter = false`. Test
`kmer_bound_is_never_below_the_local_score` vergleicht die Grenze mit dem SIMD-Score, auch mit
N, Fehlern und Lücken; ein früherer Entwurf ohne N-Prüfung scheiterte dort.

| Messung (2026-11-26) | Ergebnis |
|---|---|
| 36.805 Kopienfenster aus 200.000 HG002-MGI-Reads, Index GRCh38 | Grenze nie unter dem SIMD-Score; 95,4 % überspringbar |
| Schnellweg allein, 1 Kern, A-B-B-A (`FAST=1 FAST_FIELD=filter seedbench`) | 2,91 / 3,10 s ohne, 2,62 / 2,62 s mit |
| `fq2bam` auf 200.000 HG002-MGI-Paaren und 300.000 simulierten Paaren (ganzes GRCh38) | BAM byte-gleich |

## Was sich ändert
- `XS` bei Reads, für die das volle Verfahren nur zufällige kurze Treffer anderswo findet
  (Score unter 0,8 des besten, meist 19 bis 50). Der Schnellweg schreibt dort `XS:i:0`.
- In sehr seltenen Fällen bei Gleichstand ein anderer von zwei gleich guten Orten eines geretteten
  Partners (1 Read von 843.128 auf dem ganzen Genom, gleicher Score).

## Prüfmodus
`--schnellweg pruefen` rechnet für jeden Read beide Wege, gibt den vollen aus und zählt in
`ALIGN.FERTIG` unter `schnellweg`, wie oft der Schnellweg denselben Bereich geliefert hätte
(`gleich_wie_voller_weg`, `anders_als_voller_weg`) und warum Reads den vollen Weg gingen
(`mehrere_fundstellen`, `zu_wenig_stuecke`, `stuecke_uneinig`, `pruefung_nicht_bestanden`,
`kopie_zu_gut`). `kopien_geprueft` zählt die Reads, die den Schnellweg trotz weiterer Kopien nehmen.
Als gleich zählt ein Read, wenn der beste Bereich des vollen Wegs derselbe ist und alle übrigen
höchstens die Grenze für weitere Kopien erreichen.

## Messung
Siehe [M18 Schnellweg](../../bench/M18-schnellweg.md).

| Prüfung | Ergebnis |
|---|---|
| HG002 chr20, 843.128 echte Reads, Index chr20: Zeilen gegen `--schnellweg aus` | MAPQ, Flags, Orte, CIGAR gleich; `XS` anders in 9 % |
| dieselben Reads, Index GRCh38: Zeilen gegen `aus` | MAPQ gleich; `XS` anders in 44 %, 1 CIGAR anders (Gleichstand) |
| mit weiteren Kopien, 421.564 HG002-Paare, Index GRCh38: Anteil Schnellweg | 70,9 % statt 59,6 % (`kopien_geprueft` 94.927, `kopie_zu_gut` 2.046) |
| dieselben Paare: Zeilen gegen `aus` | MAPQ, Flags, Orte gleich; 1 CIGAR anders wie zuvor; `XS` anders in 54 % statt 44 % |
| dieselben Paare, Prüfmodus | 597.455 von 597.455 Schnellweg-Reads gleich wie der volle Weg |
| 400.000 simulierte Reads chr20 gegen die Wahrheit (`genome-validate`) | jeder Read gleich richtig oder falsch wie ohne Schnellweg |
| Tests: leichte Reads gleich wie der volle Weg, zweite Kopie und N gehen den vollen Weg | bestanden |

## Grenzen
- Gemessen bisher auf chr20-Reads einer Probe (HG002, Illumina 148 Basen). Das ganze Genom und
  andere Geräte stehen aus; bis dahin bleibt die Vorgabe `aus`.
- Reads unter 60 Basen gehen immer den vollen Weg.
- Eine Kopie, die in jedem der sieben Stücke eine Abweichung hat, bleibt unentdeckt. Fallen zwei
  dieser Abweichungen auf die letzte Base an beiden Enden, schneidet die lokale Ausrichtung sie ab,
  und die Kopie kann trotz sieben Abweichungen knapp über 0,8 des besten Scores liegen (`XA`).
  Das gilt für beide Fassungen; im Prüfmodus ist es bisher nicht aufgetreten.
