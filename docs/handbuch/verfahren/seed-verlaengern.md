# Seed verlängern: Smith-Waterman mit Band und affinen Lücken

Stand: Rechenkern gebaut und geprüft (`gc-align`, `extend.rs`), noch ohne SIMD. Noch nicht
gebaut: das Zusammensetzen zum Alignment eines Reads (links und rechts verlängern, CIGAR, Ende des
Reads abschneiden oder nicht), Paarung und MAPQ.

## Aufgabe
Eine Kette sagt, wo der Read ungefähr liegt. Das Alignment entsteht, indem man vom Seed aus nach
links und nach rechts verlängert und dabei Fehltreffer und Lücken zulässt. Der Rechenkern
entspricht `ksw_extend2` aus bwa (Li 2013).

## Rechnung
Verankert am Seed mit dessen Score `h0`. Für Zeile i (Referenz) und Spalte j (Read):

```text
M(i,j) = H(i-1,j-1) + s(i,j)      falls H(i-1,j-1) > 0, sonst tot
E(i,j) = max(M(i-1,j) - o_del - e_del, E(i-1,j) - e_del)      Deletion
F(i,j) = max(M(i,j-1) - o_ins - e_ins, F(i,j-1) - e_ins)      Insertion
H(i,j) = max(M, E, F, 0)
```

- Ein Pfad endet, sobald sein Score auf 0 fällt; er beginnt nirgends neu (anders als beim
  lokalen Smith-Waterman).
- Lücken öffnen nur aus M. Eine Insertion direkt nach einer Deletion entsteht so nicht.
- Bewertung wie BWA-MEM: gleiche Base +1, ungleiche -4, N gegen irgendetwas -1, Lücke öffnen 6,
  verlängern 1 (für Insertion und Deletion getrennt einstellbar).

Gemeldet werden der höchste Score mit den dafür verbrauchten Basen von Read und Referenz, und
getrennt der höchste Score, bei dem der Read bis zu seinem Ende ausgerichtet ist. Aus dem
Vergleich beider entscheidet später das Alignment, ob das Read-Ende abgeschnitten wird.

## Band und Abbruch
- **Band:** Gerechnet wird nur bis `w` Spalten (Vorgabe 100) neben der Diagonale. Das Band wird
  zusätzlich auf die längste Lücke begrenzt, die sich bei dieser Read-Länge noch lohnen kann,
  und nach jeder Zeile auf die noch lebenden Spalten eingeengt.
- **Abbruch (Z-Drop):** Fällt der beste Score einer Zeile um mehr als 100 unter das bisherige
  Maximum, nachdem der Versatz zur Diagonale als Lücke angerechnet ist, endet die Rechnung.

Beides kann nur Ergebnisse wegnehmen, nie bessere erfinden.

- **Frühes Ende:** Nach jeder Zeile, in der das Read-Ende schon einmal erreicht war, wird
  geprüft, ob noch etwas zu gewinnen ist. Jeder Pfad in spätere Zeilen kommt aus einer lebenden
  Zelle der fertigen Zeile und kann von dort höchstens noch für jede übrige Base des Reads einen
  Treffer bekommen. Liegt diese Obergrenze für alle Zellen nicht über dem bisherigen Maximum und
  unter dem besten Score am Read-Ende, ändert sich keine gemeldete Größe mehr, und die Rechnung
  endet. bwa rechnet in dem Fall weiter bis zum Z-Drop und kommt zum selben Ergebnis.

## Rechnung in Zeilen
Passen die Zahlen in 16 Bit (bei Reads bis einige tausend Basen immer), wird jede Zeile in
Blöcken von acht Zellen gerechnet (`extend_rows`), sonst Zelle für Zelle (`extend_cells`). In
einer Zeile läuft nur die Insertion `F` von Zelle zu Zelle weiter: Sie ist das Maximum über alle
Zellen links davon, je Schritt um die Verlängerungskosten gesunken. Mit einem Ausgleich von
`k` mal Verlängerungskosten wird daraus ein laufendes Maximum; alles andere hängt nur von der
vorigen Zeile ab und lässt sich für acht Zellen zugleich rechnen. Bei echten Reads sind die
Zeilen schmal (rund 10 Zellen bei rund 160 Zeilen je Aufruf), der Gewinn ist deshalb begrenzt:
28 % weniger ausgeführte Befehle in der Verlängerung, in der Laufzeit mit echten Reads kein
messbarer Unterschied, mit erzeugten Reads rund 5 % ([M6 Tempo](../../bench/M6-tempo.md)).

Seit M18 rechnen die Blöcke mit Vektoren aus `gc-sys` (`simd.rs`, SSE2 oder NEON) statt mit
Feldern, die der Übersetzer nur zum Teil in Vektorbefehle verwandelte. Auch das laufende Maximum
läuft jetzt im Vektor: je Block drei Schritte (um 1, 2 und 4 Zellen verschoben) und der Übertrag
aus dem Block davor. Zusammen mit dem lokalen Alignment der Partnersuche sinkt die Zeit fürs
Ausrichten auf dem ganzen GRCh38 von 31,1 auf 25,8 s (echte Reads, 4 Kerne), bei gleicher BAM
([M18](../../bench/M18-schnellweg.md)). `GENOME_SIMD=0` rechnet wieder mit Feldern.

## Gleiche Verlängerungen nur einmal rechnen
Ein Read aus einer Wiederholung (Alu, L1, Satelliten) passt an viele fast gleiche Kopien und hat
dann viele Ketten. Oft verlängern mehrere davon exakt dasselbe: gleiches Read-Stück, gleiche
Referenzfolge, gleicher Startwert `h0`. Bei PJ (1 Mio. Paare, MGI 2x150, Schnellweg an,
Partnersuche für 20 Bereiche) sind 36,7 % der Verlängerungen innerhalb eines Reads doppelt, mit
26,6 % der Zellen. Bei 7,8 % der Paare mit vielen Bereichen kommen auf einen Read rund 98 Ketten,
sonst rund 5.

`extend_gemerkt` (`region.rs`) merkt sich je Kern die Verlängerungen des aktuellen Reads und
rechnet eine Verlängerung mit genau demselben Inhalt nicht noch einmal. Der Schlüssel ist der
vollständige Inhalt, kein Hash allein: Bewertung (alle Punkte, Bandbreite, Z-Drop), `end_bonus`,
`h0`, Read-Stück und Referenzstück. Vom Referenzstück zählen nur die ersten
`Read-Länge + Bandbreite + 1` Basen. In Zeile `i` beginnt das Band bei Spalte `i - w`; ab
`i = Read-Länge + w` ist es leer, der Zeilenhöchstwert ist 0 und die Rechnung endet dort. Der Rest
der Referenz kann das Ergebnis nicht beeinflussen. Ein Treffer liefert deshalb dasselbe Ergebnis
wie die Rechnung. Das Gedächtnis wird vor jedem Read geleert (`neuer_read`) und hält höchstens
8.192 Einträge je Kern.

`fq2bam --verlaengerung-merken aus` rechnet jede Verlängerung neu, wie bwa.

| Messung (Mac M5, 1 Mio. PJ-Paare, A-B-B-A mit 3 min Pause) | ohne | mit |
|---|---|---|
| Wandzeit, Mittel | 34,82 s | 32,93 s (-5,4 %) |
| Rechenzeit aller Kerne, Mittel | 288,7 s | 273,4 s (-5,3 %) |
| BAM (`samtools view`, md5) | 1a02646f… | gleich |

Inhaltsgleiche Fenster der Partnersuche gibt es dagegen kaum (1,5 %); dort lohnt das Merken nicht.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Ohne Band und ohne Abbruch gegen eine unabhängige Rechnung mit drei vollen Matrizen (M, E, F getrennt geführt): Score, Endpunkt, Score am Read-Ende und dessen Zeile; 3.000 Zufallspaare mit Fehltreffern, Insertionen, Deletionen und N | gleich |
| Mit engem Band und Abbruch: Score nie höher als ohne, nie unter `h0` (1.000 Paare) | bestanden |
| Gesetzte Fälle: gleiche Folgen, Fehltreffer am Ende, leere Referenz | bestanden |
| Frühes Ende und Rechnung in Zeilen gegen die Rechnung Zelle für Zelle ohne frühes Ende: 20.000 Zufallspaare mit vier Bandbreiten, vier Abbruchgrenzen, zwei Punktzahlen je Treffer, gleichen Folgen, wenigen Fehlern, Lücken und N | alle sechs gemeldeten Größen gleich |
| Echte Reads (2.000 Paare, chr20): ganze SAM-Ausgabe vorher und nachher | gleich |
| Merken an und aus: alle Bereiche von Reads aus einer Referenz mit 12 Kopien eines Elements (Fehltreffer, Insertion, Deletion), Gedächtnis nachweislich benutzt | gleich |
| Gemerkte Verlängerung mit gekürztem Referenzstück gegen volle Rechnung: 300 Zufallspaare, vier Bandbreiten, Referenz bis 400 Basen länger als erreichbar | gleich |
| 1 Mio. PJ-Paare, ganze BAM mit und ohne Merken | md5 gleich |

## Grenzen
- Kein direkter Vergleich mit `ksw_extend2` aus bwa; die Funktion ist dort nicht von außen
  aufrufbar. Ob Band und Abbruch in jedem Einzelfall wie bei bwa greifen, zeigt erst der Vergleich
  der fertigen Alignments.
- Die Punkte je Referenzbase stehen vorab in einer Tabelle. Tempo: [M3 Tempo](../../bench/M3-tempo.md).

## Literatur
Li 2013; Gotoh 1982 für affine Lücken, siehe [Literatur](../literatur.md).
