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
Zeilen schmal (rund 10 Zellen), der Gewinn ist deshalb begrenzt. Die Zeit ist noch nicht
gemessen.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Ohne Band und ohne Abbruch gegen eine unabhängige Rechnung mit drei vollen Matrizen (M, E, F getrennt geführt): Score, Endpunkt, Score am Read-Ende und dessen Zeile; 3.000 Zufallspaare mit Fehltreffern, Insertionen, Deletionen und N | gleich |
| Mit engem Band und Abbruch: Score nie höher als ohne, nie unter `h0` (1.000 Paare) | bestanden |
| Gesetzte Fälle: gleiche Folgen, Fehltreffer am Ende, leere Referenz | bestanden |
| Frühes Ende und Rechnung in Zeilen gegen die Rechnung Zelle für Zelle ohne frühes Ende: 20.000 Zufallspaare mit vier Bandbreiten, vier Abbruchgrenzen, zwei Punktzahlen je Treffer, gleichen Folgen, wenigen Fehlern, Lücken und N | alle sechs gemeldeten Größen gleich |
| Echte Reads (2.000 Paare, chr20): ganze SAM-Ausgabe vorher und nachher | gleich |

## Grenzen
- Kein direkter Vergleich mit `ksw_extend2` aus bwa; die Funktion ist dort nicht von außen
  aufrufbar. Ob Band und Abbruch in jedem Einzelfall wie bei bwa greifen, zeigt erst der Vergleich
  der fertigen Alignments.
- Die Punkte je Referenzbase stehen vorab in einer Tabelle. Tempo: [M3 Tempo](../../bench/M3-tempo.md).

## Literatur
Li 2013; Gotoh 1982 für affine Lücken, siehe [Literatur](../literatur.md).
