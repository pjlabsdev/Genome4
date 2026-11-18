# Suffix-Array nach SA-IS

Stand: gebaut und geprüft (`gc-index`, `sais.rs`).

## Aufgabe
Zu einem Text das Suffix-Array berechnen: die Anfangspositionen aller Suffixe in sortierter
Reihenfolge. Voraussetzung: das letzte Zeichen ist ein Endzeichen, das genau einmal vorkommt und
kleiner ist als alle anderen.

## Verfahren (Nong, Zhang und Chan 2009)

1. **Typen:** jedes Suffix ist vom Typ S (kleiner als das folgende) oder L (größer). Ein
   LMS-Suffix ist ein S-Suffix, dessen Vorgänger vom Typ L ist.
2. **LMS-Teilstücke sortieren:** LMS-Suffixe ans Ende der Fächer ihres ersten Zeichens legen, dann
   in zwei Durchgängen ableiten: von links alle L-Suffixe, von rechts alle S-Suffixe.
3. **Benennen:** gleiche LMS-Teilstücke bekommen denselben Namen. Die Namen in Textreihenfolge
   bilden einen verkürzten Text (höchstens halb so lang).
4. **Rekursion:** sind nicht alle Namen verschieden, wird das Suffix-Array des verkürzten Textes
   mit demselben Verfahren berechnet.
5. **Endgültig ableiten:** LMS-Suffixe in ihrer nun bekannten Reihenfolge einlegen und wie in
   Schritt 2 ableiten.

Laufzeit linear in der Textlänge. Speicher: 4 Bytes je Zeichen für das Suffix-Array, 1 Bit je
Zeichen für die Typen, dazu die Rekursion (höchstens die Hälfte).

## Grenze
Positionen sind 32 Bit breit: höchstens 4.294.967.294 Zeichen, also eine Referenz bis rund
2,1 Mrd. Basen (beide Stränge). Für das ganze Humangenom dient der
[blockweise Bau](blockweiser-bau.md).

## Prüfung

- **Tests:** Vergleich mit einfachem Sortieren aller Suffixe auf bekannten Wörtern und 400
  Zufallstexten mit erzwungenen Wiederholungen (Alphabetgrößen 1 bis 20); ein DNA-artiger Text mit
  30.000 Basen langer Wiederholung und Homopolymer.
- **`verify`:** prüft jedes Suffix-Array in linearer Zeit, unabhängig vom Bau: jede Position
  kommt genau einmal vor, und für benachbarte Zeilen gilt (erstes Zeichen, Rang des folgenden
  Suffixes) aufsteigend. Das ist die Prüfung nach Burkhardt und Kärkkäinen 2003. Der Test
  `verify_detects_errors` zeigt, dass vertauschte und doppelte Einträge auffallen.
- **Messung:** chr20, 128.888.335 Zeilen, vollständig geprüft ([M2](../../bench/M2-fm-index.md)).

## Literatur
Nong, Zhang & Chan 2009; Burkhardt & Kärkkäinen 2003. Siehe [Literatur](../literatur.md).
