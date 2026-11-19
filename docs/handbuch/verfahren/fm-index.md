# FM-Index

Stand: gebaut und geprüft (`gc-index`, `fm.rs`).

## Text
Wie bei BWA deckt der Index beide Stränge ab. Der Text ist der Vorwärtsstrang aller Contigs
hintereinander, gefolgt von seinem reversen Komplement, gefolgt vom Endzeichen. Bei einer Referenz
mit L Basen hat die BWT 2L + 1 Zeilen.

**N-Strecken:** N und Mehrdeutigkeits-Codes haben keinen 2-Bit-Code. Sie werden mit festen
Pseudozufallsbasen gefüllt (xorshift64, fester Startwert, also auf jedem System gleich) und als
Strecken vermerkt. Das Prinzip stammt von BWA; die gefüllten Basen selbst sind andere.

## Bestandteile

| Teil | Inhalt |
|---|---|
| BWT | je Zeile die Base vor dem Suffix, 2 Bit, 32 Basen je 64-Bit-Wort |
| Zählblöcke | vor je `b` Basen der BWT der Zählstand von A, C, G, T bis dorthin (4 × 64 Bit) |
| Primärzeile | Zeile, deren Suffix an Position 0 beginnt; dort steht in der BWT das Endzeichen (gespeichert als A, beim Zählen abgezogen) |
| `c_before` | Zahl der Zeilen vor dem ersten Suffix, das mit A, C, G, T beginnt |
| Suffix-Array-Stichproben | die Textposition jeder `s`-ten Zeile, 40 Bit |
| Referenz | Vorwärtsstrang 2-Bit-gepackt, Contig-Tabelle, N-Strecken |

## Varianten

| Variante | `b` (Zählblock) | `s` (Stichprobe) | Bits je Referenzbase |
|---|---|---|---|
| kompakt | 128 | 32 | 12,5 |
| schnell | 64 | 8 | 24 |

Die Zahl gilt für beide Stränge samt gepackter Referenz. Rechnung für kompakt, je Zeile der BWT:
Zählblock (256 + 256 Bit auf 128 Zeilen) 4 Bit, Stichproben 40 / 32 = 1,25 Bit; mal zwei Zeilen
je Referenzbase, plus 2 Bit Referenz.

## Suche
`occ(c, i)` zählt die Base c in den Zeilen 0 bis i-1 der BWT: Zählstand des Blocks plus Auszählen
der Wörter bis i (Bitmaske und Popcount). `occ4(i)` liefert die Zählung für alle vier Basen in
einem Durchgang durch denselben Block (wie `bwt_occ4` in bwa): je Wort drei Popcounts, A ergibt
sich als Rest. Die Seed-Suche braucht immer alle vier. Rückwärtssuche über das Muster von hinten nach vorn:

```text
lo = 0; hi = Zeilen
für c im Muster, von der letzten Base zur ersten:
    lo = c_before[c] + occ(c, lo)
    hi = c_before[c] + occ(c, hi)
```

`[lo, hi)` ist der Bereich aller Suffixe, die mit dem Muster beginnen; `hi - lo` die rohe
Trefferzahl auf beiden Strängen.

## Ortung
Von einer Zeile aus mit der LF-Abbildung (`Zeile ← c_before[c] + occ(c, Zeile)`, c = BWT der
Zeile) rückwärts gehen, bis eine Stichprobenzeile oder die Primärzeile erreicht ist. Textposition =
Stichprobe + Zahl der Schritte.

## Von der Textposition zum Treffer
Eine Fundstelle der Länge m an Textposition p liegt auf dem Vorwärtsstrang, wenn p + m ≤ L, und
auf dem Gegenstrang, wenn p ≥ L (linkes Ende auf dem Vorwärtsstrang: 2L − p − m). Verworfen wird
sie, wenn sie über die Grenze zwischen den Strängen oder zwischen zwei Contigs reicht oder eine
N-Strecke berührt. Übrig bleiben die gefilterten Treffer.

## Prüfung
Siehe [Werkzeug index](../werkzeuge/index.md), Abschnitt Validierung. Tests: Treffer gegen naive
Suche in beiden Varianten, `occ` und `occ4` gegen direktes Zählen an jeder Zeile, Datei-Rundreise,
beschädigte und abgeschnittene Dateien, Vollprüfung erkennt eine geänderte Base und eine geänderte
Stichprobe.

## Literatur
Burrows & Wheeler 1994; Ferragina & Manzini 2000; Li & Durbin 2009. Siehe [Literatur](../literatur.md).
