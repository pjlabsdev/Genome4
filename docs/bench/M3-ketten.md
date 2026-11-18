# M3, Teil 2: Seeds orten und verketten

Stand: 18.11.2026. Verfahren: [Seeds und Ketten](../handbuch/verfahren/seeds-und-ketten.md).

## Aufbau
- Referenz und Reads wie in [M3 Seed-Suche](M3-seed-suche.md): chr20, 200.000 erzeugte Reads zu
  150 Basen, 1 % Austauschfehler, wahrer Ort im Namen.
- Befehl: `genome index ketten chr20.kompakt.gcfm --fastq sim150.fq --wahrheit-im-namen [--max-orte N] [--fehlende datei]`
- Eine Kette gilt als am wahren Ort, wenn Contig und Strang stimmen und das aus ihrem ersten Seed
  gerechnete linke Ende des Reads höchstens 100 Basen vom wahren abweicht.
- Rechner: Cloud-Container, Intel Xeon 2,3 GHz, ein Kern. Ein Lauf je Zeile.

## Ergebnis

| | höchstens 500 Orte je Treffer (Vorgabe) | ohne Ausdünnen |
|---|---|---|
| Seeds | 7.844.037 | 8.663.704 |
| Ketten | 4.848.788 | 5.391.867 |
| Reads ohne Kette | 0 | 0 |
| Reads mit einem Seed genau am wahren Ort | 199.926 | 199.961 |
| Reads mit einer Kette am wahren Ort | 199.909 (99,95 %) | 199.944 |
| schwerste Kette am wahren Ort, eindeutig | 192.211 (96,1 %) | 192.261 |
| schwerste Kette am wahren Ort, Gleichstand mit anderen | 6.442 (3,2 %) | 6.863 |
| Zeit | 43,5 s (4.594 Reads je Sekunde) | 45,3 s |

## Die Reads ohne Seed am wahren Ort
Mit der Vorgabe sind es 74 von 200.000. Jeder Read hat am wahren Ort einen fehlerfreien Abschnitt
von mindestens 19 Basen (nachgezählt), die Seed-Suche müsste ihn also grundsätzlich finden
können. `bench/m3/explain_missing.py` prüft jeden Fall mit der einfachen Zählung im Index:

- **40 Reads:** Jeder fehlerfreie Abschnitt lässt sich im Read über die fehlerhafte Nachbarbase
  hinaus verlängern und kommt so anderswo im Genom vor. Der Abschnitt am wahren Ort ist dann kein
  längster Treffer des Reads; gemeldet wird der längere Treffer, der woanders liegt. Das ist das
  Verhalten des Verfahrens in Wiederholungen, kein Fehler der Umsetzung.
- **34 Reads:** Der Abschnitt kommt öfter als 500-mal vor, und das Ausdünnen hat den wahren Ort
  nicht erwischt.
- **0 Reads** unerklärt.

Weitere 17 Reads haben einen Seed am wahren Ort, aber nach dieser Zählung keine Kette dort. Die
Ursache ist nicht einzeln geprüft. Naheliegend ist, dass der Seed als enthalten in eine Kette
einging, deren erster Seed auf einer anderen Diagonale liegt (Tandem-Wiederholungen); das ist
eine Vermutung.

## Einordnung
- Bei 1.256 Reads (0,63 %) liegt die schwerste Kette nicht am wahren Ort, obwohl es dort eine
  Kette gibt. Das Gewicht der Kette ist noch kein Alignment-Score; die Entscheidung fällt erst im
  Alignment. Die Zahl ist deshalb keine Fehlerrate des fertigen Programms.
- Wie bwa bei denselben Reads abschneidet, ist nicht gemessen.
- Die dritte Runde hat die Zahl der Reads mit Seed am wahren Ort in diesem Datensatz nicht
  verändert (199.926 mit und ohne), aber 2,0 Mio. Seeds hinzugefügt und die Zeit von 32,8 s auf
  rund 44 s erhöht. Sie bleibt an, weil bwa sie in der Vorgabe nutzt und der Vergleich mit bwa sonst
  nicht fair wäre.

Rohdaten: `bench/m3/ergebnisse/ketten_*.json`, `reads_ohne_seed_am_wahren_ort.txt`.
