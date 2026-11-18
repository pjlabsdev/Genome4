# Blockweiser Index-Bau für das ganze Genom

Stand: gebaut (`gc-index`, `blockwise.rs`). Messung:
[Index ganzes Genom](../../bench/M2-index-ganzes-genom.md).

## Warum
Beide Stränge des Humangenoms sind 6,2 Mrd. Zeichen. Das volle Suffix-Array bräuchte mit
40-Bit-Positionen rund 31 GB RAM. Ziel ist ein Bau mit festem, einstellbarem RAM-Bedarf, der auch
auf einem Rechner mit 16 GB läuft.

## Idee (Kärkkäinen 2007)
Die Suffixe werden nach ihren ersten 12 Basen in 16,8 Mio. Fächer geteilt. Benachbarte Fächer
werden zu Blöcken von höchstens `--block-mio` Millionen Suffixen zusammengefasst. Jeder Block wird
für sich sortiert. Weil die Blöcke in lexikographischer Reihenfolge aufeinander folgen, entsteht
die BWT fortlaufend, ohne dass je das ganze Suffix-Array im Speicher liegt.

## Das Problem langer Wiederholungen
Zwei Suffixe direkt zu vergleichen dauert so lange wie ihr gemeinsamer Anfang. In
Segmentduplikationen und Tandem-Wiederholungen ist der Zehntausende Basen lang. Lösung ist eine
**Differenzüberdeckung** (Burkhardt und Kärkkäinen 2003):

- D ist eine Menge von 40 Resten modulo v = 1024 mit der Eigenschaft: zu jedem Abstand r gibt es
  ein a in D, sodass auch a + r (mod v) in D liegt. Konstruktion nach Colbourn und Ling 2000; die
  Eigenschaft wird beim Bau nachgeprüft.
- **Stichprobe:** alle Textpositionen, deren Rest modulo v in D liegt (3,9 % der Positionen).
  Ihre Suffixe werden vorab vollständig sortiert; jede bekommt einen Rang.
- **Vergleich zweier beliebiger Suffixe i und j:** höchstens v Basen direkt (32 Basen je
  Wortvergleich). Sind sie dann gleich, gibt es einen Versatz d < v, bei dem i + d und j + d beide
  in der Stichprobe liegen; deren Ränge entscheiden.

Damit kostet jeder Vergleich höchstens 32 Wortvergleiche und einen Blick in die Rangliste, gleich
wie lang die Wiederholung ist.

## Schritte

1. **Text:** beide Stränge 2-Bit-gepackt, erste Base in den höchsten Bits, sodass ein Vergleich
   zweier Wörter 32 Basen lexikographisch ordnet (1,55 GB beim Humangenom).
2. **Stichprobe sortieren:** nach den ersten 12 Basen in Fächer, jedes Fach nach den ersten v
   Basen sortieren. Gruppen mit gleichen ersten v Basen werden durch Verdoppeln getrennt (Larsson
   und Sadakane): der Rang der Position v weiter, die wieder zur Stichprobe gehört, trennt nach 2v
   Basen, dann 4v und so weiter, bis jede Position einen eigenen Rang hat (0,97 GB).
3. **Fächer zählen** über alle Suffixe, Blöcke bilden.
4. **Je Block:** Text einmal durchgehen und alle Suffixe des Blocks einsammeln (16 Bytes je
   Suffix: die nächsten 32 Basen als Schlüssel, Position, Base davor); jedes Fach sortieren nach
   Schlüssel, dann direktem Vergleich, dann Rangliste; Zeilen an den Index anfügen.
5. Zeile 0 ist das Suffix, das nur aus dem Endzeichen besteht.

Endet ein Suffix vor der Vergleichstiefe, ist das kürzere das kleinere (das Endzeichen ist kleiner
als jede Base). Fächer und Schlüssel füllen hinter dem Textende mit A auf; der direkte Vergleich
stellt die richtige Reihenfolge her.

## Mehrere Kerne
Zählen, Einsammeln und Sortieren laufen über alle Kerne. Beim Einsammeln zählt jeder Kern erst,
wie viele seiner Positionen in jedes Fach fallen; daraus bekommt er feste Plätze im Block und
schreibt dann ohne Abstimmung. Das Ergebnis hängt nicht von der Kernzahl ab, weil alle Suffixe
verschieden sind und die Sortierung deshalb eindeutig ist.

## RAM-Bedarf (Humangenom, Variante kompakt)

| Teil | Größe |
|---|---|
| Text | 1,55 GB |
| Rangliste der Stichprobe | 0,97 GB (beim Sortieren kurzzeitig das Doppelte) |
| Block | 16 MB je Million Suffixe, bei `--block-mio 96` 1,5 GB |
| Index im Aufbau | bis 4,8 GB |
| Referenz (gepackt) | 0,78 GB |

Gemessener Spitzenwert: siehe Messdokument. Der Index wird noch im RAM aufgebaut; ein Bau, der
fortlaufend auf die Platte schreibt, würde den Bedarf weiter senken (offen).

## Prüfung

- **Tests:** in 60 Zufallsfällen (bis drei Contigs, exakte Wiederholungen auch als reverses
  Komplement, Tandems mit kurzer Periode, Homopolymere, N) liefert der blockweise Bau denselben
  Index wie der direkte, mit kleinen Stellgrößen (v ab 4, Präfix ab 1 Base, Blöcke ab 1 Suffix,
  1 bis 4 Kerne), damit jeder Weg durchlaufen wird. Die Überdeckung wird für v von 4 bis 2048 geprüft.
- **chr20:** blockweise und direkt gebaute Datei sind Byte für Byte gleich.
- **Ganzes Genom:** Stichprobenprüfung und Vollprüfung, siehe Messdokument.

## Literatur
Kärkkäinen 2007; Burkhardt & Kärkkäinen 2003; Colbourn & Ling 2000; Larsson & Sadakane 2007.
Siehe [Literatur](../literatur.md).
