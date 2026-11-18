# Ausrichten auf mehreren Kernen

Stand: gebaut und geprüft (`gc-align`, `batch.rs`; Prüfbefehl `genome index sam --kerne N`).
Messung: [M3 Kerne und BAM](../../bench/M3-kerne-bam.md).

## Grundsatz
Die Ausgabe darf nicht von der Zahl der Kerne abhängen. Derselbe Aufruf mit denselben Daten gibt
auf jedem Rechner dieselben Zeilen in derselben Reihenfolge. Das ist die Voraussetzung dafür,
dass ein Ergebnis nachgerechnet werden kann.

Bei bwa gilt das nur mit festem `-K`: Ohne diese Option wächst der Block mit der Zahl der Kerne
(10 Millionen Basen je Kern), und weil die Insertgröße je Block geschätzt wird, ändern sich
einzelne Paare. Genome 4.0 bildet die Blöcke immer gleich (`--block-basen`, Standard
10.000.000 Basen).

## Ablauf je Block

Einzelne Reads:

1. Jeder Read wird für sich ausgerichtet (Seeds, Ketten, Bereiche, Rangfolge, Zeilen).
2. Die Zeilen werden in der Reihenfolge der Eingabe geschrieben.

Paare ([Read-Paare](read-paare.md)):

1. Bereiche je Read, auf allen Kernen.
2. Insertgröße aus dem ganzen Block, auf einem Kern. Sie braucht die Bereiche aller Paare
   des Blocks und ist in Millisekunden gerechnet.
3. Je Paar Partnersuche, Bewertung des Paares und Zeilen, auf allen Kernen.
4. Die Zeilen werden in der Reihenfolge der Eingabe geschrieben.

## Verteilung der Arbeit
Die Kerne holen sich die Reads eines Blocks in Portionen von 64 aus einem gemeinsamen Zähler.
Ein fester Anteil je Kern wäre schlechter: Reads aus Wiederholungen kosten ein Vielfaches,
und ein Kern mit einem teuren Abschnitt würde die anderen warten lassen. Jede Portion merkt sich
ihre Stelle im Block; am Ende werden die Portionen nach dieser Stelle geordnet.

Kein Schritt hängt davon ab, welcher Kern einen Read rechnet oder wann:

- Die laufende Nummer eines Reads oder Paares (sie entscheidet zwischen gleich guten Orten)
  kommt aus seiner Stelle in der Eingabe.
- Es gibt keinen gemeinsamen veränderlichen Zustand zwischen den Reads; Index und Optionen
  werden nur gelesen.
- Es gibt keine Zufallszahlen.

Ein Absturz in einem Kern beendet den ganzen Lauf mit der ursprünglichen Meldung; es wird kein
unvollständiges Ergebnis geschrieben, das wie ein fertiges aussieht.

## Protokoll
Am Ende steht `ALIGN.FERTIG` mit Reads, Zeilen, Blöcken, Kernen, Laufzeit und Reads je Sekunde.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Ergebnisse kommen für 1, 2, 3, 8 und 64 Kerne in der Reihenfolge der Eingabe zurück | bestanden (`batch.rs`) |
| Ein Absturz in einem Kern erreicht den Aufrufer | bestanden (`batch.rs`) |
| `genome index sam` mit 1 und 3 Kernen: gleiche Zeilen, einzeln und als Paare, mehrere Blöcke | bestanden (`crates/genome/tests/ausrichten.rs`) |
| chr20, 200.000 Reads und 100.000 Paare mit 1, 2 und 4 Kernen | Ausgabe Byte für Byte gleich, siehe Messung |

## Grenzen
- Gemessen bis 4 Kerne (mehr hat der Messrechner nicht). Wie gut es auf 10 oder mehr Kernen
  skaliert, ist offen.
- Lesen und Entpacken der FASTQ und das Schreiben von SAM-Text laufen auf einem Kern. Bei vielen
  Kernen wird das die Grenze; die BAM-Ausgabe packt auf allen Kernen.
- Ein Block liegt ganz im Arbeitsspeicher (Reads, Bereiche, Zeilen).

## Literatur
Li 2013 (Blöcke und `-K` in bwa), siehe [Literatur](../literatur.md).
