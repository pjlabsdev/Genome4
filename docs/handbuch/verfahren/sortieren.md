# Sortieren nach Ort mit festem Speicher

Stand: gebaut und geprüft (`gc-bam`, `sort.rs` und `header.rs`; Befehl
[`genome sort`](../werkzeuge/sort.md)). Messung: [M5 Sortieren](../../bench/M5-sortieren.md).

## Aufgabe
Der Aligner schreibt die Datensätze in der Reihenfolge der Reads. Alles Weitere (Duplikate
markieren, Index, Varianten suchen) braucht sie nach dem Ort geordnet. Ein ganzes Genom mit
30-facher Abdeckung hat rund 600 bis 900 Millionen Datensätze und unkomprimiert weit über
100 GB. Das passt nicht in den Speicher eines Arbeitsplatzrechners, also wird in Teilen sortiert
und zusammengeführt (externes Mischsortieren, Knuth 1998, Abschnitt 5.4).

## Reihenfolge
Dieselbe wie bei `samtools sort` ohne Optionen:

1. Contig in der Reihenfolge des Kopfes. Datensätze ohne Contig (`refID` = -1) stehen am Ende.
2. Position (linkes Ende, wie im Feld `pos`).
3. Vorwärtsstrang vor Gegenstrang (Bit 16 im Feld `flag`).
4. Bei völligem Gleichstand bleibt die Reihenfolge der Eingabe (stabil).

Ein Read ohne eigenen Ort, dessen Partner einen Ort hat, trägt nach der SAM-Spezifikation Contig
und Position des Partners und wird dort einsortiert. Das ist gewollt: So stehen beide beieinander.

Die drei Kriterien stehen in einer Zahl von 64 Bit: 31 Bit Contig (ohne Contig der größte Wert),
32 Bit Position plus 1, 1 Bit Strang. Verglichen wird nur diese Zahl.

## Ablauf
1. Die Datensätze werden roh gelesen (so wie sie in der BAM stehen, mit Längenfeld) und nicht
   in Felder zerlegt. Gelesen werden nur Contig, Position und Flag an ihren festen Stellen.
   Damit kann das Sortieren keinen Inhalt verändern.
2. Sie sammeln sich in einem Puffer. Zu jedem merkt sich eine Liste Schlüssel und Stelle im Puffer
   (16 Bytes).
3. Würde der nächste Datensatz die Grenze `--speicher-mb` überschreiten, wird die Liste sortiert
   (nach Schlüssel, dann Stelle, also stabil) und der Inhalt in dieser Reihenfolge als Teil auf
   die Platte geschrieben ([Format](../formate/gcsort.md)). Der Puffer ist danach leer.
4. Am Ende wird der Rest im Speicher sortiert. Gibt es keine Teile, wird er direkt geschrieben.
   Sonst laufen alle Teile und der Rest gemeinsam durch eine Vorrangschlange: Es kommt immer der
   Datensatz mit dem kleinsten Schlüssel, bei Gleichstand der aus dem früheren Teil. Weil die
   Teile in der Reihenfolge der Eingabe entstanden sind, bleibt das Ganze stabil.

Folge: Das Ergebnis hängt nicht von `--speicher-mb`, nicht von der Zahl der Teile und nicht von
der Zahl der Kerne ab.

## Sicherheit der Daten
- Jeder Teil endet mit Zahl der Datensätze, Zahl der Bytes und Prüfsumme (CRC-32). Beim
  Zusammenführen wird alles drei geprüft. Ein abgeschnittener oder veränderter Teil führt zum
  Abbruch und nicht zu einer lückenhaften Ausgabe.
- Die Ausgabe trägt bis zum letzten Byte den Namen `<ziel>.unfertig`.
- Die Teile liegen in einem eigenen Unterordner und werden am Ende und im Fehlerfall gelöscht.

## Kopf
`@HD` bekommt `SO:coordinate`; ein vorhandenes Feld `SS` (Untersortierung) wird entfernt, weil es
die alte Ordnung beschreibt. Eine neue `@PG`-Zeile nennt Programm, Version und Aufruf und verweist
mit `PP` auf die letzte vorhandene. Der Rest des Kopfes bleibt Byte für Byte erhalten.

## Abweichungen von samtools
Keine in der Reihenfolge (gemessen, siehe Messung). Unterschiede im Drumherum:

- samtools legt die Teile als BAM ab, Genome 4.0 in einem eigenen, einfacheren Format mit
  Prüfsumme über den ganzen Teil.
- samtools verteilt den Speicher je Kern (`-m` je Thread), hier gilt eine Grenze für den Lauf.

## Literatur
Knuth 1998 (Sorting and Searching); Li et al. 2009 (SAM/BAM und samtools); SAM-Spezifikation,
Abschnitt 1.3 (Sortierung) und 4.2 (BAM-Datensatz). Siehe [Literatur](../literatur.md).
