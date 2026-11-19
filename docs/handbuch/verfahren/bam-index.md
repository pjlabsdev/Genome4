# BAM-Index

Stand: gebaut und geprüft (`gc-bam`, `index.rs`; Befehl
[`genome bam-index`](../werkzeuge/bam-index.md)). Messung: [M5 Index](../../bench/M5-index.md).

## Aufgabe
Eine BAM für ein ganzes Genom hat Dutzende Gigabyte. Wer die Reads eines Gens sehen will, darf
nicht die ganze Datei lesen müssen. Der Index sagt zu jedem Bereich des Genoms, an welchen
Stellen der Datei seine Datensätze stehen.

## Stellen in der Datei
Eine BAM besteht aus unabhängig gepackten Blöcken von höchstens 64 kB (BGZF). Eine Stelle heißt
deshalb "virtuelle Position": Dateistelle des Blocks mal 65.536 plus Stelle im entpackten Block.

## Aufbau (SAM-Spezifikation, Abschnitt 5)
Je Contig zwei Teile:

- Fächer. Das Contig ist in sechs Ebenen geteilt: ein Fach für alles (512 Mbp), dann 8 Fächer zu
  64 Mbp, 64 zu 8 Mbp, 512 zu 1 Mbp, 4.096 zu 128 kbp und 32.768 zu 16 kbp. Jeder Datensatz
  gehört in das kleinste Fach, das ihn ganz enthält. Zu jedem Fach stehen die Dateistücke (Anfang
  und Ende als virtuelle Positionen), in denen seine Datensätze liegen.
- Linearer Index. Für jedes Fenster von 16.384 Basen die Stelle des ersten Datensatzes, der das
  Fenster berührt. Damit überspringt eine Abfrage die großen Fächer, soweit ihre Datensätze vor
  dem gesuchten Bereich enden.

Eine Abfrage sammelt die Fächer, die den Bereich schneiden können, liest deren Dateistücke ab
der Stelle aus dem linearen Index und prüft jeden Datensatz selbst.

## Bau
Die Datei wird einmal von vorn nach hinten gelesen. Für jeden Datensatz zählen Contig, Anfang,
Ende auf der Referenz (aus dem CIGAR; ohne Ort oder ohne CIGAR eine Base) und die virtuelle
Position vor und hinter ihm. Der Bau folgt htslib (`hts_idx_push`, `hts_idx_finish`), damit das
Ergebnis mit dem von samtools verglichen werden kann:

1. Aufeinanderfolgende Datensätze desselben Fachs bilden ein Dateistück.
2. Der lineare Index merkt sich je Fenster den ersten Datensatz; leere Fenster zeigen auf den
   nächsten dahinter.
3. Am Ende gehen kleine Fächer im übergeordneten auf, wenn ihre Stücke weniger als 65.536
   gepackte Bytes überspannen und das übergeordnete Fach schon Einträge hat. Das spart Sprünge
   in der Datei.
4. Stücke eines Fachs, die im selben gepackten Block aneinanderstoßen, werden verbunden.
5. Ein Sonderfach (Nummer 37450) trägt je Contig Anfang und Ende seiner Daten und die Zahl der
   Reads mit und ohne Ort. Daraus liest `samtools idxstats` seine Zahlen.
6. Am Schluss steht die Zahl der Datensätze ohne Contig.

Geprüft wird beim Bau, dass die Datei wirklich sortiert ist: Positionen steigen, jedes Contig
kommt in einem Stück, Datensätze ohne Contig stehen am Ende. Sonst bricht der Bau ab.

## Literatur
SAM-Spezifikation, Abschnitt 5 (Indexing BAM); Li et al. 2009; Kent et al. 2002 (das Schema der
Fächer). Siehe [Literatur](../literatur.md).
