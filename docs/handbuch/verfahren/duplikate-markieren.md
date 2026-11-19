# Duplikate markieren

Stand: gebaut und geprüft (`gc-bam`, `markdup.rs`; Befehl
[`genome markdup`](../werkzeuge/markdup.md)). Messung: [M5 Duplikate](../../bench/M5-duplikate.md).

## Aufgabe
Bei der Vorbereitung der Probe werden die Moleküle vervielfältigt. Liest das Gerät zwei Kopien
desselben Moleküls, sind das keine zwei unabhängigen Beobachtungen. Ein Lesefehler in der
Vervielfältigung sähe sonst aus wie eine Variante mit zwei Belegen. Deshalb bleibt von jeder
Gruppe von Kopien eine Lesung, die anderen bekommen Flag 1024.

Kopien erkennt man am Ort: Zwei Paare aus demselben Molekül beginnen und enden an denselben
Stellen des Genoms. Dass zwei unabhängige Moleküle das zufällig auch tun, kommt vor und wird in
Kauf genommen.

## Vorgehen
Das Verfahren ist das von `samtools markdup` (Danecek et al. 2021) in der Vorgabe, Schritt für
Schritt nachgebaut, damit jede Entscheidung mit samtools verglichen werden kann.

### Was verglichen wird
Für jeden Read zählt das Ende, an dem das Lesen begann, einschließlich abgeschnittener Basen
(`S` und `H` im CIGAR): auf dem Vorwärtsstrang der Anfang, auf dem Gegenstrang das rechte Ende.
So bleibt der Vergleich richtig, wenn bei einer Kopie die ersten Basen nicht passten.

- Read mit Partner: Schlüssel aus eigenem Contig und Ende, Contig und Ende des Partners,
  Richtung des Paares (FR, RF, FF, RR) und der Angabe, ob dieser Read der linke ist. Den Ort
  des Partners kennt der Read aus seinen eigenen Feldern und aus `MC` (CIGAR des Partners).
- Read ohne Partner (auch: Partner ohne Ort): Schlüssel aus Contig, Ende und Strang.

### Wer bleibt
- Unter Paaren mit gleichem Schlüssel bleibt das mit der höchsten Summe der Basenqualitäten ab
  15 über beide Reads. Die Summe des Partners steht in `ms`. Bei Gleichstand bleibt das Paar,
  dessen Name im Alphabet früher steht. Beide Reads eines Paares kommen so für sich zum selben
  Urteil, ohne dass sie einander sehen.
- Liegt ein Read ohne Partner am selben Ende wie ein Read eines Paares, ist der ohne Partner
  das Duplikat, gleich wie gut er ist.
- Unter Reads ohne Partner bleibt der mit der höheren Qualitätssumme, bei Gleichstand der
  frühere.

### Was nicht bewertet wird
Reads ohne Ort, zusätzliche Zeilen (Flags 256 und 2048) und Reads mit Flag 512 werden weder
verglichen noch markiert. Das entspricht samtools ohne `-S`, wie Genome 3.x es aufruft.

### Fenster
Die Datei wird einmal von vorn nach hinten gelesen. Ein Datensatz bleibt im Speicher, bis die
Eingabe mehr als `--max-laenge` Basen hinter seinem Ende ist; danach kann kein gleicher mehr
kommen, und er wird geschrieben. Der Speicher hängt deshalb von der Abdeckung ab und nicht von
der Größe der Datei (gemessen: höchstens 26 Datensätze gleichzeitig bei rund 0,6-facher
Abdeckung).

Folge der Bauart: Die Reihenfolge der Datensätze ändert sich nicht, und am Datensatz ändert sich
nur ein Bit im Flag.

## Geschätzte Größe der Bibliothek
Aus der Zahl der Paare und der doppelten Paare lässt sich schätzen, wie viele verschiedene
Moleküle die Bibliothek enthielt (Lander und Waterman 1988, in der Form von Picard
`DuplicationMetrics`): gesucht ist `X` mit `C / X = 1 - exp(-N / X)`, `N` Paare, `C`
verschiedene Paare. Gelöst wird durch Halbieren des Bereichs in 40 Schritten, wie in samtools.

## Was die Messung über die Grenzen zeigt
- Position und Richtung sind das einzige Merkmal. Verschiebt eine Lücke nahe am Ende das Ende
  des Reads, wird die Kopie nicht erkannt.
- In Wiederholungen wählt der Aligner unter gleich guten Orten nach der laufenden Nummer des
  Reads. Kopien landen dann an verschiedenen Orten und bleiben unerkannt. Solche Reads haben
  MAPQ 0 und zählen bei der Variantensuche ohnehin nicht.

## Literatur
Danecek et al. 2021 (samtools); Lander und Waterman 1988; Picard (Broad Institute),
`MarkDuplicates` und `DuplicationMetrics`. Siehe [Literatur](../literatur.md).
