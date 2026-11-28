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

## Optische Duplikate
Nicht jede Kopie stammt aus der PCR. Liegt eine Kopie auf der Fläche des Geräts direkt neben
dem Original, hat das Gerät dasselbe Molekül zweimal gelesen (Illumina: Nachbarcluster;
DNBSEQ: Nachbar-DNB). Solche optischen Duplikate werden genauso markiert, zählen aber bei der
Größe der Bibliothek nicht mit, denn sie sagen nichts über die Vervielfältigung.

Ein Duplikat aus einem Paar gilt als optisch, wenn es mit dem behaltenen Paar derselben Gruppe
Flowcell, Lane und Kachel bzw. Bildfeld teilt und nahe daneben liegt:

- Illumina (`…:LANE:KACHEL:X:Y`): Abstand in x und in y höchstens `--optisch-pixel` (100, wie
  Picard).
- MGI/DNBSEQ (`FLOWCELL L<lane> C<spalte> R<reihe> <laufende Nummer>`): Abstand der laufenden
  Nummer höchstens `--optisch-mgi-abstand` (1000).

Warum bei MGI nicht wie samtools und Genome 3.x: Dort werden `C` und `R` als x und y gelesen und
mit `-d 300` verglichen. `C` und `R` sind aber die Nummern des Bildfelds (bei PJ 2.833 Felder,
`C` 1 bis 36, `R` 1 bis 42, im Mittel rund 126.000 Reads je Feld). Innerhalb einer Lane liegen
alle Felder näher als 300 beieinander, so wurde jedes Duplikat optisch genannt (an einem
Ausschnitt von PJ 3.438 von 3.438). Gemessen an diesem Ausschnitt liegen 47 % der Duplikate im
selben Feld wie ihr Original, zufällig wären es 2,6 %. Bei 79 % von ihnen unterscheidet sich die
laufende Nummer um genau 1, bei 97 % um höchstens 1000; zufällig wären es unter 0,04 %. Die
Nummer gibt also die Lage im Feld wieder.

Der Ausschnitt umfasst nur 40 der 2.833 Felder. Nachbarn im selben Feld liegen immer beide darin,
PCR-Kopien meist nicht; deshalb ist ihr Anteil dort viel zu hoch. Am ganzen Genom von PJ sind
129.700 von 5.391.316 Duplikaten aus Paaren optisch (2,4 %); die geschätzte Bibliothek wächst
dadurch von 21,96 auf 22,50 Milliarden Moleküle (Mac M5, 695 s, 364 MB).

Die Grenze 1000 ist aus der Verteilung am Ausschnitt abgeleitet (1.749 Duplikat-Moleküle,
beobachtet gegen zufällig erwartet im selben Feld):

| Abstand der Nummer höchstens | beobachtet | zufällig erwartet |
|---|---|---|
| 1 | 646 | 0,000 |
| 10 | 651 | 0,000 |
| 100 | 651 | 0,002 |
| 1.000 | 797 | 0,018 |
| 10.000 | 797 | 0,18 |

Zwischen 1.000 und 10.000 kommt nichts hinzu; die übrigen 25 Duplikate im selben Feld liegen im
Bereich des Zufalls.

Das ist eine statistisch begründete Regel, keine gegen eine Wahrheit geprüfte Einstufung: Es gibt
noch keinen Datensatz, in dem für jedes Duplikat bekannt ist, ob es optisch oder aus der PCR
stammt. Bis dahin gilt die Zahl der optischen Duplikate als Schätzung; die markierten Reads hängen
nicht davon ab.

Eingestuft wird erst, wenn die Gruppe das Fenster verlässt: Jedes Duplikat wird mit dem Paar
verglichen, das am Ende behalten wird. So hängt das Ergebnis nicht von der Reihenfolge gleich
liegender Reads ab. Grenze: Picard vergleicht zusätzlich die Duplikate untereinander; bei Gruppen
aus zwei Paaren, dem häufigsten Fall, ist das gleich. Mit `--ohne-optisch` wird nichts eingestuft; dann sind alle Zahlen die von samtools ohne `-d`.

## Geschätzte Größe der Bibliothek
Aus der Zahl der Paare und der doppelten Paare lässt sich schätzen, wie viele verschiedene
Moleküle die Bibliothek enthielt (Lander und Waterman 1988, in der Form von Picard
`DuplicationMetrics`): gesucht ist `X` mit `C / X = 1 - exp(-N / X)`, `N` Paare, `C`
verschiedene Paare. Gelöst wird durch Halbieren des Bereichs in 40 Schritten, wie in samtools.
Wie bei Picard zählen optische Duplikate nicht zu `N`.

## Was die Messung über die Grenzen zeigt
- Position und Richtung sind das einzige Merkmal. Verschiebt eine Lücke nahe am Ende das Ende
  des Reads, wird die Kopie nicht erkannt.
- In Wiederholungen wählt der Aligner unter gleich guten Orten nach der laufenden Nummer des
  Reads. Kopien landen dann an verschiedenen Orten und bleiben unerkannt. Solche Reads haben
  MAPQ 0 und zählen bei der Variantensuche ohnehin nicht.

## Tempo
Das Fenster fragt je Satz zwei Tabellen ab (Paare und Einzelne). Seit 2026-11-28 streuen sie
die Schlüssel mit einer schnellen Multiplikation (Art FxHash) statt mit SipHash, und die
Puffer geschriebener Sätze werden wiederverwendet. Die Tabellen werden nie durchlaufen, die
Ausgabe bleibt Satz für Satz gleich. Profil am Mac (HG002 30x): SipHash war rund 25 % des
Hauptfadens, Anlegen und Freigeben der Satzpuffer rund 10 %. `GC_MARKDUP_SIPHASH=1` schaltet
zurück auf SipHash.

## Literatur
Danecek et al. 2021 (samtools); Lander und Waterman 1988; Picard (Broad Institute),
`MarkDuplicates` und `DuplicationMetrics`. Siehe [Literatur](../literatur.md).
