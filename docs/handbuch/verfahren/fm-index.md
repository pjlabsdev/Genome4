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
| standard (Vorgabe) | 128 | 8 | 20 |
| kompakt | 128 | 32 | 12,5 |
| schnell | 64 | 8 | 24 |

Die Zahl gilt für beide Stränge samt gepackter Referenz. Rechnung für kompakt, je Zeile der BWT:
Zählblock (256 + 256 Bit auf 128 Zeilen) 4 Bit, Stichproben 40 / 32 = 1,25 Bit; mal zwei Zeilen
je Referenzbase, plus 2 Bit Referenz. Bei standard sind es 40 / 8 = 5 Bit für die Stichproben,
zusammen 20 Bit.

Am ganzen Genom (GRCh38 ohne ALT) ist die Datei mit standard 7,75 GB groß, mit kompakt 4,84 GB
und mit schnell 9,3 GB (gerechnet). Der Index liegt beim Ausrichten ganz im Arbeitsspeicher.
`genome fq2bam` braucht am ganzen Genom damit in der Spitze 8,07 GB statt 5,24 GB (gemessen).
Standard ist die Vorgabe, weil die dichteren Stichproben das Ausrichten am ganzen Genom auf rund
0,87-mal die Zeit von kompakt bringen (siehe [M6 Tempo](../../bench/M6-tempo.md)); die Ausgabe
ist mit jeder Variante dieselbe. Wer wenig Arbeitsspeicher hat, baut mit `--variante kompakt`.
Jede Datei trägt ihre Abstände im Kopf, ältere Indexdateien bleiben also lesbar.

### Ausdünnen beim Laden

`FmIndex::read_thinned(datei, n)` behält beim Laden nur jede n-te Stichprobe (n eine
Zweierpotenz, der Abstand wächst höchstens bis 1024). Die Stichprobe zur Zeile `r` steht an Stelle
`r / s`; behält man jede n-te, steht die Stichprobe zur Zeile `r` mit `r` durch `n·s` teilbar an
Stelle `r / (n·s)`. Die Textpositionen sind also dieselben, nur das Rückwärtsgehen bis zur
nächsten Stichprobe dauert im Mittel n-mal so lange. Gelesen wird stückweise, die vollen
Stichproben liegen nie zugleich im Speicher. Am ganzen GRCh38 (standard, 421.564 echte Paare,
4 Kerne) brauchte das Ausrichten 31,6 s mit allen Stichproben, 32,4 s mit n = 2 und 38,9 s mit
n = 4; die BAM war in allen drei Fällen gleich (ohne `@PG`). Mit n = 2 sinkt der Speicher für die
Stichproben von 4,0 auf 2,0 GB. Das hilft auf Rechnern mit 16 GB, auf denen der volle Index samt
Sortierpuffer sonst ausgelagert wird (`genome fq2bam --sa-ausduennen`).

## Suche
`occ(c, i)` zählt die Base c in den Zeilen 0 bis i-1 der BWT: Zählstand des Blocks plus Auszählen
der Wörter bis i. `occ4(i)` liefert die Zählung für alle vier Basen in einem Durchgang durch
denselben Block (wie `bwt_occ4` in bwa): gezählt werden C, G und T, A ergibt sich als Rest. Die
Seed-Suche braucht immer alle vier, und immer an beiden Enden eines Bereichs: `occ4_pair(k, l)`
(wie `bwt_2occ4` in bwa) zählt, wenn beide Enden im selben Wort der BWT liegen, nur bis `k` und
das kurze Stück bis `l` dazu. Nach wenigen Basen sind die Bereiche der Seed-Suche eng, das ist
dann der häufigste Fall. Liegt `i` in der hinteren Hälfte eines Blocks, zählt `occ4` vom
Zählstand des nächsten Blocks rückwärts (Zählstand minus die Basen von `i` bis zum Blockende,
das Endzeichen wie vorn berücksichtigt); so sind es höchstens halb so viele Basen. Der letzte
Block zählt immer vorwärts. Getestet an jeder Zeile in allen drei Varianten.

**Zählen.** Je Base wird ein Bit gebildet (das untere ihrer zwei), gesetzt, wo die gesuchte Base
steht. Die Wörter eines Blocks werden nicht einzeln fertig gezählt: Jedes Wort wird nur zu Summen
je Halbbyte verdichtet, die Summen werden addiert, und erst am Ende wird einmal ganz
zusammengezählt. Je Halbbyte kommen höchstens 2 dazu, sieben Wörter passen also in eine Summe;
ein Block hat höchstens vier. Das braucht keinen besonderen Befehl des Prozessors und läuft auf
jedem Rechner gleich.

**Zählen bei Blöcken zu 128 Basen.** Für den häufigsten Abstand (`b` = 128, Varianten
`standard` und `kompakt`) hat `occ4` einen eigenen Weg: Vorn wie hinten werden genau zwei Wörter
mit einer Maske gezählt, ob vorwärts oder rückwärts, ist eine Auswahl von Werten statt eines
Sprungs (welche Blockhälfte, wechselt von Zugriff zu Zugriff und ist kaum vorherzusagen). Wo der
Prozessor einen Befehl zum Zählen gesetzter Bits hat (alle 64-Bit-ARM wie Apple M, x86 mit
`popcnt`), zählt dieser statt der Halbbyte-Summen. Ohne den Befehl bringt der Weg nichts, mit
ihm braucht die Seed-Suche allein in der Cloud (x86 mit `-C target-cpu=native`, 200.000
HG002-Reads, ganzes GRCh38, ein Kern, abwechselnd je dreimal) 7,5 bis 7,8 s statt 8,8 bis 8,9 s,
mit denselben Treffern. Schalter zurück: `FmIndex::set_occ_plain(true)`; Messprogramm
`examples/seedbench` mit `OCC_PLAIN=1`. Der letzte Block zählt hinten weiter mit dem alten Weg.

**Abstände.** `b` und `s` sind Zweierpotenzen, `b` höchstens 128. Block und Stelle im Block
ergeben sich durch Schieben und Maskieren. Eine echte Division je Zugriff kostete mehr als das
Zählen selbst ([M6 Tempo](../../bench/M6-tempo.md)).

Rückwärtssuche über das Muster von hinten nach vorn:

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

`locate_many` ortet viele Zeilen auf einmal, mit demselben Ergebnis wie `locate` je Zeile. Die
Wege von bis zu 32 Zeilen laufen im Gleichschritt: In jeder Runde wird erst für alle offenen
Zeilen die Base der BWT gelesen, dann für alle gezählt. Die Zugriffe einer Runde hängen nicht
voneinander ab, der Prozessor wartet also auf mehrere Stellen des Arbeitsspeichers zugleich statt
auf eine nach der anderen. Beim ganzen Genom liegt fast jeder Schritt außerhalb des
Zwischenspeichers. Der Aligner ortet so alle Seeds eines Reads gemeinsam.

### Ortsspeicher
Reads aus Wiederholungen (Alu, L1) finden ihre Seeds immer wieder in denselben Zeilen: Jeder
Alu-Read ortet bis zu 500 Fundstellen je Seed, und es sind meist dieselben. `locate_many` merkt
sich deshalb Zeile und Position in einem Speicher mit `2^b` Plätzen (`fq2bam --ortsspeicher-mb`,
Vorgabe 64 MB, 0 = aus). Platz = untere `b` Bits der Zeile; im Platz steht in einem Wort ein
Belegt-Bit, die oberen Bits der Zeile und 33 Bit Position. Ein Treffer liefert genau den Wert, den
der Weg durch die BWT ergäbe; ein neuer Wert überschreibt den alten Platz. Alle Kerne teilen sich
den Speicher und lesen und schreiben ohne Sperre, denn ein Platz ist ein einziges atomares Wort.
Geht nur bis 2^33 Zeilen (Referenz bis 4,2 Mrd. Basen), sonst bleibt er aus.

| Messung (421.564 HG002-Paare, GRCh38, 4 Kerne, `--schnellweg an`) | ohne | 16 MB | 64 MB | 128 MB |
|---|---|---|---|---|
| Orten je Read aus Wiederholung | 121–124 µs | 75 µs | 75–79 µs | 81 µs |
| Ausrichten gesamt | 25,3–26,5 s | 23,0 s | 22,4–24,0 s | 23,7 s |

Die BAM ist byte-gleich mit und ohne Speicher (844.609 von 844.609 Sätzen).

## Von der Textposition zum Treffer
Eine Fundstelle der Länge m an Textposition p liegt auf dem Vorwärtsstrang, wenn p + m ≤ L, und
auf dem Gegenstrang, wenn p ≥ L (linkes Ende auf dem Vorwärtsstrang: 2L − p − m). Verworfen wird
sie, wenn sie über die Grenze zwischen den Strängen oder zwischen zwei Contigs reicht oder eine
N-Strecke berührt. Übrig bleiben die gefilterten Treffer.

## Referenzstücke lesen
Für das Verlängern und die Partnersuche holt `text_range` Stücke der Referenz, auf dem
Rückwärtsstrang als reverses Komplement, N als Code 4. Die Referenz liegt 2-Bit-gepackt, 32 Basen
je Wort. `Packed::extend_range` entpackt Wort für Wort statt jede Base einzeln zu suchen; ein
Stück des Rückwärtsstrangs wird als Vorwärtsstück gelesen, umgedreht und komplementiert. Am Mac
kostete das Lesen Base für Base 3,3 % der Rechenzeit von `fq2bam`. Die Stücke sind dieselben
(Tests `extend_range_equals_get` und `text_range_wordwise_equals_base_by_base` über alle
Bereiche, auch über Wort- und Stranggrenzen), die BAM der Simulation (300.000 Paare) hat
dieselben Sätze. `GC_TEXT_WORTWEISE=0` liest wie bisher Base für Base.

## Prüfung
Siehe [Werkzeug index](../werkzeuge/index.md), Abschnitt Validierung. Tests: Treffer gegen naive
Suche in beiden Varianten, `occ` und `occ4` (auf beiden Wegen) gegen direktes Zählen an jeder Zeile, `occ4_pair` gegen
zweimal `occ4` für jeden Anfang und 140 Längen, `locate_many` gegen `locate` für alle Zeilen und
für zufällige Mengen mit Wiederholungen, Datei-Rundreise,
beschädigte und abgeschnittene Dateien, Vollprüfung erkennt eine geänderte Base und eine geänderte
Stichprobe.

## Literatur
Burrows & Wheeler 1994; Ferragina & Manzini 2000; Li & Durbin 2009. Siehe [Literatur](../literatur.md).
