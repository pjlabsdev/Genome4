# Read-Paare

Stand: gebaut und geprüft (`gc-align`, `pair.rs` und `local.rs`; Prüfbefehl
`genome index sam --fastq R1 --fastq2 R2`). Messung: [M3 Paare gegen bwa](../../bench/M3-paare.md).
Mehrere Kerne: [Ausrichten auf mehreren Kernen](mehrere-kerne.md). Noch nicht gebaut:
ALT-Contigs.

## Worum es geht
Beim Sequenzieren wird ein DNA-Stück von beiden Enden gelesen. Read 1 und Read 2 liegen im Genom
deshalb in bekannter Richtung zueinander (bei Illumina und MGI: aufeinander zu) und in einem
Abstand, der um die Stücklänge streut. Das nutzt das Ausrichten dreifach:

1. Ein Read, der allein an mehrere Stellen passt, bekommt durch seinen Partner einen Ort.
2. Ein Read, der allein keinen Ort findet (zu viele Fehler für einen Seed), wird neben seinem
   Partner gezielt gesucht.
3. Paare, die nicht zusammenpassen, bleiben erkennbar. Das ist später der Hinweis auf
   Strukturvarianten.

Das Verfahren folgt `bwamem_pair.c` aus BWA-MEM (Li 2013), Schritt für Schritt.

## Ablauf für einen Block von Paaren

### 1. Bereiche je Read
Jeder Read wird zuerst allein behandelt wie bisher: Seeds, Ketten, Bereiche, Bereinigung
([Bereiche bewerten](bereiche-bewerten.md)). Die Rangfolge wird noch nicht festgelegt.

### 2. Insertgröße schätzen (`mem_pestat`)
Aus den Paaren des Blocks, bei denen beide Reads einen deutlich besten Ort auf demselben Contig
haben (der zweitbeste Score höchstens 80 % des besten), wird der Abstand der Read-Anfänge
gesammelt, getrennt nach vier Richtungen (FF, FR, RF, RR; F = vorwärts, R = Gegenstrang, Read
mit dem kleineren Ort zuerst). Abstände über 10.000 zählen nicht.

Je Richtung mit mindestens 10 Paaren:

- 25., 50. und 75. Perzentil (p25, p50, p75).
- Mittel und Streuung nur aus den Werten zwischen p25 − 2·(p75 − p25) und p75 + 2·(p75 − p25),
  damit Ausreißer nicht mitzählen.
- Ein Paar gilt als **richtiges Paar**, wenn sein Abstand zwischen p25 − 3·(p75 − p25) und
  p75 + 3·(p75 − p25) liegt; die Grenzen werden auf mindestens Mittel ± 4 Streuungen geweitet.

Richtungen mit weniger als 10 Paaren oder weniger als 5 % der häufigsten Richtung werden nicht
benutzt. Bei üblichen Daten bleibt nur FR.

Die Schätzung steht je Block im Protokoll (`ALIGN.INSERTGROESSE`).

### 3. Partnersuche (`mem_matesw`)
Für jeden Bereich eines Reads, der höchstens 17 Punkte unter dem besten liegt (höchstens 50
Bereiche), wird geprüft, ob der Partner schon einen Bereich im erlaubten Abstand hat. Wenn nicht,
wird der Partner im Fenster, in dem er liegen müsste, lokal ausgerichtet
(Smith-Waterman, `ksw_align2` in bwa, hier `local.rs`). Ein Fund mit Score ab 19 kommt zu den
Bereichen des Partners dazu. Er merkt sich auch den zweitbesten Score im Fenster (`csub`): Liegt
dort eine fast gleich gute zweite Stelle (dicht benachbarte Wiederholung), begrenzt das später
die MAPQ.

### 4. Rangfolge je Read
Wie bei einzelnen Reads (`mem_mark_primary_se`), jetzt einschließlich der neu gefundenen Bereiche.

### 5. Das beste Paar (`mem_pair`)
Alle Bereiche beider Reads werden nach Ort geordnet. Jede Kombination aus einem Bereich von
Read 1 und einem von Read 2 in erlaubter Richtung und erlaubtem Abstand bekommt einen Score:

    Score(Paar) = Score(Read 1) + Score(Read 2) + 0,721 · ln(2 · erfc(|z| / √2))

mit z = (Abstand − Mittel) / Streuung. Der letzte Summand ist negativ: die Wahrscheinlichkeit,
dass ein Stück mindestens so weit vom Mittel abweicht, umgerechnet in Score-Punkte
(0,721 = 1 / ln 4). Das beste Paar gewinnt. Bei Gleichstand entscheidet ein Streuwert aus der
laufenden Nummer des Paares, damit die Wahl nicht von der Lage im Genom abhängt und
wiederholbar ist.

### 6. Paar oder Einzelorte
Das Paar wird mit den beiden besten Einzelorten verglichen, denen 17 Punkte abgezogen werden
(`-U`). Ist das Paar besser, werden seine Orte gemeldet, sonst die besten Einzelorte.

MAPQ, wenn das Paar gewinnt:

- MAPQ des Paares: 6,02 · (Score des Paares − Score der nächstbesten Möglichkeit), vermindert bei
  mehreren fast gleich guten Paaren, höchstens 60, gedämpft um den Anteil sehr häufiger Seeds.
  Die nächstbeste Möglichkeit ist das zweitbeste Paar oder die Einzelorte mit Abzug.
- Je Read gilt die größere von eigener MAPQ und MAPQ des Paares, aber höchstens eigene MAPQ + 40.
- Obergrenze 6,02 · (Score − `csub`) aus der Partnersuche.

Hat ein Read mehrere Teile an verschiedenen Orten (geteilter Read), wird das Paar wie zwei
einzelne Reads gemeldet, mit den Feldern des Partners.

### Partnersuche gesammelt rechnen
Vor dem Paaren plant `planned_rescues` für jedes Paar die Fenster, die die Partnersuche mit den
Bereichen vor der Suche ausrichten würde. Alle Fenster des Blocks werden dann gemeinsam gerechnet
(heute auf allen Kernen, später auch auf der Grafikeinheit), und `pair_lines_with` nimmt beim
Paaren die fertigen Ergebnisse. Ändert die Suche selbst die Bereiche so, dass ein Fenster
dazukommt, wird es einzeln gerechnet; ein überflüssiges bleibt ungenutzt. Die Ausgabe ist in jedem
Fall dieselbe. Auf 421.564 HG002-Paaren (GRCh38) waren 218.556 Fenster geplant und 218.563
gebraucht, die BAM war byte-gleich, das Tempo gleich. `PairOptions::batch_rescue = false` schaltet
das Sammeln ab.

**Je zwei Fenster verzahnt.** Die gesammelten Fenster werden paarweise ausgerichtet
(`local_align_pair`): Die Vorwärtsläufe zweier Fenster rechnen Zeile für Zeile in derselben
Schleife. Die Hauptschleife eines Laufs ist eine Kette, in der jeder Schritt auf den vorigen
wartet; zwei unabhängige Ketten nebeneinander lassen die Rechenwerke weniger leer laufen. Jeder
Lauf rechnet genau seine Schritte wie allein, ein Test vergleicht beide Wege auf zufälligen
Reads und Fenstern (gleich und verschieden lang, beide Breiten, mit und ohne Anfangssuche).
Auf 421.564 HG002-Paaren war die BAM byte-gleich. Schalter zurück:
`PairOptions::rescue_two = false`; Messprogramm `examples/swbench` mit `PAIR=1`. Auf dem Mac
(M5) ist der Rechenkern allein so 32 % schneller, im ganzen Lauf auf 1 Mio. PJ-Paaren aber nur
rund 1 %: Rechnen alle Kerne gleichzeitig, schrumpft der Gewinn auf etwa 17 %, und unter
Dauerlast begrenzt die Wärme des lüfterlosen Geräts den Takt.

**Weniger Bereiche.** Die Partnersuche läuft für bis zu `PairOptions::max_matesw` Bereiche
eines Reads (bwa `-m`, 50), die fast so gut sind wie der beste. Bei Reads aus Wiederholungen
sind das viele Fenster, von denen rund 70 % nichts finden. `fq2bam --partnersuche-bereiche`
setzt die Zahl; weniger spart Rechenzeit, ändert aber einzelne Zeilen. Auf 421.564
HG002-Paaren (GRCh38, Schnellweg an):

| Bereiche | Fenster | Zeilen anders (von 844.609) |
|---|---|---|
| 50 (bwa) | 218.556 | 0 |
| 10 | 83.343 | 1.070 |
| 3 | 32.974 | 1.415 |
| 1 | 12.232 | 1.882 |

Bei 10 sind die Unterschiede zu 72 % nur ein anderer MAPQ (meist etwas höher, weil weniger
gerettete Nebenstellen den zweitbesten Score drücken), 21 % andere Zusatzfelder, 5 % ein anderer
Ort, davon fast alle bei Reads mit MAPQ 0.

## Was in den SAM-Zeilen dazukommt

| Feld | Bedeutung |
|---|---|
| Flag 1 | Read gehört zu einem Paar |
| Flag 2 | richtiges Paar (Richtung und Abstand passen) |
| Flag 8 | Partner ohne Ort |
| Flag 32 | Partner auf dem Gegenstrang |
| Flag 64, 128 | Read 1, Read 2 |
| Spalte 7, 8 | Contig (`=` wenn derselbe) und Position des Partners |
| Spalte 9 (TLEN) | Länge des Stücks von äußerem Ende zu äußerem Ende; positiv beim linken Read, negativ beim rechten |
| `MC` | CIGAR des Partners |
| `MQ` | MAPQ des Partners |

Ein Read ohne Ort, dessen Partner einen hat, steht mit Flag 4 an Contig und Position des
Partners, ohne CIGAR. So bleiben die beiden beim Sortieren nebeneinander.

## Blöcke: wovon das Ergebnis abhängt
Die Insertgröße wird je Block geschätzt. Ein Block endet, sobald 10 Millionen Basen gelesen
sind (`--block-basen`). Zwei Folgen:

- Das Ergebnis eines Paares hängt von den anderen Paaren seines Blocks ab. Bei bwa hängt die
  Blockgröße außerdem an der Zahl der Kerne (10 Millionen Basen je Kern), wenn man `-K` nicht
  setzt; dieselben Daten geben dort mit 8 und mit 16 Kernen leicht verschiedene Ergebnisse.
  Genome 4.0 bildet die Blöcke unabhängig von der Zahl der Kerne.
- Die laufende Nummer des Paares zählt über die Blöcke weiter. Sie entscheidet zwischen gleich
  guten Orten.

## Feinheiten, die aus bwa übernommen sind
- **Lokales Alignment in Spuren.** bwa rechnet `ksw_align2` mit Vektorbefehlen in gestreifter
  Anordnung (Farrar 2007). Diese Anordnung weicht in seltenen Fällen von der Lehrbuch-Rechnung
  ab (eine Deletion direkt nach einer Insertion gibt es nicht). `local.rs` rechnet dieselben
  Schritte in derselben Anordnung und liefert in 1.000 Prüffällen jeden Wert wie bwa (Score,
  Enden, zweitbester Score, Anfänge). Die Spuren sind feste Felder aus 16 Bytes oder 8 Werten
  mit 16 Bit, die Schritte gesättigte Additionen, Subtraktionen und Maxima. Die Vektoren kommen
  aus `gc-sys` (`simd.rs`): SSE2 auf Intel und AMD, NEON auf ARM, sonst gewöhnliche Felder. Ein
  Test prüft jede Operation gegen die Felder. Vorher sollte der Übersetzer die Felder selbst in
  Vektorbefehle verwandeln; das gelang nur zum Teil. Mit den Vektoren ist das lokale Alignment
  2,6-mal so schnell (148 gegen 1.000 Basen) und das Ausrichten auf chr20 10 bis 15 % schneller,
  bei gleicher BAM. Schalter `GENOME_SIMD=0` rechnet wieder mit den Feldern. Das Ergebnis ist
  auf jedem System dasselbe. Die erste Fassung, die jeden Wert einzeln mit 32 Bit
  rechnet, bleibt als Prüfmaßstab im Test: 3.000 Zufallsfälle, beide Breiten, jeder Wert gleich.
- **Bei der Partnersuche gefundene Bereiche** haben keinen Anteil häufiger Seeds (er bleibt 0)
  und werden für den CIGAR mit schmalem Band ausgerichtet, wie in bwa.
- **Paarnummer mit 32 Bit.** bwa führt die Nummer des Paares bei der Wahl zwischen gleich guten
  Paaren als 32-Bit-Zahl. Das ist nachgebildet, damit die Wahl auch bei sehr vielen Paaren
  dieselbe ist.
- **Fehlerfunktion `erfc`.** Sie kommt aus der Rust-Bibliothek `libm` (auf jedem System
  dieselbe Rechnung); bwa nimmt die des Betriebssystems. Abweichungen in der letzten Stelle
  wären möglich, im Vergleich ist keine aufgetreten.
- **`MC` in ergänzenden Zeilen** zeigt die abgeschnittenen Enden des Partners als harten Schnitt,
  wie in bwa.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Lokales Alignment: 1.000 Fälle mit den Ergebnissen von `ksw_align2` aus bwa, alle sieben Werte | alle gleich (`tests/lokal_wie_bwa.rs`) |
| Lokales Alignment: 8-Bit- und 16-Bit-Rechnung gleich, Score nie über der Lehrbuch-Rechnung | bestanden |
| Lokales Alignment: schnelle Fassung gegen die Rechnung Wert für Wert, 3.000 Zufallsfälle mit N, Abbruch und beiden Breiten | alle gleich |
| Richtung und Abstand zweier Orte | bestanden |
| Insertgröße aus 40 Paaren mit Stücken von 380 bis 420; keine Schätzung unter 10 Paaren | bestanden |
| Gewöhnliches Paar: Flags 99 und 147, Partnerfelder, TLEN ±Stücklänge, `MC`, `MQ` | bestanden |
| Partner ohne Seed (jede 12. Base falsch) wird neben dem Partner gefunden | bestanden |
| Partner ohne Ort: Flags 73 und 133, steht beim Ort des Partners | bestanden |
| Read mit zwei gleich guten Orten: der Partner entscheidet, MAPQ über 20, anderer Ort in `XA` | bestanden |
| Ganze Zeilen gegen bwa, chr20, 100.000 Paare | siehe [Messung](../../bench/M3-paare.md) |
| Richtigkeit gegen den wahren Ort, mit und ohne Partner | siehe [Messung](../../bench/M3-paare.md) |

## Grenzen
- Verglichen mit dem in pybwa eingebauten bwa (0.7.18). Dessen Ausgabe trägt die Flag-Bits 2, 64
  und 128 nicht; diese Bits sind nur durch eigene Tests und den Quelltext von bwa belegt. Der
  Vergleich mit dem Programm bwa 0.7.19 steht aus.
- **N-Strecken.** An N-Stellen der Referenz rechnet Genome 4.0 mit "unbekannt"
  ([Ausgabe als SAM](sam-ausgabe.md)), bwa mit seinen Zufallsbasen. Fenster der Partnersuche
  und Verlängerungen, die an eine N-Strecke stoßen, können sich deshalb um einzelne Basen und
  Punkte unterscheiden. Betrifft nur Orte direkt an N-Strecken.
- **Sehr häufige Seeds.** Kommt ein Seed mehr als 500-mal vor, wird nur jede n-te Fundstelle
  benutzt. Welche das sind, kann von bwa abweichen; dann unterscheiden sich in langen
  Wiederholungen (Zentromer) einzelne Bereiche und damit die Wahl zwischen gleich guten Orten.
  Die Ursache ist eingegrenzt, aber nicht bewiesen, siehe Messung.

## Literatur
Li 2013; Farrar 2007; SAM-Spezifikation, siehe [Literatur](../literatur.md).
