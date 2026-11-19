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
  mit 16 Bit, die Schritte gesättigte Additionen, Subtraktionen und Maxima. Der Übersetzer macht
  daraus Vektorbefehle (SSE2 auf Intel und AMD, NEON auf ARM), ohne systemnahen Code; das
  Ergebnis ist auf jedem System dasselbe. Die erste Fassung, die jeden Wert einzeln mit 32 Bit
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
