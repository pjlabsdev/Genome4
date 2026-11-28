# Strobemer-Weg

Stand: Vorgabe seit 2026-11-28 (`fq2bam --strobeweg auto`: an ab 14 GiB Arbeitsspeicher; pjs
Entscheidung nach dem GIAB-Vergleich). Die Ausgabe ist nicht bwa-gleich und wird gegen die Wahrheit
geprüft, nicht gegen bwa; `--strobeweg aus` gibt den bwa-gleichen Weg.

## Aufgabe
Der volle Weg (Seed-Suche mit Nachsuche, Ketten, Verlängern) kostet bei Reads, die der
[Schnellweg](schnellweg.md) nicht nimmt, den größten Teil der Rechenzeit. Viele dieser Reads
passen trotzdem eindeutig an eine Stelle, nur mit zwei bis fünf Abweichungen, einer kleinen Lücke
oder einer Abweichung am Ende. Der Strobemer-Weg findet diese Stelle direkt und gibt den Read nur
dann weiter, wenn die Stelle eindeutig ist.

## Randstrobes
Wie bei strobealign (Sahlin 2022). Aus dem Vorwärtsstrang werden offene Syncmere gewählt
(k = 20, s = 16: das kleinste der fünf s-mere eines k-mers steht an Stelle 2, rund jedes fünfte
k-mer). Jedes Syncmer wird mit einem der 2. bis 9. folgenden verbunden, dem mit der kleinsten
Bitdifferenz der Hashes. Ein Randstrobe überspannt 25 bis etwa 70 Basen und bleibt erhalten, wenn
eine Abweichung zwischen den beiden Strobes liegt.

Der Index (`gc_index::StrobeIndex`) hält je Ort 6 Bytes: 16 Bit des Hashes und die Position.
Die obersten 24 Bit des Hashes wählen das Fach, zusammen unterscheiden 40 Bit. Randstrobes mit
mehr als 1000 Orten werden nicht gespeichert, nur ihr Hash in einer Sperrliste (Test mit einer
Wiederholung aus 1200 Kopien). Gebaut wird beim Start aus der Referenz des FM-Index, in Stücken
von 4 Mb mit 10 kb Überhang auf allen Kernen; die Fächer werden in vier Durchgängen gefüllt,
sortiert und verdichtet, so braucht der Bau neben dem fertigen Index nur ein Viertel der Einträge
zusätzlich. Das Ergebnis hängt nicht von der Zahl der Kerne ab (Test). N-Strecken liefern keine
Syncmere. GRCh38: 576,7 Mio. Orte, 5.049 gesperrte Randstrobes, 3,5 GB (vorher mit 8 Bytes je
Ort und allen Orten 4,8 GB), rund 90 s auf 4 Kernen.

Mit `--strobe-datei PFAD` wird der Index einmal gebaut und gespeichert (Kennung `GCSTRB01`,
Fingerabdruck der Referenz und der Einstellungen, dann die vier Felder in Little-Endian), danach
nur geladen. Gehört die Datei zu einer anderen Referenz, wird neu gebaut und ersetzt; eine
gekürzte oder beschädigte Datei bricht ab (Test). Geschrieben wird erst in `PFAD.teil`, dann
umbenannt, damit ein abgebrochener Lauf keine halbe Datei hinterlässt. GRCh38 in der Cloud auf 4
Kernen: bauen 117 s, laden 6,1 s; die BAM ist in beiden Fällen dieselbe.

## Ablauf je Read
1. Randstrobes des Reads und seines reversen Komplements nachschlagen. Ein Randstrobe mit mehr als
   1000 Orten gilt als gesperrt. Sind mehr als 25 % gesperrt, geht der Read den vollen Weg
   (`GC_STROBE_GESPERRT` setzt den Anteil): Sonst könnte ein gleich guter Ort in der Wiederholung
   verborgen sein. Bei höchstens 25 % findet der übrige Read jede nahe Kopie mit.
   Die Randstrobes beider Stränge werden zuerst alle berechnet, dann ihre Fächer in zwei Stufen
   vorab in den Cache geholt (erst der Fachanfang, dann Schlüssel und Positionen). Die Abfragen
   hängen nicht voneinander ab, der Prozessor wartet so auf viele Speicherzeilen gleichzeitig.
   Dieselben Treffer, auf 400.000 leichten Simulationspaaren 5,7 % weniger Ausrichtezeit
   (Cloud, zwei Läufe je Fassung, Sätze gleich); am Mac (10 Mio. HG002-Paare, je zwei Läufe)
   mit 397 und 465 s gegen ohne 441 und 467 s, im Mittel −5 %, aber die Läufe streuen dort
   stärker als der Unterschied (Auslagerung wächst von Lauf zu Lauf). Ebenso werden die Referenzfenster der Kandidaten vorab geholt,
   bevor der erste ausgerichtet wird (Cloud weitere rund 5 %). Schalter zurück für beides:
   `GC_STROBE_VORAB=0`.
2. Treffer nach Diagonale gruppieren (±10), die 16 größten Gruppen nehmen.
3. Je Gruppe zuerst ohne Lücke entlang der Diagonale vergleichen (bester Abschnitt, Kadane).
   Liegt er höchstens 30 unter dem vollen Score, lassen sich beide Enden außerhalb ohne Lücke für
   höchstens den Abschneide-Abzug anhängen und erreicht die 5-mer-Schranke nicht mehr als dieser
   Score (oder er ist ohnehin der neue Beste), entfällt die lokale Ausrichtung. Schalter zurück:
   `GC_STROBE_HAMMING=0`. Sonst das Fenster mit 12 Basen Rand lokal ausrichten (Smith-Waterman
   mit Lücken). Eine
   Gruppe, deren beweisbare 5-mer-Schranke ([Schnellweg](schnellweg.md#kopien-ohne-ausrichtung-ausschließen))
   höchstens den besten Score minus 15 erreicht, wird nicht ausgerichtet; ihre Schranke zählt als
   Zweitbester (vorsichtig: die MAPQ wird dadurch höchstens kleiner).
   Vor jeder Annahme ohne Lücke (ganzer Read beim Vergleich entlang der Diagonale, Enden beim
   Verlängern in Schritt 4) wird geprüft, ob eine einzelne Lücke von 1 bis 3 Basen besser passt
   (Vor- und Nachsummen auf den Nachbardiagonalen, Kosten wie beim vollen Weg). Dann richtet die
   lokale Ausrichtung den Read aus, oder beim Verlängern geht er den vollen Weg. Ohne diese
   Prüfung wurden Indels nahe am Read-Ende (`144M1D6M`) als `150M` mit Abweichungen geschrieben:
   GIAB HG002 chr20 Indel F1 0,7347 statt 0,7357. Simulation: Reads, deren Indel ohne Strobemer-Weg
   erkannt wird und mit ihm nicht, 62 vor, 2 nach der Prüfung; Ausrichten rund 8 % langsamer.
   Schalter zurück: `GC_STROBE_LUECKENPRUEFUNG=0`.
   Seit 1 bis 15 Basen geprüft werden, gehen auch 4 bis 7 Basen lange Indels in Wiederholungen
   (GIAB: `137M4D14M` als `151M`, 167 Reads auf chr20 mit 1 bis 3 Basen) nicht mehr verloren.
   Geprüft werden nur Längen `d`, mit denen eine Lücke überhaupt mehr Punkte als ohne Lücke
   erreichen kann (Deletion höchstens alle Basen minus `o_del + d·e_del`, Insertion `d` Basen
   weniger minus `o_ins + d·e_ins`); ein Read mit höchstens einer Abweichung braucht gar keine
   Prüfung. Die Antwort bleibt dieselbe, die Arbeit schrumpft. Reicht das Fenster (Rand 12) für die
   nötige Länge nicht, richtet die lokale Ausrichtung aus. Simulation (600.000 Reads): Hauptdatensätze
   gleich wie mit 3 Basen, Zeit im Rauschen. Schalter zurück: `GC_STROBE_LUECKE_BIS=3`.
   Ist der Zweitbeste schon weniger als 15 unter dem Besten, bricht der Read ab und geht den
   vollen Weg, ohne die übrigen Gruppen auszurichten (Schalter zurück: `GC_STROBE_FRUEH=0`).
4. Abgeschnittene Enden ohne Lücke bis zum Read-Ende verlängern, wenn das je Seite höchstens den
   Abzug fürs Abschneiden kostet (5, wie bwa `-L 5,5`). Sonst geht der Read den vollen Weg.
5. Genommen wird der Read, wenn der Score höchstens 30 unter dem vollen liegt und der
   Zweitbeste mindestens 15 schlechter ist (`GC_STROBE_ABSTAND`). Er bekommt genau einen Bereich
   mit dem Zweitbesten als `sub`; Partnersuche, Paarbewertung, MAPQ und Zeilen laufen danach
   unverändert.

## Messung
Simulation über das ganze GRCh38 (chr1 bis 22, X, Y), 300.000 Paare, 0,4 % Abweichungen und
0,04 % Lücken je Base, Wahrheit im Namen; 4 Kerne, `--schnellweg an`:

| | Ausrichten | falsch bei MAPQ ≥ 20 | MAPQ 60 (falsch) |
|---|---|---|---|
| `--strobeweg aus` | 30,6 s | 10 | 566.678 (0) |
| `--strobeweg an` | 23,1 s | 10 | 567.233 (0) |

Der Weg nimmt 249.000 der 292.000 Reads, die der Schnellweg nicht nimmt (85 %). Zu nahe Zweitorte
zuzulassen (Abstand 5 statt 15) spart weitere 2 s, erhöht aber die falschen Reads bei MAPQ ≥ 20
auf 37; verworfen. Gesperrte Randstrobes bis 50 % zuzulassen ergab 16 falsche; verworfen.

Vergleich ohne Lücke zuerst (gleiche Simulation, 4 Kerne in der Cloud): Ausrichten 22,4 s statt
25,3 s, lokale Ausrichtungen 247.991 statt 570.132. Zähler und Treffergenauigkeit gleich (falsch
bei MAPQ ≥ 20: 10, MAPQ 60: 567.233 ohne Fehler), 32 von 600.000 Sätzen anders (gleich gute
Stellen, Anfang aus der Diagonale). Früher Abbruch bei zu nahem Zweiten dazu: Ausrichten 20,4 s, lokale Ausrichtungen 102.272,
genommen 248.450 statt 248.972, falsch bei MAPQ ≥ 20 weiter 10, MAPQ 60 ohne Fehler.
Auf echten HG002-Reads kostete die lokale Ausrichtung samt
Rückwärtslauf drei Viertel der Zeit des Strobemer-Wegs (Profil am Mac).

Zähler im Protokoll (`ALIGN.FERTIG`, Feld `strobeweg`): `versucht`, `genommen`, `gesperrt`,
`ohne_treffer`, `abgeschnitten`, `zu_schwach`, `zweiter_zu_nah`, `ohne_luecke` (Kandidaten ohne
lokale Ausrichtung), `lokal` (Kandidaten mit lokaler Ausrichtung).

## Offen
- Vergleich gegen GIAB HG002 (Varianten), getrennt nach schwierigen Regionen.
- Speicher: am Mac (16 GB) passte die erste Fassung (4,8 GB) neben FM-Index und Sortierpuffer nicht ohne Auslagern; mit 3,5 GB und `--sa-ausduennen 8` neu zu prüfen.

## Literatur
Sahlin K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment.
Genome Biology 23, 260 (2022).

## Versuch: Wiederholungs-Reads ohne vollen Weg (2026-11-28, Vorgabe aus)
Reads mit zu nahem Zweitem gehen heute den vollen Weg. Der Versuch `GC_STROBE_NAH_NEHMEN=2`
gibt stattdessen alle fast gleich guten Strobe-Orte (bis `pen_unpaired` unter dem besten) an
die Paarung. `GC_STROBE_NAH_VORSICHT` legt fest, wann doch der volle Weg geht (`gesperrt`,
`gruppen`, `0`, sonst beides). `GC_STROBE_NAH_GESPERRT_NULL=1` nimmt Reads mit gesperrten
Randstrobes trotzdem, aber mit MAPQ 0, und ein Partner auf so einem Ort stützt das Paar nicht.
`GC_STROBE_NAH_NEHMEN=1` (ein Ort, `GC_STROBE_NAH_CSUB=voll`) ist die erste, verworfene Form.
Simulation über GRCh38, 2 Mio. Paare, 4 Kerne:

| Variante | Ausrichten | falsch mit MAPQ ≥ 20 | MAPQ 0 |
|---|---|---|---|
| Vorgabe (voller Weg) | 200 s | 74 | 94.708 |
| `NAH_NEHMEN=1` | 114 s | 10.928 | 99.054 |
| `NAH_NEHMEN=2`, `VORSICHT=0` | 140 s | 557 | 93.534 |
| `NAH_NEHMEN=2` (Vorsicht beides) | 205 s | 79 | 94.596 |
| `NAH_NEHMEN=2`, `VORSICHT=gesperrt` | 182 s | 149 | 94.210 |
| dazu `GESPERRT_NULL=1` | 130 s | 189 | 158.377 |
| `VORSICHT=gruppen`, `GESPERRT_NULL=1` | 150 s | 139 | 158.602 |

Die Zeit steckt fast ganz in Reads aus sehr häufigen Wiederholungen (gesperrte Randstrobes);
dort fehlen dem Strobeweg Kopien. Ob die MAPQ-0-Variante Varianten kostet, klärt GIAB.

## Versuch: Strobe-Saat (2026-11-28, Vorgabe aus)
`GC_STROBE_SAAT=1`: Reads mit zu nahem Zweitem und ohne gesperrte Randstrobes bekommen ihre
Seeds aus den Strobe-Treffern (je Diagonale die Strecken gleicher Basen ab `min_seed_len`)
statt aus der SMEM-Suche; Ketten und Verlängern laufen wie beim vollen Weg. Simulation wie
oben: 199 s statt 245 s Ausrichten (der Vergleichslauf war gestört), 93 statt 74 falsch mit
MAPQ ≥ 20. Nicht weiter verfolgt.

## Versuch: Strobe-Paarweg (2026-11-28, Vorgabe aus)
`GC_STROBE_PAARWEG=1` (Quelltext `gc-align/src/paarweg.rs`). Erst Schnell- und Strobeweg je
Read wie bisher. Bleibt bei einem Paar ein Read ohne Bereich, werden beide Reads gemeinsam aus
ihren Strobe-Ankern verortet:

1. Treffer je Read nach Diagonale bündeln (Anker, mit Zahl der Treffer).
2. Ankerpaare: ein Anker vorwärts, der des Partners auf dem Gegenstrang dahinter, Abstand bis
   `GC_STROBE_PAAR_MAX` (1500). Geordnet nach der Summe der Treffer.
3. Alle Ankerpaare ab der Hälfte der besten Trefferzahl ausrichten, dazu je Read die drei
   besten Einzelanker (`GC_STROBE_PAAR_EINZELN`). Sind es mehr als `GC_STROBE_PAAR_TOP` (16)
   Ankerpaare, geht das Paar den bisherigen Weg.
4. Ausrichten ohne Lücke bei höchstens `GC_STROBE_PAAR_FEHLER` (2) Fehlpaarungen, sonst lokal;
   abgeschnittene Enden ohne Lücke anhängen wie beim Strobeweg. Fast deckungsgleiche Orte
   (95 % des kürzeren, wie bwa) zählen einmal; Orte, die sich nur verschoben überlappen
   (kurze Tandemwiederholungen), bleiben getrennt.
5. Jeder fast beste Ort (bis `pen_unpaired` unter dem besten) eines Reads braucht einen
   Partnerort in Paarlage. Fehlt er, könnte dort ein gleich gutes Paar liegen, das nur die
   Partnersuche fände: bisheriger Weg. Ausnahme: Ist der Partner für sich eindeutig
   (Strobeweg genommen), entfällt die Prüfung für den anderen Read
   (`GC_STROBE_PAAR_EINDEUTIG=0` schaltet das ab).
6. Wiederholungsanteil wie bwa `frac_rep`: Anteil gesperrter Randstrobes
   (`GC_STROBE_PAAR_WDH=0` setzt ihn auf 0).
7. Genommene Paare laufen durch die übliche Paarbewertung, aber ohne Partnersuche. Ein schon
   verorteter Read behält seinen Bereich und bekommt die Orte seiner Anker dazu
   (`GC_STROBE_PAAR_DAZU=0`).

Weitere Schalter für Versuche: `GC_STROBE_PAAR_UNSICHER=N` nimmt Paare aus Schritt 3 und 5
selbst, mit MAPQ höchstens N; dazu `GC_STROBE_PAAR_UNSICHER_VOLL=1`: voller Weg ohne
Partnersuche. `GC_STROBE_PAAR_SUCHE=1`: Partnersuche statt vollem Weg (nur ohne gesperrte
Randstrobes, `GC_STROBE_PAAR_SUCHE_SAUBER=0` immer).

Simulation über GRCh38, 2 Mio. Paare, 4 Kerne, Vergleich 214 s Ausrichten, 74 falsch mit
MAPQ ≥ 20, 3.891.169 richtig mit MAPQ ≥ 20:

| Fassung | Ausrichten | falsch, MAPQ ≥ 20 | richtig, MAPQ ≥ 20 |
|---|---|---|---|
| erste Fassung (ohne Schritt 3 bis 6) | 75 s | 7.673 | |
| Schritt 3 und 4 | 153 s | 2.694 | |
| dazu Schritt 5 (ohne Ausnahme) | 195 s | 286 | |
| dazu Schritt 6 | 207 s | 202 | 3.892.221 |
| dazu Ausnahme in Schritt 5 (Vorgabe des Versuchs) | 195 s | 222 | 3.892.218 |
| `UNSICHER=0` | 80 s | 179 | 3.821.565 |
| `UNSICHER=0`, `UNSICHER_VOLL=1` | 177 s | 987 | |
| `SUCHE=1` (ohne Schritt 6) | 155 s | 2.076 | |

Ergebnis: Die Zeit steckt in den etwa 3 % der Paare, bei denen beide Reads in Wiederholungen
liegen. Sie kosten im vollen Weg rund 0,8 ms je Read (ein Strobe-Read etwa 0,02 ms) und
brauchen den vollen Weg samt Partnersuche, sonst steigen die Fehler oder es gehen richtige
sichere Orte verloren.

## Zeitmessung je Paarart (`GC_ZEITEN=1`)

Für die Arbeit am Tempo misst `GC_ZEITEN=1` die Rechenzeit je Teil und je Paarart und schreibt sie
unter `zeiten` in `ALIGN.FERTIG`: schneller Versuch je Read, voller Weg je Read, Partnersuche und
Zeilen je Paarart (0 = beide Reads über den Strobeweg, 1 = einer über den vollen Weg, 2 = beide).
Ohne Schalter misst sie nichts; die Sätze bleiben in jedem Fall gleich. Cloud, 2 Mio.
Simulationspaare (Rechenzeit aller Kerne):

| Teil | Reads oder Paare | Zeit | je Stück |
|---|---|---|---|
| Strobeweg-Versuch | 4.000.000 Reads | 115 s | 29 µs |
| Voller Weg | 278.372 Reads (7 %) | 415 s | 1,49 ms |
| Partnersuche, Art 2 | 108.743 Paare | 40 s | 371 µs |
| Partnersuche, Art 1 | 60.886 Paare | 7,8 s | 127 µs |
| Zeilen, Art 0 | 1.830.371 Paare | 5,3 s | 2,9 µs |
| Zeilen, Art 2 | 108.743 Paare | 16,4 s | 151 µs |

Der volle Weg für 7 % der Reads kostet zwei Drittel der Zeit.

## Versuch Strobe-Vollweg: mehr Kandidaten, billig geprüft (2026-11-28, Vorgabe aus)

Im Versuch `GC_STROBE_NAH_NEHMEN=2` fehlte bei mehr als 16 gleich guten Kopien (Segment-
Duplikationen, etwa chr1 q21) oft die richtige Kopie in der Kandidatenliste: Die Gruppen sind
nach Trefferzahl und dann nach Ort geordnet, bei Gleichstand gewinnen die vorderen Orte. Dann
legt die Paarung das ganze Paar auf eine falsche Kopie, und weil das gleich gute zweite Paar
fehlt, ist die MAPQ zu hoch. `GC_STROBE_KANDIDATEN=N` setzt die Zahl der Kandidaten (Vorgabe 16).
`GC_STROBE_BILLIG_AB=N` prüft Kandidaten ab dem N-ten nur ohne Lücke (ohne Lückenprüfung und
ohne lokale Ausrichtung); passt einer so nicht, zählt nur seine 5-mer-Schranke als Zweiter.
Auf echten Daten (Mac, 10 Mio. HG002-Paare) war ohne diese Begrenzung die lokale Ausrichtung
dreimal so häufig (15,0 statt 4,7 Mio.), und das Ausrichten wurde langsamer (432 s statt 380 s)
bei gleicher Qualität (MAPQ ≥ 20: 93,96 % statt 93,93 %, richtig gepaart 98,61 %).

Simulation über GRCh38, 2 Mio. Paare, 4 Kerne, alle mit `NAH_NEHMEN=2`, `NAH_VORSICHT=0`:

| Fassung | Ausrichten | falsch, MAPQ ≥ 20 | Reads im vollen Weg |
|---|---|---|---|
| Vorgabe (ohne Versuch) | 155 s | 74 | 278.372 |
| 16 Kandidaten | 88 s | 557 | 58.458 |
| 64 Kandidaten | 112 s | 247 | |
| 256 Kandidaten | 170 s | 101 | 57.783 |
| 256 Kandidaten, billig ab 1 | 102 s | 105 | 57.785 |
| 64 Kandidaten, billig ab 4 | 92 s | 252 | 57.977 |
