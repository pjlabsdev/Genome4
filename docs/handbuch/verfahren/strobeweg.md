# Strobemer-Weg

Stand: Versuch (`fq2bam --strobeweg an`, Vorgabe aus). Die Ausgabe ist nicht bwa-gleich und wird
gegen die Wahrheit geprüft, nicht gegen bwa.

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

Der Index (`gc_index::StrobeIndex`) hält je Randstrobe 8 Bytes: 32 Bit des Hashes und die
Position. Die obersten 24 Bit wählen das Fach. Gebaut wird beim Start aus der Referenz des
FM-Index, in Stücken von 4 Mb mit 10 kb Überhang auf allen Kernen; das Ergebnis hängt nicht von
der Zahl der Kerne ab (Test). N-Strecken liefern keine Syncmere. GRCh38: 586,6 Mio. Einträge,
4,7 GB, rund 60 s auf 4 Kernen.

## Ablauf je Read
1. Randstrobes des Reads und seines reversen Komplements nachschlagen. Ein Randstrobe mit mehr als
   1000 Orten gilt als gesperrt. Sind mehr als 25 % gesperrt, geht der Read den vollen Weg
   (`GC_STROBE_GESPERRT` setzt den Anteil): Sonst könnte ein gleich guter Ort in der Wiederholung
   verborgen sein. Bei höchstens 25 % findet der übrige Read jede nahe Kopie mit.
2. Treffer nach Diagonale gruppieren (±10), die 16 größten Gruppen nehmen.
3. Je Gruppe das Fenster mit 12 Basen Rand lokal ausrichten (Smith-Waterman mit Lücken). Eine
   Gruppe, deren beweisbare 5-mer-Schranke ([Schnellweg](schnellweg.md#kopien-ohne-ausrichtung-ausschließen))
   höchstens den besten Score minus 15 erreicht, wird nicht ausgerichtet; ihre Schranke zählt als
   Zweitbester (vorsichtig: die MAPQ wird dadurch höchstens kleiner).
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

Zähler im Protokoll (`ALIGN.FERTIG`, Feld `strobeweg`): `versucht`, `genommen`, `gesperrt`,
`ohne_treffer`, `abgeschnitten`, `zu_schwach`, `zweiter_zu_nah`.

## Offen
- Vergleich gegen GIAB HG002 (Varianten), getrennt nach schwierigen Regionen.
- Index einmal auf die Platte schreiben statt bei jedem Start bauen.
- Speicher: zusammen mit dem FM-Index am ganzen Genom rund 15 GB; für 16 GB zu knapp.

## Literatur
Sahlin K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment.
Genome Biology 23, 260 (2022).
