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
   **Versuch (Vorgabe aus, `GC_STROBE_MEHRERE=1`):** Ist der Zweitbeste zu nah, gibt der
   Strobemer-Weg alle ausgerichteten Orte als Bereiche ab (bester zuerst, einander nicht
   überlappend), statt den vollen Weg zu gehen. Primären Ort, Zweitbesten, Zahl fast gleich guter
   Orte und die Paarbildung bestimmen dann dieselben Schritte wie nach dem vollen Weg. Das gilt
   nur, wenn kein Randstrobe gesperrt war und hinter den 16 ausgerichteten Gruppen keine mehr kam,
   die mindestens ein Drittel so groß ist wie die größte; sonst könnten weitere Kopien fehlen und
   die MAPQ zu hoch werden. Simulation (600.000 Reads, ganzes GRCh38): Ausrichten 22,0 s statt
   23,1 s, falsch mit MAPQ ≥ 20 gleich (10), mit MAPQ 1 bis 19 153 statt 151, mit MAPQ 0 85
   richtige Orte mehr. Ohne die Bedingung zu gesperrten Randstrobes stiegen zwei Paare in
   Wiederholungen nahe am Zentromer von MAPQ 0 bis 5 auf 24 bis 59 (falsch). Prüfung gegen GIAB
   steht aus.
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
