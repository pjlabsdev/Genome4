# Trimmen von Read-Paaren (Ersatz für fastp)

Stand: gebaut (`gc-io`, Modul `trim`), in `genome fq2bam --trimmen` eingebaut und Read für Read gegen fastp 1.4.0 geprüft: an je 1 Mio. Paaren NovaSeq und MGI keine Abweichung ([Ergebnis](../../../bench/m9/ERGEBNIS.md)). Die Seite beschreibt, was fastp 1.4.0 mit den
Einstellungen der alten Kette tut. `genome fq2bam --trimmen` soll genau das tun (Absicht vom
19.11.2026, pj: Trimmen ist eine Option, ohne sie bleiben die Reads unverändert).

Quelle: fastp 1.4.0, Commit `ee125a52fddae9d8189e13cd8a33ad04988a1351`, Aufruf der alten Kette
`fastp -i R1 -I R2 -o … -O … --detect_adapter_for_pe --thread 8`, sonst Vorgaben. Jede Regel
unten nennt die Datei in `src/`, aus der sie stammt.

## Vorgaben, die gelten

| Einstellung | Wert | Quelle |
|---|---|---|
| Adapter entfernen | an | `main.cpp` (`disable_adapter_trimming` fehlt) |
| Adapter für Paare suchen | an | `--detect_adapter_for_pe` |
| Überlappung: Mindestlänge | 30 | `overlap_len_require` |
| Überlappung: höchstens Unterschiede | 5 | `overlap_diff_limit` |
| Überlappung: höchstens Anteil Unterschiede | 20 % | `options.cpp` |
| Überlappung mit Lücke | aus | `allow_gap_overlap_trimming` |
| Adapter-Dimer bis Länge | 2 | `dimer_max_len` |
| PolyG kürzen | automatisch, siehe unten; Mindestlänge 10 | `main.cpp`, `evaluator.cpp` |
| PolyX, Kürzen nach Qualität (cut_front, cut_tail, cut_right) | aus | Vorgabe |
| feste Kürzung vorn oder hinten, Höchstlänge | aus (0) | Vorgabe |
| Base gilt als gut ab Qualität | 15 | `qualified_quality_phred` |
| höchstens Anteil schlechter Basen | 40 % | `unqualified_percent_limit` |
| höchstens N | 5 | `n_base_limit` |
| mittlere Qualität | keine Forderung | `average_qual` 0 |
| Mindestlänge | 15 | `length_required` |
| Komplexitätsfilter, Duplikate entfernen, Korrektur, Zusammenführen | aus | Vorgabe |

## Vor dem Lauf: zwei Entscheidungen aus dem Anfang der Dateien

**PolyG an oder aus** (`Evaluator::isTwoColorSystem`): Nur der Name des allerersten Reads von
R1 zählt. Beginnt er mit `@FS`, `@MN`, `@SH`, `@NS`, `@NB`, `@NDX`, `@VL`, `@VH`, `@A`, `@NA` oder
`@LH` (Geräte mit zweifarbiger Chemie), wird PolyG gekürzt, sonst nicht. MGI-Namen (`@V3…`,
`@E…`, `@DP8…`) fallen nicht darunter. Wurden die Reads umbenannt (zum Beispiel beim Herunterladen
aus dem SRA), entscheidet der neue Name.

**Adapterfolge je Read-Seite** (`Evaluator::evalAdapterAndReadNum`), getrennt für R1 und R2:

1. Die ersten Reads lesen, bis 262.144 Reads oder 151 × 262.144 Basen erreicht sind. Weniger als
   10.000 Reads: kein Adapter.
2. Bekannte Adapter (`knownadapters.h`, Liste mit Namen) prüfen, an höchstens 100.000 Reads:
   Ein Adapter trifft einen Read, wenn er an irgendeiner Stelle `pos` (0 bis Länge minus 9) mit
   höchstens einem Unterschied je 16 verglichene Basen passt; je Read und Adapter zählt nur die
   erste passende Stelle. Abbruch, sobald ein Adapter mehr als 1.000 Treffer hat; Adapter mit weniger als einem
   Zehntel der bisher meisten Treffer werden ab 20 Treffern übersprungen. Gewinner ist der mit den
   meisten Treffern, bei Gleichstand der in alphabetischer Reihenfolge erste. Er gilt, wenn er
   mehr als 1/50 der geprüften Reads trifft, oder mehr als 1/200 und die Summe seiner
   Unterschiede kleiner als die Zahl der geprüften Reads ist.
3. Sonst: Häufigste Teilfolgen der Länge 10 ab Position 20 (letzte Base ausgenommen) zählen,
   arme Folgen (eine Base 6-mal oder öfter, 8 oder mehr G/C, Beginn mit GGGG, AAAAAAAAAA)
   auslassen, die zehn häufigsten der Reihe nach prüfen (mindestens 10 Treffer und 20-fach über
   dem Durchschnitt) und um die häufigste Fortsetzung in beide Richtungen verlängern
   (`getAdapterWithSeed`, `NucleotideTree`). Passt das Ergebnis zu einem bekannten Adapter, gilt
   dieser, sonst die gefundene Folge, wenn die Verlängerung bis zum Ende eindeutig war.
4. Eigenheit: fastp kürzt eine längere Folge mit `resize(0, 60)`, das setzt die Länge auf 0. Ein
   erkannter bekannter Adapter über 60 Basen (139 der 234 bekannten) gilt deshalb als "kein
   Adapter", ohne dass die zweite Suche läuft. Genau so nachgebaut.

## Je Paar, in dieser Reihenfolge (`PairEndProcessor::processPairEnd`)

1. **PolyG**, wenn an (`PolyX::trimPolyG`), für jeden Read einzeln: vom Ende her Basen ansehen;
   jede Nicht-G-Base ist ein Fehler. Abbruch bei mehr als 5 Fehlern, oder ab der 10. Base, sobald
   die Fehler mehr als ein Achtel der angesehenen Basen sind (ganzzahlig). Endet die Suche an
   Stelle 10 oder später (vom Ende gezählt ab 0, die Abbruchbase eingeschlossen) oder erst am
   Read-Anfang, wird der Read ab dem am weitesten vorn gesehenen G abgeschnitten. Wurde gar kein
   G gesehen, fällt nur die letzte Base weg (Anfangswert in fastp).
   Folge davon: Die Suche läuft oft ein paar Basen über das reine PolyG hinaus und schneidet am
   letzten G davor. Beispiel `…ACGT` gefolgt von 15 G: fastp schneidet 17 Basen, nicht 15.
2. **Überlappung** (`OverlapAnalysis::analyze`): R2 wird rückwärts komplementiert und gegen R1
   verschoben. Erst Verschiebungen 0, 1, 2 … (R2 beginnt weiter rechts), solange die Verschiebung
   kleiner als Länge R1 minus 30 ist; dann 0, −1, −2 … (Insert kürzer als der Read), solange größer
   als −(Länge R2 minus 30). Erlaubte Unterschiede: das Kleinere von 5 und 20 % der
   Überlappungslänge (abgerundet). Die erste passende Verschiebung gilt. Wichtig: Bei mehr als 50
   Basen Überlappung entscheiden nur die ersten 50 über den Treffer; Unterschiede dahinter zählen
   nicht. Basen werden als Zeichen verglichen, N gegen N ist gleich.
3. **Adapter über die Überlappung** (`AdapterTrimmer::trimByOverlapAnalysis`): Nur bei negativer
   Verschiebung. Beide Reads werden auf die Überlappungslänge gekürzt (vorher vorn gekürzte Basen
   des Partners werden dazugezählt, hier immer 0).
4. **Adapter über die Folge**, nur wenn Schritt 3 nichts gekürzt hat und für diese Seite eine
   Folge erkannt wurde (`AdapterTrimmer::trimBySequence`, mindestens 4 passende Basen):
   - Ohne Lücke: Stellen ab −4 (Adapter ab 16 Basen; ab 12: −3, ab 8: −2; die erste Base fehlt bei
     Dimeren oft) bis Leselänge minus 5; erlaubt ist ein Unterschied je 8 verglichene Basen. Der
     erste Treffer gilt. Treffer an negativer Stelle: der Read wird leer.
   - Sonst mit einer Einfügung, dann mit einer Auslassung (`Matcher::matchWithOneInsertion`),
     erlaubt ein Unterschied je 8 Basen weniger eins. Eigenheit von fastp 1.4.0: In beiden
     Schleifen wird immer ab dem Anfang des Reads verglichen, nicht ab der Stelle `pos`; nur die
     Vergleichslänge hängt von `pos` ab. Gekürzt wird bei Treffer an `pos`. Das wird genau so
     nachgebaut und im Handbuch als Eigenheit geführt, nicht verbessert, solange das Ziel gleiche
     Ergebnisse sind.
5. **Adapter-Dimer**: Wurde in Schritt 3 oder 4 etwas gekürzt und sind danach beide Reads höchstens
   2 Basen lang, fällt das Paar als Dimer weg.
6. **Filter je Read** (`Filter::passFilter`): leer → Länge; mehr schlechte Basen (Qualität unter
   15) als 40 % der Länge → Qualität; mehr als 5 N → N; kürzer als 15 → Länge. Die Prüfungen
   laufen in dieser Reihenfolge, der erste Grund zählt.
7. **Paar behalten** nur, wenn beide Reads bestehen. Sonst fällt das ganze Paar weg (die alte
   Kette schreibt keine Einzel-Reads).

## Unterschied beim Lesen
`genome` schreibt Basen beim Lesen groß und liest `.` als `N` ([FASTQ-Prüfung](fastq-pruefung.md)),
fastp lässt sie, wie sie sind. Bei üblichen Daten (nur `ACGTN` in Großbuchstaben) ändert das
nichts.

## Umsetzung
`gc-io`, Modul `trim`: `two_color`, `detect_adapter`, `Trimmer::trim_pair`, Zähler `Counts`. Die
Liste bekannter Adapter ist aus `knownadapters.h` erzeugt (MIT-Lizenz, Quelle im Kopf der Datei
`known_adapters.rs`). Tests prüfen jede Regel an kleinen Beispielen, darunter die Eigenheiten.

## In `fq2bam`
Mit `--trimmen` (nur für Paare) liest `fq2bam` vor dem Laden des Index den Anfang beider Dateien
für die Adaptersuche, wie fastp. Dann wird jedes Paar beim Lesen getrimmt; verworfene Paare werden
nicht ausgerichtet und erscheinen nicht in der BAM (wie bei fastp, das sie nicht schreibt). Die
Blöcke für die Schätzung der Insertgröße zählen die Basen nach dem Trimmen. Zwischendateien
entstehen keine. Die Adaptersuche läuft, während der Index lädt; Lesen und Trimmen des nächsten
Blocks laufen in einem eigenen Faden neben dem Ausrichten. Beide Seiten werden gleichzeitig untersucht, die bekannten Adapter blockweise auf allen Kernen (ein übersprungener Adapter bleibt übersprungen, daher genügt je Block die Menge der noch aktiven; Test `known_adapters_parallel_equals_sequential`), die Keime der zweiten Suche gleichzeitig (gewählt wird wie bei fastp der erste der Reihe nach, Test `unknown_adapter_is_found_by_seed`), und eine Vorprüfung der ersten
8 Basen jedes bekannten Adapters spart Vergleiche, ohne ein Ergebnis zu ändern (Test
`differing_bytes_counts_like_plain_loop`; an HG002 chr20 gleiche BAM, Vorbereitung 6,3 s statt
18,3 s bei 2 x 262.144 Reads). Protokoll: `TRIM.ADAPTER`, `TRIM.POLYG`, `TRIM.FERTIG`; die Zähler stehen auch
im Bericht unter `trimmen` (ohne die Option `null`).

## Was `genome` zusätzlich tut
Jede Entscheidung wird gezählt: Paare gesamt, mit PolyG gekürzt (je Seite), über Überlappung
gekürzt, über Folge gekürzt (je Seite), Dimer, verworfen nach Grund, behalten. Die Zahlen gehen
im Protokoll auf. Die beiden Entscheidungen vor dem Lauf (PolyG, Adapterfolgen samt Name des
bekannten Adapters) stehen mit Begründung im Protokoll.

## Validierung
Read für Read gegen fastp 1.4.0 ([Ergebnis](../../../bench/m9/ERGEBNIS.md)): HG001 NovaSeq
(PolyG an) und PJ MGI (PolyG aus), je 1.000.000 Paare. Alle behaltenen Reads gleich in Folge und
Qualität, keine Abweichung im Verbleib, gleiche Zahl gekürzter Adapterbasen. Die Laufzeiten dieses
Ausschnitts taugen nicht als Tempovergleich (siehe Ergebnis): fastp braucht für das ganze HG001
14 min 17 s, rund 2,4 s je Million Paare. Offen: Sätze, in denen eine
Adapterfolge erkannt wird (beide hier: keine), und die Wirkung am Ausrichten mit erzeugten Reads.
