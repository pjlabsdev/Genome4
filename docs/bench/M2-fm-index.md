# M2: Suffix-Array und FM-Index

Stand 18.11.2026. Teil 1 von Meilenstein M2 (Plan Abschnitt 0.5): Verfahren, Format und Prüfung
an einem Chromosom. Der Bau für das ganze Genom folgt als Teil 2 (siehe "Offen").

## Was gebaut ist

- `gc-index::sais`: Suffix-Array nach SA-IS (lineare Zeit), dazu eine unabhängige Prüfung in
  linearer Zeit, die für jedes Suffix-Array beweist, dass es richtig sortiert ist.
- `gc-index::reference`: Referenz als 2-Bit-Folge mit Contig-Tabelle und Liste der N-Strecken.
  N wird wie bei bwa mit festen Pseudozufallsbasen gefüllt, die Strecken bleiben vermerkt.
- `gc-index::fm`: FM-Index über beide Stränge (Vorwärtsstrang + reverses Komplement, wie bwa).
  Zwei Varianten desselben Formats:

  | Variante | Zählstand | Suffix-Array-Stichprobe | Bits je Referenzbase | ganzes Genom (3,1 Mrd. Basen), gerechnet |
  |---|---|---|---|---|
  | kompakt | alle 128 Basen | jede 32. Zeile | 12 | rund 4,7 GB |
  | schnell | alle 64 Basen | jede 8. Zeile | 22 | rund 8,5 GB |

  Die Zahlen enthalten beide Stränge und die gepackte Referenz. Die Größe für das ganze Genom ist
  hochgerechnet, nicht gemessen.
- Datei `*.gcfm`: little-endian, mit Formatversion und CRC-32; halbe oder beschädigte Dateien
  werden beim Lesen erkannt. Geschrieben wird erst in eine Nachbardatei, dann umbenannt.
- `gc-log`: jeder Bau schreibt ein Lauf-Protokoll `<index>.log.jsonl` (Plan Abschnitt 15).

## Prüfung auf chr20 (GRCh38, 64.444.167 Basen, davon 499.910 N in 90 Strecken)

Rechner: Linux-Container, 4 Kerne, 15 GB. Referenz: GIAB `GCA_000001405.15_GRCh38_no_alt_analysis_set`.

| Prüfung | Ergebnis |
|---|---|
| Suffix-Array (128.888.335 Zeilen) vollständig geprüft | richtig sortiert, jede Position genau einmal |
| 1.000.000 Zufallsabfragen (Länge 1 bis 150, ein Drittel mit einer Änderung) gegen binäre Suche im vollen Suffix-Array | 0 Abweichungen; 17.939.021 Fundstellen einzeln verglichen |
| 1.000.000 Abfragen (Länge 20 bis 100, beide Stränge, ein Viertel mit einer Änderung) gegen bwa (`bwa aln`, exakt, Tag X0) | 999.935 gleiche Trefferzahl, 65 um genau 1 verschieden |
| Positionen von 4.988 dieser Abfragen gegen bwa | alle gleich |
| Die 65 Abweichungen gegen direkte Suche in der FASTA | alle erklärt, siehe unten |

**Die 65 Abweichungen:** Verglichen wurde die rohe Trefferzahl (Breite des Suffix-Array-Bereichs).
Sie zählt auch Fundstellen, die in eine N-Strecke hineinreichen. Dort stehen Zufallsbasen, bei bwa
andere als bei uns. Bei häufigen Wiederholungen direkt neben einer Lücke passt so ein Stück
gelegentlich zufällig. Gegen die direkte Suche in der FASTA: bwa zählt in 43 Fällen eine Fundstelle
zu viel, Genome roh in 24 Fällen, keiner zählt je zu wenig. Die **gefilterten** Genome-Treffer
(`hits`, ohne Fundstellen in N-Strecken) stimmen in allen 65 Fällen mit der FASTA überein.

Zeiten (ein Lauf, nicht als Tempo-Vergleich gedacht): Bau kompakt ohne Selbstprüfung rund 40 s,
davon Suffix-Array 22 s und FASTA lesen 9 s. 1 Mio. Zählabfragen: Genome 10,8 s auf einem Kern
inklusive Laden des Index, bwa aln 9,9 s auf vier Kernen. Die Suche ist noch nicht optimiert.

## Nachrechnen

```bash
genome index build chr20.fa chr20.gcfm --variante kompakt --pruefen 1000000
python bench/m2/compare_bwa.py   target/release/genome chr20.gcfm chr20.fa 1000000 vergleich.json
python bench/m2/explain_diffs.py target/release/genome chr20.gcfm chr20.fa 1000000 erklaert.json
```

bwa kommt über das Python-Paket `pybwa` 2.3.2 (enthält den Original-C-Code von bwa), weil in der
Arbeitsumgebung kein conda erreichbar ist. Ergebnisse und Lauf-Protokoll: `bench/m2/ergebnisse/`.

## Offen (Teil 2 von M2)

- **Ganzes Genom:** Die Positionen sind in dieser Formatversion 32 Bit breit; beide Stränge des
  Humangenoms (6,2 Mrd. Zeichen) passen nicht hinein. Außerdem bräuchte das volle Suffix-Array im
  Speicher rund 30 GB. Nötig sind 40-Bit-Positionen und ein blockweiser Bau, der mit festem
  RAM-Budget auskommt (Ziel: Bau auf einem Rechner mit 16 GB). Das Programm meldet den Fall heute
  mit klarer Fehlermeldung.
- Laden per Speicherabbild über `gc-sys` statt vollständigem Einlesen.
- Schnelle Suche für die Seeds (Vorabruf, mehrere Abfragen verschränkt) gehört zu M3.
