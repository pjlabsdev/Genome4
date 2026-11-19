# M5: Duplikate markieren

Stand: 19.11.2026. Verfahren: [Duplikate markieren](../handbuch/verfahren/duplikate-markieren.md),
[Partner-Angaben](../handbuch/verfahren/sam-ausgabe.md). Werkzeug:
[`genome markdup`](../handbuch/werkzeuge/markdup.md).

## Fragen
1. Schreibt `genome index sam --partner-angaben` dieselben Angaben wie `samtools fixmate -m`?
2. Markiert `genome markdup` dieselben Datensätze wie `samtools markdup`?
3. Findet das Verfahren die Kopien, und bleibt das beste Paar?

## Aufbau
- Referenz chr20 aus GRCh38, Index `chr20.schnell.gcfm`. Prüfer: samtools 1.24 aus pysam 0.24.1.
- Frage 1: die 100.000 Paare aus [M3 Paare](M3-paare.md).
- Fragen 2 und 3: `bench/m5/make_duplicates.py` macht aus diesen Paaren 130.268: jedes Paar
  bekommt zufällige Basenqualitäten, 15 % bekommen 1 bis 3 Kopien (30.268 Kopien). Von den Kopien
  behält ein Viertel die Qualitäten (Gleichstand), ein Viertel hat bis zu drei Lesefehler je
  Read, bei einem Zehntel ist der zweite Read durch Zufallsbasen ersetzt. Die Reihenfolge ist
  gemischt.
- `bench/m5/compare_markdup.py`, Ergebnis `bench/m5/markdup_vergleich.json`.

## Ergebnis

### 1. Partner-Angaben gegen `samtools fixmate -m`
200.000 Zeilen: Alle Zusatzangaben sind gleich (als Menge je Zeile; fixmate stellt `MQ` und `MC`
ans Ende). Alle anderen Felder sind gleich bis auf TLEN in 1.018 Zeilen:

| Lage des Paares | Zeilen | Unterschied im Betrag |
|---|---|---|
| Reads zeigen voneinander weg (RF) | 1.012 | fixmate 2 kleiner |
| beide auf demselben Strang | 6 | fixmate 1 kleiner |

Keines dieser Paare gilt als richtiges Paar (Flag 2). bwa und samtools rechnen die Länge für
solche Lagen verschieden; Genome 4.0 behält hier den Wert von bwa. Für das Markieren der
Duplikate spielt TLEN keine Rolle.

### 2. Gegen `samtools markdup`

| Vergleich | Zeilen | gleich | markiert Genome 4.0 | markiert samtools |
|---|---|---|---|---|
| A: dieselbe sortierte BAM in beide Programme | 260.536 | 260.536 (alle Felder) | 56.067 | 56.067 |
| B: ganze Kette (`fixmate -m`, `sort`, `markdup` von samtools) | 260.536 | 260.536 (ohne TLEN und Reihenfolge der Zusatzangaben) | 56.067 | 56.067 |

Statistik in A (samtools `-f --json` gegen `genome.markdup.v1`), alle gleich: gelesen 260.536,
geschrieben 260.536, nicht bewertet 4.292, bewertet 256.244, mit Partner 251.952, ohne Partner
4.292, Duplikate mit Partner 52.830, ohne Partner 3.237, zusammen 56.067, geschätzte
Bibliotheksgröße 256.745.

In B unterscheidet sich TLEN in 1.318 Zeilen (wie unter 1).

### 3. Gegen die Wahrheit
Betrachtet werden Paare, bei denen beide Reads einen Ort haben.

| | Gruppen |
|---|---|
| Paare ohne Kopie | 84.638, davon markiert: 0 |
| Gruppen mit Kopien | 14.390 |
| genau ein Paar bleibt | 14.047 (97,6 %) |
| davon bleibt ein Paar mit der höchsten Qualitätssumme | 14.047 |
| mehrere bleiben | 342 |
| keines bleibt | 1 |

Die 342 Gruppen, in denen mehr als ein Paar unmarkiert bleibt:

- 337: Die Paare liegen in Wiederholungen (MAPQ 0), und die Kopien wurden an verschiedene, gleich
  gute Orte gelegt.
- 3: verschiedener Ort bei MAPQ über 0 (Kopien mit Lesefehlern).
- 2: gleicher Ort, aber eine Lücke von 3 Basen kurz vor dem Ende des Reads. Die Kopie ohne
  Lesefehler ist dort abgeschnitten (`144M6S`), die andere hat die Lücke (`141M3D9M`); das Ende
  liegt 3 Basen weiter.

Die eine Gruppe ohne unmarkiertes Paar: An denselben Enden liegt ein weiteres, unabhängig
erzeugtes Paar mit besserer Qualität. Das Verfahren kann die beiden Moleküle nicht unterscheiden.

## Zeiten
Nur zur Einordnung (kleine Datei, ein Lauf, samtools innerhalb von Python): `genome markdup`
1,8 s, `samtools markdup` 4,9 s für 260.536 Datensätze.

## Grenzen
- Erzeugte Reads und erzeugte Kopien auf einem Chromosom. Echte Daten mit echten Duplikaten,
  zusätzlichen Zeilen in größerer Zahl und Paaren über zwei Chromosomen sind nicht gemessen.
- Die Abdeckung ist gering (rund 0,6-fach). Speicher und Laufzeit bei 30-facher Abdeckung sind
  nicht gemessen.
- Optische Duplikate und `-S` sind nicht gebaut und nicht verglichen.
- Vergleich nur mit samtools, nicht mit Picard `MarkDuplicates`.

## Wiederholen
```bash
python bench/m5/make_duplicates.py pairs_1.fq pairs_2.fq d_1.fq d_2.fq wahrheit.tsv
python bench/m5/compare_markdup.py target/release/genome chr20.schnell.gcfm \
  d_1.fq d_2.fq wahrheit.tsv arbeit bench/m5/markdup_vergleich.json
```
