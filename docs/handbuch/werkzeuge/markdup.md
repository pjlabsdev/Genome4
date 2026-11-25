# `genome markdup`

Stand: gebaut und geprüft ([Messung](../../bench/M5-duplikate.md)). Baustein von
[`fq2bam`](fq2bam.md), einzeln aufrufbar.

## Zweck
Markiert Duplikate: Reads oder Paare, die Kopien desselben Moleküls sind (aus der
Vervielfältigung bei der Vorbereitung der Probe). Sie tragen danach Flag 1024 und werden von
Varianten-Aufrufern nicht mitgezählt. Kein Datensatz wird entfernt oder sonst verändert.
Ersetzt `samtools markdup` in der Einstellung, die Genome 3.x benutzt.

## Aufruf
```bash
genome markdup sortiert.bam markiert.bam [--max-laenge 300] [--kerne 0] [--log lauf.jsonl] \
  [--ohne-optisch] [--optisch-pixel 100] [--optisch-mgi-abstand 1000]
```

Die Eingabe muss nach Ort sortiert sein ([`genome sort`](sort.md)) und bei Paaren die
Partner-Angaben `MC` und `ms` tragen. Die schreibt `genome index sam --partner-angaben`
([Partner-Angaben](../verfahren/sam-ausgabe.md)) oder `samtools fixmate -m`.

| Option | Bedeutung | Vorgabe |
|---|---|---|
| `--max-laenge` | längster erwarteter Read einschließlich abgeschnittener Basen (samtools `-l`) | 300 |
| `--ohne-optisch` | optische Duplikate nicht zählen; Zahlen und Bibliotheksgröße dann genau wie samtools ohne `-d` | aus |
| `--optisch-pixel` | Illumina: höchster Abstand in x und in y (Pixel) zum behaltenen Paar, damit ein Duplikat optisch zählt | 100 |
| `--optisch-mgi-abstand` | MGI/DNBSEQ: höchster Abstand der laufenden Nummer im selben Bildfeld | 1000 |
| `--kerne` | Kerne für das Entpacken und Packen der BAM (0 = alle) | 0 |
| `--log` | Lauf-Protokoll als JSON-Zeilen | nur lesbar auf der Fehlerausgabe |

## Ausgaben
- `markiert.bam`: dieselben Datensätze in derselben Reihenfolge, Duplikate mit Flag 1024. Kopf
  unverändert bis auf eine neue `@PG`-Zeile `ID:genome-markdup`. Die Datei entsteht als
  `markiert.bam.unfertig` und bekommt ihren Namen erst am Ende.
- JSON, Schema `genome.markdup.v1`: `gelesen`, `geschrieben`, `nicht_bewertet`, `bewertet`,
  `mit_partner`, `ohne_partner`, `duplikate_mit_partner`, `duplikate_ohne_partner`, `duplikate`,
  `anteil_duplikate` (an den bewerteten Reads), `optische_duplikate_mit_partner` und
  `anteil_optisch` (an den Duplikaten aus Paaren), `geschaetzte_bibliotheksgroesse` (verschiedene
  Moleküle, Formel aus Picard, optische Duplikate abgezogen; `null`, wenn nicht schätzbar),
  `fenster_spitze_saetze`, `max_laenge`, `optisch` (die Regel oder `null`), `sekunden`, `saetze_je_sekunde`, `lauf`. Gezählt werden Reads, nicht Paare, wie
  in der Statistik von samtools.
- Protokoll: Stufe `markdup.durchlauf`, Ereignis `MARKDUP.FERTIG`, Fehler
  `MARKDUP.FEHLGESCHLAGEN` ([Ereignisse und Fehlercodes](../ereignisse-und-fehlercodes.md)).

## Verfahren
[Duplikate markieren](../verfahren/duplikate-markieren.md).

## Validierung
130.268 Paare auf chr20 mit 30.268 erzeugten Kopien (260.536 Datensätze):

- Gegen `samtools markdup` 1.24 auf derselben Eingabe: alle 260.536 Zeilen gleich, 56.067
  markiert, alle zehn Zahlen der Statistik gleich, auch die geschätzte Bibliotheksgröße.
- Gegen die ganze Kette `samtools fixmate -m`, `sort`, `markdup`: alle Zeilen gleich bis auf
  TLEN in 1.318 Zeilen (bewusst der Wert von bwa, siehe
  [Partner-Angaben](../verfahren/sam-ausgabe.md)).
- Gegen die Wahrheit: Von 14.390 Gruppen von Kopien bleibt in 14.047 genau ein Paar, und immer
  das mit der besten Qualität. Von 84.638 Paaren ohne Kopie ist keines markiert. Die übrigen 343
  Gruppen sind auf der Messseite erklärt.

Optische Einstufung (Regel: [Verfahren](../verfahren/duplikate-markieren.md#optische-duplikate)):
Bauteiltests für MGI- und Illumina-Namen und dafür, dass die Flags mit und ohne Einstufung
gleich bleiben. Ganzes Genom von PJ (360 Mio. Paare, MGI): 129.700 von 5.391.316 Duplikaten aus
Paaren optisch (2,4 %), 695 s, 364 MB auf dem Mac M5; mit `--ohne-optisch` 736 s, 354 MB und
dieselbe Bibliotheksgröße wie zuvor (21,96 Milliarden). Die Sätze beider Ausgaben und der
Eingabe sind gleich (md5 über `samtools view`), der Kopf unterscheidet sich nur in der `@PG`-Zeile.

## Grenzen
- Welche Reads Flag 1024 bekommen, ist das Verfahren von samtools in der Vorgabe; die
  optische Einstufung ändert kein Flag und setzt keine Angabe `dt`. Nicht gebaut: Markieren zusätzlicher Zeilen und des Partners ohne Ort (`-S`), Strichcodes und UMI,
  getrennte Bewertung je Lesegruppe oder Bibliothek, Entfernen statt Markieren.
- Eine Lücke nahe am Ende eines Reads kann das Ende um wenige Basen verschieben; dann wird eine
  Kopie nicht erkannt (2 von 14.390 Gruppen). Das gilt für jedes Verfahren, das nach dem Ort
  urteilt.
- Kopien in Wiederholungen (MAPQ 0) landen oft an verschiedenen, gleich guten Orten und werden
  dann nicht erkannt (337 Gruppen).
- Ganzes Genom nicht gemessen (Laufzeit, Speicher bei sehr tiefer Abdeckung).
- Vorhandene Markierungen (Flag 1024 in der Eingabe) bleiben stehen, wie bei samtools ohne `-c`.
