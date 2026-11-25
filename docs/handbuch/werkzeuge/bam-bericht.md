# `genome bam-bericht`

Stand: gebaut und gegen samtools 1.24 geprüft. Baustein von [`fq2bam`](fq2bam.md), einzeln
aufrufbar.

## Zweck
Ein Überblick über eine nach Ort sortierte BAM: was in der Datei steckt (Kopf, Arten der
Datensätze, Zusatzfelder) und wie gut die Daten sind (Kartierung, MAPQ, Insertgröße, Fehlerrate,
Abdeckung). Als HTML-Seite zum Lesen, umschaltbar zwischen Deutsch und Englisch, und als
JSON-Datei für Programme. Ersetzt für diesen Zweck `samtools flagstat`, `samtools stats`,
`samtools coverage` und die Abdeckungsverteilung von mosdepth.

## Aufruf
```bash
genome bam-bericht probe.bam [--protokoll probe.bam.log.jsonl] [--aus probe.html] \
  [--kerne 0] [--log lauf.jsonl]
```

| Option | Bedeutung | Vorgabe |
|---|---|---|
| `--protokoll` | Lauf-Protokoll von `fq2bam`; daraus kommen Trimmen, Laufzeit, Rechner und Aufruf | `<bam>.log.jsonl`, wenn vorhanden |
| `--aus` | HTML-Seite; die JSON-Datei liegt daneben mit der Endung `.bamqc.json` | `<bam ohne .bam>.html` |
| `--kerne` | Kerne für das Entpacken der BAM (0 = alle) | 0 |
| `--log` | eigenes Lauf-Protokoll dieses Befehls | nur lesbar auf der Fehlerausgabe |

Beispiel für eine BAM, die schon fertig ist:

```bash
genome bam-bericht /Volumes/ExtremeSSD/genome4-bench/pj-voll/PJ.g4.bam
# schreibt PJ.g4.html und PJ.g4.bamqc.json neben die BAM
```

## Eingaben
Eine nach Ort sortierte BAM. Ein Index ist nicht nötig. Ohne Protokoll fehlen auf der Seite nur
die Abschnitte zum Lauf (Trimmen, Laufzeit, Duplikate markieren).

## Ausgaben
- `<name>.html`: eine einzige Datei ohne Netzverbindung, Schrift des Systems, hell und dunkel nach
  Systemeinstellung, druckbar. Die Sprache schaltet ein Knopf oben rechts um; die Wahl bleibt im
  Browser gespeichert, Vorgabe ist Deutsch. Die Seite enthält Kennzahlen, das Genom auf einen
  Blick (jedes Chromosom als Leiste in Feldern zu 1 Mio. Basen, gefärbt nach der Abdeckung im
  Verhältnis zum Mittel der Autosomen), die Abdeckung je Chromosom im Verhältnis zu den
  Autosomen (1,0 = zwei Kopien), jedes Chromosom im Detail (Kennzahlen, Abdeckung je 100 kb,
  Anteil MAPQ 0, Anteil Duplikate, Verteilung der Abdeckung; Werte beim Darüberfahren), Datei und Kopf (`@HD`,
  `@SQ` nach Art der Contigs, `@RG`, `@PG`), Datensätze nach Art, Zusatzfelder mit Bedeutung,
  Reads und Basen, MAPQ, Insertgröße, Abdeckung entlang des Genoms, Verteilung der Abdeckung,
  Tabelle je Chromosom, Geschlechtschromosomen und, mit Protokoll, wie die BAM entstanden ist.
- `<name>.bamqc.json`: alle Zahlen der Seite, Schema
  [`genome.bamqc.v1`](../formate/genome-bamqc-v1.md). Die Seite wird nur aus dieser Datei gesetzt.
- Beide Dateien entstehen unter vorläufigem Namen und werden am Ende umbenannt.
- JSON auf der Standardausgabe, Schema `genome.bam-bericht.v1`: `html`, `json`, `schema`,
  `mit_protokoll`, `lauf`.
- Protokoll: Stufe `bericht.lesen`, Ereignis `BERICHT.FERTIG`, Fehler `BERICHT.FEHLGESCHLAGEN`
  ([Ereignisse und Fehlercodes](../ereignisse-und-fehlercodes.md)).

## Verfahren
[Überblick über eine BAM](../verfahren/bam-ueberblick.md): ein Durchlauf über alle Datensätze,
die Abdeckung je Stelle ohne Feld über das ganze Genom.

## Validierung
HG002 chr20:10–12 Mb (GIAB, novoalign, 844.606 Datensätze) gegen samtools 1.24 (pysam):

| Zahl | samtools | Genome |
|---|---|---|
| `coverage`: Reads, abgedeckte Stellen, mittlere Abdeckung chr20 | 838.472 · 2.000.052 · 1,92157 | gleich |
| `depth -a`: Stellen mit ≥ 1, 10, 20, 30x | 2.000.052 · 1.999.228 · 1.997.578 · 1.992.359 | gleich |
| `depth -a`: Summe der Abdeckung | 123.834.008 | gleich (ausgerichtete Basen M) |
| `stats`: bases mapped (cigar) | 123.870.250 | 123.834.008 M + 36.242 I, gleich |
| `stats`: error rate (NM je M+I) | 0,3619755 % | gleich (`rate_wie_samtools`) |
| `stats`: average quality | 35,6 | 35,6 |
| `flagstat` | | gleich (gleiche Zählung wie [flagstat](flagstat.md)) |

Unit-Tests in `crates/gc-bam/src/qc.rs`: Abdeckung je Stelle mit Deletion, Duplikat und leerem
Contig, Fehlpaarungen aus NM, Tags, Insertgröße, Erkennung unsortierter Eingaben.

Laufzeit: auf 844.606 Datensätzen braucht der Überblick rund 0,15 s Rechenzeit mehr als der
Index allein, also grob 150 ns je Datensatz einschließlich Schreiben der Seite. Für eine BAM mit
30x (rund 730 Mio. Datensätze) sind das geschätzt 1 bis 2 Minuten; am Mac nicht gemessen.

## Grenzen
- Insertgröße: Paare mit TLEN über 2.000 stehen nur gezählt im letzten Fach; Mittel und
  Streuung gelten bis 2.000. `samtools stats` schneidet anders ab, die Werte weichen deshalb
  leicht ab (chr20-Ausschnitt: 555,5 ± 151,6 gegen 550,6 ± 143,8).
- Abdeckung zählt nur M, = und X; Deletionen zählen nicht als abgedeckt (wie `samtools depth`
  ohne `-J`). Ausgerichtete Stücke länger als 4 Mio. Basen würden gekürzt und gezählt
  (`gekuerzte_stuecke`); bei kurzen Reads kommt das nicht vor.
- Geschlechtschromosomen: nur eine Deutung aus der Abdeckung (X/A < 0,75 und Y/A > 0,2 passt zu
  XY; X/A > 0,8 und Y/A < 0,05 passt zu XX). Kein Ersatz für eine Untersuchung auf
  Aneuploidien.
- Kontamination, Verwandtschaft und Identität prüft der Überblick nicht (geplant in `gqc`).
- CIGAR mit mehr als 65.535 Operationen (Tag `CG`) werden nicht ausgewertet.

## Fehler und was zu tun ist
`BERICHT.FEHLGESCHLAGEN`: BAM nicht lesbar oder nicht nach Ort sortiert. Mit `samtools
quickcheck` prüfen, sonst erst `genome sort`. In `fq2bam` ist der Fehler nur eine Warnung; die
BAM bleibt gültig und der Überblick lässt sich mit diesem Befehl nachholen.

## Bericht
Die HTML-Seite ist der eigenständige Bericht zur BAM. Teile für den Gesamtbericht
(`genome.report.v1`) entstehen später aus der JSON-Datei.

## Literatur
Danecek et al. 2021 (samtools). Siehe [Literatur](../literatur.md).
