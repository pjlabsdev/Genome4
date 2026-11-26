# `genome fq2bam`

Stand: erste lauffähige Fassung, auf chr20 ([Messung](../../bench/M5-fq2bam.md)) und auf dem
ganzen Genom mit erzeugten Reads ([Messung](../../bench/M5-ganzes-genom.md)) und mit echten
Reads auf chr20 ([Messung](../../bench/M5-echte-reads.md)) geprüft. Ein Lauf mit echten Reads in
voller Abdeckung, Trimmen, CRAM und die weitere Qualitätsauswertung sind offen.

## Zweck
Ein Programm vom FASTQ bis zur sortierten, duplikatmarkierten, indizierten BAM mit ersten
Qualitätszahlen. Ersetzt in Genome 3.x die Kette aus `bwa mem` 0.7.19 und samtools (fixmate,
sort, markdup, index, flagstat).

## Aufruf
```bash
genome fq2bam <index.gcfm> --fastq R1.fq[.gz] [--fastq2 R2.fq[.gz]] --aus probe.bam \
  [--lesegruppe '@RG\tID:lauf1\tSM:probe\tPL:ILLUMINA'] [--kerne 0] [--speicher-mb 2048] \
  [--temp ORDNER] [--block-basen 10000000] [--max-laenge 300] [--log DATEI] \
  [--trimmen] [--sa-ausduennen 0] [--schnellweg aus] [--ortsspeicher-mb 64] [--grafikeinheit aus] [--ohne-bericht] [--ohne-optisch] \
  [--optisch-pixel 100] [--optisch-mgi-abstand 1000]
```

| Option | Bedeutung | Vorgabe |
|---|---|---|
| `--fastq`, `--fastq2` | Reads; mit `--fastq2` als Paare | |
| `--aus` | Ausgabe-BAM | |
| `--lesegruppe` | Kopfzeile `@RG` wie bei `bwa mem -R`, Tabulatoren als `\t` | keine |
| `--kerne` | Kerne (0 = alle); die Ausgabe hängt nicht davon ab | 0 |
| `--speicher-mb` | Speicher für das Sortieren; die Ausgabe hängt nicht davon ab | 2048 |
| `--temp` | Ordner für die Teile des Sortierens | Ordner der Ausgabe |
| `--block-basen` | Paare: Basen je Block, aus dem die Insertgröße geschätzt wird | 10.000.000 |
| `--max-laenge` | längster erwarteter Read für das Markieren der Duplikate | 300 |
| `--log` | Lauf-Protokoll | `<aus>.log.jsonl` |
| `--trimmen` | Paare vor dem Ausrichten trimmen wie fastp 1.4.0 in der alten Kette (Adapter, PolyG, Qualität, Länge); verworfene Paare werden nicht ausgerichtet, siehe [Trimmen](../verfahren/trimmen.md). Nur für Paare | aus |
| `--ohne-optisch`, `--optisch-pixel`, `--optisch-mgi-abstand` | optische Duplikate zählen wie bei [`markdup`](markdup.md); ändert keine Flags | an, 100, 1000 |
| `--ohne-bericht` | keinen Überblick schreiben (wie vor dem Überblick) | aus |
| `--sa-ausduennen` | nur jede n-te Suffix-Array-Stichprobe des Index laden (Zweierpotenz). 0 wählt nach dem Arbeitsspeicher: die kleinste Zahl aus 1, 2 und 4, mit der Index und Sortierpuffer in die Hälfte des Arbeitsspeichers passen (`FQ2BAM.SPEICHERPLAN` im Protokoll). Mit 2 braucht der Index von GRCh38 (standard) rund 2 GB weniger Arbeitsspeicher und das Ausrichten etwa 2 % mehr Rechenzeit, mit 4 gut 20 % mehr. Die Ausgabe bleibt byte-gleich, siehe [FM-Index](../verfahren/fm-index.md#ausdünnen-beim-laden) | 0 |
| `--schnellweg` | `an`: leichte Reads ohne Seed-Suche ausrichten, auch neben schwachen weiteren Kopien, MAPQ gleich, `XS` kann sich ändern; `einfach`: wie `an`, aber nur Reads, deren Stücke je höchstens eine Fundstelle haben; `pruefen`: beide Wege rechnen und zählen, siehe [Schnellweg](../verfahren/schnellweg.md) | aus |
| `--ortsspeicher-mb` | Zwischenspeicher für das Orten der Seeds im Index, in MB; 0 schaltet ihn aus. Reads aus Wiederholungen treffen immer wieder dieselben Zeilen des Index, deren Position dann nicht neu durch die BWT gesucht wird. Die Ausgabe bleibt byte-gleich, siehe [FM-Index](../verfahren/fm-index.md#ortsspeicher) | 64 |
| `--grafikeinheit` | `an`: Die Grafikeinheit desselben Geräts rechnet die lokalen Ausrichtungen der Partnersuche, gesammelt je Block; die Ausgabe bleibt byte-gleich. Ohne nutzbare Grafikeinheit rechnet der Prozessor (`GPU.FEHLT`). Siehe [Grafikeinheit](../verfahren/grafikeinheit.md) | aus |

## Ablauf
1. Index laden ([FM-Index](../verfahren/fm-index.md)). Mit `--trimmen` läuft gleichzeitig die
   Adaptersuche am Anfang der FASTQ.
2. Ausrichten, Block für Block auf allen Kernen ([mehrere Kerne](../verfahren/mehrere-kerne.md),
   [Read-Paare](../verfahren/read-paare.md)), mit Lesegruppe und Partner-Angaben
   ([SAM-Ausgabe](../verfahren/sam-ausgabe.md)). Die Datensätze gehen ohne Zwischendatei in den
   Sortierer ([Sortieren](../verfahren/sortieren.md)); ist der Speicher voll, entsteht ein
   sortierter Teil auf der Platte. Drei Schritte laufen gleichzeitig: Ein Faden liest (und
   trimmt) den nächsten Block, alle Kerne richten den aktuellen aus, ein Faden schreibt den
   vorigen in den Sortierer. Die Blöcke sind dieselben wie nacheinander, die Ausgabe hängt
   nicht davon ab. `ALIGN.FERTIG` nennt die Sekunden je Schritt und wie lange das Ausrichten auf
   das Lesen gewartet hat.
3. Zusammenführen der Teile. Der sortierte Strom läuft durch das Markieren der Duplikate
   ([Duplikate](../verfahren/duplikate-markieren.md)) und die Zählung nach Flags
   ([flagstat](flagstat.md)) und wird als BAM geschrieben.
4. Index `.bai` in einem eigenen Durchlauf über die fertige BAM ([BAM-Index](../verfahren/bam-index.md)).
   Im selben Durchlauf entstehen die Zahlen für den Überblick
   ([Überblick über eine BAM](../verfahren/bam-ueberblick.md)); danach schreibt `fq2bam` die
   HTML-Seite und die JSON-Datei, mit Trimmen und Laufzeit aus dem eigenen Protokoll. Schlägt nur
   der Überblick fehl, bleibt die BAM gültig (Warnung `BERICHT.FEHLGESCHLAGEN`).

Auf der Platte liegen während des Laufs nur die Teile des Sortierens. Genome 3.x reicht die
Daten zwischen vier Programmen durch Pipes weiter und kodiert sie dabei mehrfach neu.

## Ausgaben
- `<aus>`: BAM, nach Ort sortiert, Duplikate mit Flag 1024, Kopf mit `SO:coordinate`, `@RG` und
  einer `@PG`-Zeile. Sie entsteht als `<aus>.unfertig` und bekommt ihren Namen erst am Ende.
- `<aus>.bai`.
- `<aus ohne .bam>.html` und `<aus ohne .bam>.bamqc.json`: Überblick über die BAM, Deutsch und
  Englisch umschaltbar ([`bam-bericht`](bam-bericht.md)). Fehlt mit `--ohne-bericht`.
- `<aus>.log.jsonl`: Lauf-Protokoll mit den Stufen `fq2bam.index_laden`, `fq2bam.ausrichten`,
  `fq2bam.zusammenfuehren`, `bamindex.lesen`, den Ereignissen `ALIGN.INSERTGROESSE` je Block,
  `SORT.TEIL` je Teil auf der Platte, `ALIGN.FERTIG`, `BAMINDEX.FERTIG`, `FQ2BAM.FERTIG` und im Fehlerfall `FQ2BAM.FEHLGESCHLAGEN`
  ([Ereignisse und Fehlercodes](../ereignisse-und-fehlercodes.md)).
- JSON, Schema `genome.fq2bam.v1`: `bam`, `bai`, `reads`, `saetze`, `paare`, `kerne`;
  `sortieren` (Bytes, Teile); `duplikate` (Felder wie bei [`markdup`](markdup.md)); `flagstat`
  (Felder wie bei [`flagstat`](flagstat.md)); `index`; `sekunden` und `sekunden_schritte`
  (`index_laden`, `ausrichten_und_teile`, `zusammenfuehren_markieren_schreiben`, `bam_index`);
  `reads_je_sekunde`; `bericht` (Pfade, siehe [`bam-bericht`](bam-bericht.md)); `lauf`.

## Validierung
chr20, 130.268 Paare mit 30.268 erzeugten Kopien (260.536 Datensätze):

- Dieselben Zeilen in derselben Reihenfolge wie die einzelnen Schritte `index sam`, `sort`,
  `markdup`; ebenso mit 1 statt 4 Kernen und mit 8 MB Speicher (10 Teile auf der Platte).
- Gegen die Kette `samtools fixmate -m`, `sort`, `markdup` 1.24: alle Zeilen gleich bis auf TLEN
  in 1.318 Zeilen ([Partner-Angaben](../verfahren/sam-ausgabe.md)), dieselben 56.067 Datensätze
  markiert.
- `samtools quickcheck` bestanden, `idxstats` mit dem eigenen Index gleich wie mit dem von
  samtools, Zahlen im Bericht gleich wie `samtools flagstat` und die Statistik von
  `samtools markdup`.

Die Übereinstimmung des Ausrichtens selbst mit bwa steht in den Messungen zu M3
([Verzeichnis](../messungen.md)).

## Grenzen
- Nur chr20 gemessen. Für das ganze Genom fehlen: Laufzeit und Speicherspitze, ALT-Contigs,
  Platzbedarf der Teile.
- Kein Trimmen von Adaptern (Genome 3.x: fastp), keine CRAM-Ausgabe, keine optischen Duplikate.
- Qualitätsauswertung: Zählung nach Flags, MAPQ, Insertgröße, Fehlerrate und Abdeckung im
  Überblick; Kontamination und Identität fehlen.
- Der Index braucht einen eigenen Durchlauf über die BAM.
- Kein Fortsetzen nach einem Abbruch.
- Eine Lesegruppe je Lauf, eine FASTQ-Datei je Richtung.

## Meilensteine

| M | Inhalt | Stand |
|---|---|---|
| M1 | FASTQ und BAM lesen und schreiben, Basismessung der heutigen Kette | gebaut; Auswertung der Basismessung folgt |
| M2 | Suffix-Array und FM-Index, Protokoll-Grundlage | gebaut und geprüft |
| M3 | Seeds, Ketten, Erweiterung, einzelne Reads und Paare auf chr20 | gebaut und gegen bwa geprüft |
| M4 | ergänzende Alignments im Vergleich, ALT, ganzes Genom | offen |
| M5 | Sortieren, Duplikate, Index, QC, Trimmen; erster Werkzeugbericht; Schattenbetrieb in 3.x | Sortieren, Duplikate, Index, Zählung nach Flags und der eine Befehl gebaut und geprüft; übrige QC, Trimmen, Bericht offen |
| M6 | eingefrorene Version, volle Testrunde, Forschungsbericht | offen |

## Vergleichsmaßstab
`bwa mem` 0.7.19 mit samtools 1.23.1, wie Genome 3.x sie heute aufruft
([Bestand](../../bench/M1-Bestand-Genome3-FASTQ-BAM.md)). Messziele: Plan Abschnitt 0.4.

## Literatur
Li 2013; Li & Durbin 2009; Vasimuddin et al. 2019; Danecek et al. 2021. Siehe
[Literatur](../literatur.md).
