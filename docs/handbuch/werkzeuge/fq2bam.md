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
  [--lesegruppe '@RG\tID:lauf1\tSM:probe\tPL:ILLUMINA'] [--kerne 0] [--speicher-mb 0] \
  [--temp ORDNER] [--block-basen 10000000] [--max-laenge 300] [--log DATEI] \
  [--trimmen] [--sa-ausduennen 0] [--schnellweg aus] [--ortsspeicher-mb 64] [--grafikeinheit aus] [--partnersuche-bereiche 20] [--verlaengerung-merken an] [--bam-stufe 3] [--strobeweg auto] [--strobe-datei PFAD] [--ohne-bericht] [--ohne-optisch] \
  [--optisch-pixel 100] [--optisch-mgi-abstand 1000]
```

| Option | Bedeutung | Vorgabe |
|---|---|---|
| `--fastq`, `--fastq2` | Reads; mit `--fastq2` als Paare | |
| `--aus` | Ausgabe-BAM | |
| `--lesegruppe` | Kopfzeile `@RG` wie bei `bwa mem -R`, Tabulatoren als `\t` | keine |
| `--kerne` | Kerne (0 = alle); die Ausgabe hängt nicht davon ab | 0 |
| `--speicher-mb` | Speicher für das Sortieren in MB; die Ausgabe hängt nicht davon ab. 0 wählt: 2048, mit Strobemer-Weg auf Geräten unter 24 GiB Arbeitsspeicher 512 (`FQ2BAM.SORTIERPUFFER` im Protokoll). Der Puffer belegt rund das 1,4-Fache seiner Größe; auf dem MacBook Air M5 mit 16 GB liefen alle Messungen mit Strobemer-Weg bei 512 ohne Auslagern | 0 |
| `--temp` | Ordner für die Teile des Sortierens | Ordner der Ausgabe |
| `--block-basen` | Paare: Basen je Block, aus dem die Insertgröße geschätzt wird | 10.000.000 |
| `--max-laenge` | längster erwarteter Read für das Markieren der Duplikate | 300 |
| `--log` | Lauf-Protokoll | `<aus>.log.jsonl` |
| `--trimmen` | Paare vor dem Ausrichten trimmen wie fastp 1.4.0 in der alten Kette (Adapter, PolyG, Qualität, Länge); verworfene Paare werden nicht ausgerichtet, siehe [Trimmen](../verfahren/trimmen.md). Nur für Paare | aus |
| `--ohne-optisch`, `--optisch-pixel`, `--optisch-mgi-abstand` | optische Duplikate zählen wie bei [`markdup`](markdup.md); ändert keine Flags | an, 100, 1000 |
| `--ohne-bericht` | keinen Überblick schreiben (wie vor dem Überblick) | aus |
| `--sa-ausduennen` | nur jede n-te Suffix-Array-Stichprobe des Index laden (Zweierpotenz). 0 wählt nach dem Arbeitsspeicher: die kleinste Zahl aus 1, 2, 4 und 8, mit der Index, Sortierpuffer und mit `--strobeweg an` der Strobe-Index (geschätzt 1,15 Bytes je Base) in die Hälfte des Arbeitsspeichers passen (`FQ2BAM.SPEICHERPLAN` im Protokoll). Auf dem MacBook Air M5 mit 16 GB (standard, 10 Mio. HG002-Paare, ohne Strobemer-Weg, 2026-11-27) war 2 am schnellsten: 478 s gegen 538 s mit 8 und 539 s mit 4; mit 1 (9,5 GB) lagerte macOS aus und der Lauf brauchte 683 s. Die Wahl nach Arbeitsspeicher trifft dort ohne Strobemer-Weg 2. Mit 2 braucht der Index von GRCh38 (standard) rund 2 GB weniger Arbeitsspeicher und das Ausrichten etwa 2 % mehr Rechenzeit, mit 4 gut 20 % mehr. Die Ausgabe bleibt byte-gleich, siehe [FM-Index](../verfahren/fm-index.md#ausdünnen-beim-laden) | 0 |
| `--schnellweg` | `an`: leichte Reads ohne Seed-Suche ausrichten, auch neben schwachen weiteren Kopien, MAPQ gleich, `XS` kann sich ändern; `einfach`: wie `an`, aber nur Reads, deren Stücke je höchstens eine Fundstelle haben; `pruefen`: beide Wege rechnen und zählen, siehe [Schnellweg](../verfahren/schnellweg.md) | aus |
| `--ortsspeicher-mb` | Zwischenspeicher für das Orten der Seeds im Index, in MB; 0 schaltet ihn aus. Reads aus Wiederholungen treffen immer wieder dieselben Zeilen des Index, deren Position dann nicht neu durch die BWT gesucht wird. Die Ausgabe bleibt byte-gleich, siehe [FM-Index](../verfahren/fm-index.md#ortsspeicher) | 64 |
| `--grafikeinheit` | `an`: Die Grafikeinheit desselben Geräts rechnet die lokalen Ausrichtungen der Partnersuche, gesammelt je Block; die Ausgabe bleibt byte-gleich. Ohne nutzbare Grafikeinheit rechnet der Prozessor (`GPU.FEHLT`). Siehe [Grafikeinheit](../verfahren/grafikeinheit.md) | aus |
| `--partnersuche-bereiche` | Für höchstens so viele Bereiche eines Reads wird der Partner gesucht. bwa nimmt 50 (`-m`), damit ist die Ausgabe bwa-gleich. Der Standard 20 spart Rechenzeit bei Reads aus Wiederholungen (Mac, 1 Mio. PJ-Paare: 33,3 statt 38,5 s) und ruft am ganzen chr20 von HG002 gegen GIAB gleich gute Varianten, ändert aber einzelne Zeilen (dort 0,75 %, bei PJ 1,7 %), meist MAPQ und Ort von Reads in Wiederholungen. Siehe [Read-Paare](../verfahren/read-paare.md#partnersuche-gesammelt-rechnen) | 20 |
| `--verlaengerung-merken` | Gleiche Verlängerungen eines Reads nur einmal rechnen. Das Ergebnis ist byte-gleich; `aus` rechnet jede neu wie bwa (Mac, 1 Mio. PJ-Paare: 32,9 statt 34,8 s). Siehe [Seed verlängern](../verfahren/seed-verlaengern.md#gleiche-verlängerungen-nur-einmal-rechnen) | an |
| `--bam-stufe` | Kompressionsstufe der BAM, 0 (keine) bis 9 (kleinste Datei). Vorgabe 3 (von pj am 26.11.2026 so gewählt); `--bam-stufe 6` schreibt wie samtools und wie fq2bam bisher. Kleinere Stufen schreiben schneller und größer; die Sätze sind bei jeder Stufe dieselben, nur die Bytes der Datei ändern sich. Gemessen (Simulation, 300.000 Paare, 4 Kerne): Zusammenführen 1,17 s bei 6, 0,94 s bei 3 (+3,3 % Größe), 0,63 s bei 1 (+40 %). Siehe [BAM schreiben](../verfahren/bam-schreiben.md#kompressionsstufe) | 3 |
| `--strobeweg` | `auto` (Vorgabe seit 2026-11-28, pjs Entscheidung): an, wenn das Gerät mindestens 14 GiB Arbeitsspeicher hat, sonst aus (`FQ2BAM.STROBEWEG` im Protokoll mit Grund). `an`: Reads, die der Schnellweg nicht nimmt, über Randstrobes verorten, wenn sie dort eindeutig und vollständig passen; alle anderen gehen den vollen Weg. Nicht bwa-gleich, gegen die Wahrheit geprüft: GIAB HG002 chr20 (60 Mio. Paare, MacBook Air M5, Schnellweg aus, Ausdünnen 8) SNV F1 0,8416 und Indel F1 0,7356 wie zuvor, segdups und lowmap besser als der bwa-gleiche Weg, 3.320 statt 4.081 s; Simulation über GRCh38 (2 Mio. Paare): falsche Reads mit MAPQ ≥ 20 74 statt 76, Ausrichten 116 statt 232 s. Mit Strobemer-Weg werden die Zählblöcke des FM-Index dicht angeordnet (1 GB weniger). `aus`: bwa-gleicher Weg. Siehe [Strobemer-Weg](../verfahren/strobeweg.md) | auto |
| `--strobe-datei PFAD` | Strobe-Index dort speichern und beim nächsten Start laden statt neu bauen (GRCh38: Bau rund 90 s auf 4 Kernen, 3,5 GB; Laden rund 1 s). Die Datei trägt einen Fingerabdruck der Referenz (Namen, Längen, Basen) und der Einstellungen; passt er nicht, wird neu gebaut und die Datei ersetzt. Eine beschädigte Datei bricht mit Fehler ab. Ohne Angabe gilt `<Index>.gcstrobe` neben dem Index; ist der Ordner nicht schreibbar, steht `FQ2BAM.STROBE_DATEI` als Warnung im Protokoll und der Lauf geht weiter. Nur mit Strobemer-Weg wirksam | `<Index>.gcstrobe` |

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

   Das Lesen selbst ist geteilt: Jede FASTQ wird in einem eigenen Faden entpackt und in Stücke
   von 8192 Reads zerlegt, die Paare eines Stücks werden auf allen Kernen getrimmt, danach
   schneidet ein Faden die Blöcke wie nacheinander nach den getrimmten Basen. Auf HG002
   (10 Mio. Paare, Mac M5) war vorher ein einzelner Faden fast die ganze Laufzeit mit Entpacken,
   Zerlegen und Trimmen beschäftigt, und das Ausrichten wartete 140 bis 167 s auf ihn.
   `GC_LESEN_PARALLEL=0` liest wie bisher in einem Faden. Die BAM ist in beiden Fällen dieselbe,
   auch mit Stückgrenzen mitten in einem Block (Test).
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
