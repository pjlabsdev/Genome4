# `genome sort`

Stand: gebaut und geprüft ([Messung](../../bench/M5-sortieren.md)). Baustein von
[`fq2bam`](fq2bam.md), einzeln aufrufbar.

## Zweck
Ordnet die Datensätze einer BAM nach dem Ort auf dem Referenzgenom. Die Reihenfolge ist die von
`samtools sort`. Der Speicherbedarf ist fest: Was nicht in den Speicher passt, liegt sortiert in
Teilen auf der Platte.

## Aufruf
```bash
genome sort ein.bam aus.bam [--speicher-mb 2048] [--temp ORDNER] [--kerne 0] [--log lauf.jsonl]
```

| Option | Bedeutung | Vorgabe |
|---|---|---|
| `--speicher-mb` | Speicher für die Datensätze. Dazu kommen 16 Bytes je Datensatz im Speicher und die Puffer für das Lesen und Schreiben | 2048 |
| `--temp` | Ordner für die Teile | Ordner der Ausgabe |
| `--kerne` | Kerne für das Entpacken und Packen der BAM (0 = alle) | 0 |
| `--log` | Lauf-Protokoll als JSON-Zeilen | nur lesbar auf der Fehlerausgabe |

## Ausgaben
- `aus.bam`. Sie entsteht als `aus.bam.unfertig` und bekommt ihren Namen erst nach dem letzten
  Byte. Nach einem Abbruch gibt es keine Datei unter dem Zielnamen.
- Kopf: `@HD` mit `SO:coordinate` (ein Feld `SS` wird entfernt, fehlt `@HD`, kommt
  `@HD VN:1.6 SO:coordinate` an den Anfang). Alle anderen Kopfzeilen bleiben Byte für Byte
  gleich. Am Ende steht eine neue `@PG`-Zeile `ID:genome-sort` mit dem Aufruf, die mit `PP` auf
  die letzte vorhandene verweist.
- JSON, Schema `genome.sort.v1`: `saetze`, `bytes_saetze` (unkomprimiert, mit Längenfeldern),
  `teile` und `bytes_teile` (Teile auf der Platte), `speicher_mb`, `sekunden`,
  `sekunden_lesen`, `saetze_je_sekunde`, `lauf`.
- Protokoll: Stufen `sort.lesen` und `sort.schreiben`, Ereignisse `SORT.TEIL` und `SORT.FERTIG`,
  Fehler `SORT.FEHLGESCHLAGEN` ([Ereignisse und Fehlercodes](../ereignisse-und-fehlercodes.md)).

## Verfahren
[Sortieren nach Ort mit festem Speicher](../verfahren/sortieren.md). Format der Teile:
[`*.gcsort`](../formate/gcsort.md).

## Validierung
Gegen `samtools sort` 1.24 (über pysam), Datensatz für Datensatz über `samtools view`:

- 100.000 Paare auf chr20 (200.000 Datensätze), in der Reihenfolge des Ausrichtens und zufällig
  gemischt,
- 300.000 erzeugte Datensätze auf 7 Contigs mit vielen gleichen Orten, beiden Strängen und
  Datensätzen ohne Contig,

jeweils mit 2048, 16 und 2 MB Speicher (0 bis 32 Teile): in allen neun Läufen dieselbe Folge wie
samtools, dieselbe Menge wie die Eingabe, `quickcheck` bestanden.

## Grenzen
- Nur BAM als Eingabe und Ausgabe, nur Sortierung nach Ort (nicht nach Namen).
- Bisher nur an kleinen Dateien gemessen (bis 68 MB Datensätze). Ein ganzes Genom mit rund
  100 GB und der Speicherbedarf dabei sind nicht gemessen.
- Das Sortieren im Speicher läuft auf einem Kern.
- Kein Fortsetzen nach einem Abbruch: Der Lauf beginnt von vorn. Wird das Programm hart beendet
  (Stromausfall, `kill -9`), bleibt der Ordner `genome-sort-<Prozess>-<Zahl>` im Temp-Ordner
  liegen und kann gelöscht werden.
- Je Teil bleibt beim Zusammenführen eine Datei offen. Bei sehr kleinem Speicher und sehr großer
  Eingabe kann die Grenze des Betriebssystems für offene Dateien erreicht werden (Beispiel:
  100 GB mit 2048 MB sind rund 50 Teile, das ist unkritisch).
- Schreibt noch keinen Index (`.bai`, `.csi`).
