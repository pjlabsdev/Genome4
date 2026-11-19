# `genome bam-index`

Stand: gebaut und geprüft ([Messung](../../bench/M5-index.md)). Baustein von
[`fq2bam`](fq2bam.md), einzeln aufrufbar.

## Zweck
Schreibt zu einer nach Ort sortierten BAM den Index `.bai`. Mit ihm findet ein Programm die
Datensätze eines Bereichs, ohne die ganze Datei zu lesen. Ersetzt `samtools index`.

## Aufruf
```bash
genome bam-index sortiert.bam [--aus DATEI] [--kerne 0] [--log lauf.jsonl]
```

| Option | Bedeutung | Vorgabe |
|---|---|---|
| `--aus` | Name der Indexdatei | `<bam>.bai` |
| `--kerne` | Kerne für das Entpacken der BAM (0 = alle) | 0 |
| `--log` | Lauf-Protokoll als JSON-Zeilen | nur lesbar auf der Fehlerausgabe |

## Ausgaben
- `<bam>.bai` im Format der SAM-Spezifikation (Abschnitt 5.2). Die Datei entsteht unter
  vorläufigem Namen und wird am Ende umbenannt.
- JSON, Schema `genome.bam-index.v1`: `saetze`, `ohne_contig`, `contigs`, `contigs_mit_daten`,
  `faecher`, `bytes_index`, `sekunden`, `lauf`.
- Protokoll: Stufe `bamindex.lesen`, Ereignis `BAMINDEX.FERTIG`, Fehler
  `BAMINDEX.FEHLGESCHLAGEN` ([Ereignisse und Fehlercodes](../ereignisse-und-fehlercodes.md)).

## Verfahren
[BAM-Index](../verfahren/bam-index.md).

## Validierung
Gegen `samtools index` 1.24 an drei Dateien (chr20 mit 260.536 und 200.000 Datensätzen, 7
Contigs mit 300.000): Fächer mit allen Dateistücken, linearer Index, Zähler und Zahl ohne Contig
sind gleich. Je Datei 2.000 zufällige Abfragen über htslib, einmal mit jedem Index: dieselben
Datensätze (zusammen über 16 Millionen gelieferte Datensätze), `idxstats` gleich.

## Grenzen
- Nur `.bai`: Contigs bis 512 Mbp (für den Menschen ausreichend). Das Format `.csi` für längere
  Contigs ist nicht gebaut.
- Der Index entsteht in einem eigenen Durchlauf über die fertige BAM. Schneller wäre, ihn beim
  Schreiben mitzuführen; das braucht einen eigenen BGZF-Schreiber und ist offen.
- Die Datei ist inhaltlich gleich, aber nicht immer Byte für Byte: htslib schreibt die Fächer in
  der Reihenfolge seiner Hashtabelle, Genome 4.0 aufsteigend nach Nummer.
- Ganzes Genom nicht gemessen.
