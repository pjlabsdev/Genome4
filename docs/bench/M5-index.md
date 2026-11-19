# M5: BAM-Index

Stand: 19.11.2026. Verfahren: [BAM-Index](../handbuch/verfahren/bam-index.md). Werkzeug:
[`genome bam-index`](../handbuch/werkzeuge/bam-index.md).

## Fragen
1. Hat der Index denselben Inhalt wie der von `samtools index`?
2. Liefern Abfragen mit ihm dieselben Datensätze?

## Aufbau
- Prüfer: samtools 1.24 und htslib aus pysam 0.24.1.
- Dateien:
  - `idx_chr20`: die markierte BAM aus [M5 Duplikate](M5-duplikate.md), 260.536 Datensätze.
  - `idx_7contigs`: die 300.000 erzeugten Datensätze auf 7 Contigs aus
    [M5 Sortieren](M5-sortieren.md), sortiert von `genome sort`.
  - `idx_samtools_sortiert`: 200.000 Datensätze auf chr20, sortiert und gepackt von samtools
    (andere Blockgrenzen als bei Genome 4.0).
- `bench/m5/compare_index.py`: beide `.bai` zerlegen und vergleichen; dann je Datei 2.000
  zufällige Bereiche (1 Base bis 5 Mbp) über htslib lesen, einmal mit jedem Index, und die
  gelieferten Datensätze vergleichen; dazu `idxstats` und die Datensätze ohne Contig.
  Ergebnis: `bench/m5/index_vergleich.json`.

## Ergebnis

| Datei | Datensätze | Fächer | Fächer und Stücke | linearer Index | Abfragen verschieden | gelieferte Datensätze | Datei Byte für Byte |
|---|---|---|---|---|---|---|---|
| chr20, markiert | 260.536 | 384 | gleich | gleich | 0 von 2.000 | 8.524.230 | nein |
| 7 Contigs | 300.000 | 14 | gleich | gleich | 0 von 2.000 | 1.020.970 | ja |
| chr20, von samtools sortiert | 200.000 | 375 | gleich | gleich | 0 von 2.000 | 6.548.202 | nein |

In allen drei Dateien sind auch die Zähler je Contig, die Zahl ohne Contig und `idxstats` gleich.
Wo die Dateien nicht Byte für Byte gleich sind, liegt es nur an der Reihenfolge der Fächer in der
Datei (htslib: Reihenfolge der Hashtabelle; Genome 4.0: aufsteigend).

Ein Fund beim Bau: Die Stelle "hinter dem letzten Datensatz" darf nicht nach dem Dateiende
abgelesen werden, sonst liegt sie hinter dem leeren Schlussblock der BAM, 28 Bytes zu weit. Das
betraf nur das Ende des letzten Dateistücks und ist behoben.

## Grenzen
- Kleine Dateien, höchstens 7 Contigs. Ein ganzes Genom mit seinen vielen Contigs und sehr
  tiefer Abdeckung ist nicht gemessen.
- Laufzeit nicht gemessen.

## Wiederholen
```bash
python bench/m5/compare_index.py target/release/genome arbeit bench/m5/index_vergleich.json \
  idx_chr20.bam idx_7contigs.bam idx_samtools_sortiert.bam
```
