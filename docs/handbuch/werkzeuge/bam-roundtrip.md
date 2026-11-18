# `genome bam-roundtrip`

Stand: gebaut und geprüft ([Messung](../../bench/M1-bam-rundreise.md)).

## Zweck
Liest eine BAM und schreibt sie unverändert neu. Das ist kein Werkzeug für den Alltag, sondern der
Beleg, dass Lesen und Schreiben von `gc-io` verlustfrei und zu samtools gleichwertig sind.

## Aufruf
```bash
genome bam-roundtrip ein.bam aus.bam [-t KERNE]   # 0 = alle Kerne (Standard)
```

## Ausgaben
- `aus.bam`
- JSON, Schema `genome.bam-roundtrip.v1`: `records` (Zahl der Datensätze), `seconds`.

## Verfahren
[BAM lesen und schreiben](../verfahren/bam-schreiben.md).

## Validierung
GIAB HG002, Illumina 60x, chr20:10-12 Mb, 844.606 Datensätze: Kopftext identisch, alle Datensätze
identisch, `samtools quickcheck`, `samtools index` und eine Regionenabfrage ohne Befund. samtools
(über pysam) dient als unabhängiger Prüfer.

## Grenzen
- Eine Probe, ein Aligner, nur gemappte Lesungen einer Region. Offen: ungemappte Lesungen, MGI, CRAM.
- Die Zeiten im Messdokument sind kein fairer Tempo-Vergleich.
- Die komprimierte Datei ist nicht bytegleich mit der von samtools (anderer Kompressor), der
  Inhalt ist es.
