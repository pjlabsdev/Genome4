# M1: BAM-Rundreise gegen samtools

Stand 18.11.2026, Branch `m1-samtools-roundtrip`. Frage: Schreibt `gc-io` eine BAM so, dass
samtools darin exakt dieselben Daten liest wie im Original?

## Aufbau

- **Daten:** GIAB HG002, Illumina HiSeq 60x, GRCh38, Ausschnitt chr20:10.000.000-12.000.000
  (844.606 Datensätze, 78 MB), geholt mit
  `samtools view -b <GIAB-S3>/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam chr20:10000000-12000000`.
- **Prüfer:** samtools 1.24 / htslib über pysam (fremdes Werkzeug, unabhängig vom Geprüften).
  `bench/m1/compare_bam.py` bildet SHA-256 über Kopftext und jeden Datensatz als SAM-Zeile.
- **Lauf:** `genome bam-roundtrip hg002_chr20.bam rt.bam -t 4`, Linux x86-64, 4 Kerne.

## Ergebnis

| Prüfung | Ergebnis |
|---|---|
| Kopftext identisch | ja |
| Datensätze identisch | 844.606 von 844.606 |
| `samtools quickcheck` auf unserer BAM | ok |
| `samtools index` auf unserer BAM, Regionenabfrage chr20:10,5-10,6 Mb | 41.505 Datensätze, wie im Original |

## Gefundener Fehler

Der erste Lauf war **nicht** identisch: Die Datensätze stimmten, aber die Bibliothek (noodles)
baut den Kopf beim Schreiben neu auf und sortiert die Zeilen nach Art (`@HD`, `@SQ`, `@RG`,
`@PG`, `@CO`). Im Original stand die `@PG`-Kette vor den `@SQ`-Zeilen. Inhaltlich gleich, aber
nicht dieselbe Datei. Behoben: `open_raw`/`create_raw` reichen den Kopftext Byte für Byte durch;
Test `roundtrip_keeps_header_line_order` sichert das ab.

## Was diese Messung nicht zeigt

- **Tempo:** Im selben Lauf brauchte `genome` 1,8 s und `samtools view -b -@4` 2,5 s (Median aus 3).
  Das ist **kein fairer Vergleich**: das samtools aus pysam nutzt sehr wahrscheinlich zlib statt
  libdeflate, und unsere Datei ist 1,3 % größer (79,4 gegen 78,4 MB), also schwächer komprimiert.
  Die gültige Tempo-Messung kommt mit pj's samtools auf dem Mac bei gleicher Dateigröße.
- **Breite:** eine Probe, ein Aligner (novoalign), nur gemappte Lesungen einer Region. Offen:
  ungemappte Lesungen, sehr lange CIGARs, BWA-Ausgaben mit `SA`/`XA`-Tags, MGI-Daten, CRAM.
