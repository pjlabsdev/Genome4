# BAM lesen und schreiben

Stand: gebaut und geprüft (`gc-io`, `bam.rs`; [Messung](../../bench/M1-bam-rundreise.md)).

## Grundlage
Die Crate `noodles` (reines Rust) für BAM und BGZF, Kompression über `flate2` mit `zlib-rs`
(reines Rust). BGZF-Blöcke werden über alle Kerne gepackt und entpackt.

## Kopf unverändert durchreichen
Beim Neuschreiben einer BAM baut noodles den Kopf aus seinen geparsten Feldern neu auf und ordnet
dabei die Zeilen nach Art (die `@PG`-Kette landet hinter `@SQ`). Die Datensätze bleiben gleich,
die Datei nicht. Für Reproduzierbarkeit und für den Vergleich mit samtools muss der Kopf Byte für
Byte erhalten bleiben. Deshalb gibt es zwei Wege:

- `open` / `create`: geparster Kopf, für Dateien, die Genome selbst neu erzeugt.
- `open_raw` / `create_raw`: der rohe Kopftext wird mitgelesen und unverändert geschrieben
  (Kennung, Textlänge, Text, Referenzliste von Hand nach der SAM-Spezifikation).

Der Test `roundtrip_keeps_header_line_order` hält das fest.

## Abschluss
`finish` schreibt die letzten Blöcke und die BGZF-Endmarke. Ohne sie gilt eine BAM als abgeschnitten.

## Grenzen
Noch kein Index (BAI/CSI), kein CRAM, kein Sortieren. Das kommt mit `fq2bam` (M5).

## Literatur
Li et al. 2009; SAM-Spezifikation. Siehe [Literatur](../literatur.md).
