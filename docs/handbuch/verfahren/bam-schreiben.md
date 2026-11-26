# BAM lesen und schreiben

Stand: gebaut und geprüft (`gc-io`, `bam.rs`; [Messung](../../bench/M1-bam-rundreise.md);
Ausgabe des Aligners: `gc-align`, `sam.rs`, [Messung](../../bench/M3-kerne-bam.md)).

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

## Zeilen des Aligners direkt als BAM
`genome index sam --bam` schreibt die Zeilen des Aligners ohne Umweg über SAM-Text. Jede Zeile
wird in `gc-align` (`Line::write_bam`) nach Abschnitt 4.2 der SAM-Spezifikation kodiert und als
fertiger Satz an den BGZF-Schreiber gegeben (`write_encoded`). Der Kopf ist derselbe Text wie
bei der SAM-Ausgabe (`create_from_text`).

| Feld | Kodierung |
|---|---|
| Position, Partnerposition | 0-basiert; -1 ohne Ort |
| Behälter (`bin`) | aus Position und Länge auf der Referenz (Abschnitt 5.3); ohne Ort 4680 wie bei samtools |
| CIGAR | je Schritt 32 Bit: Länge und Art (M, I, D, S, H) |
| Basen | 4 Bit je Base (A=1, C=2, G=4, T=8, N=15) |
| Qualitäten | Phred-Wert (ASCII minus 33); ohne Qualitäten 0xFF |
| Zusatzangaben | Text als `Z`; ganze Zahlen im kleinsten passenden Typ, wie samtools |

Ein Read-Name über 254 Zeichen passt nicht in das Längenfeld; der Lauf endet dann mit einer
Meldung statt einer beschädigten Datei.

Geprüft wird die Kodierung von außen: samtools liest die BAM zurück, und das Ergebnis ist Zeile
für Zeile die SAM-Ausgabe desselben Laufs ([Messung](../../bench/M3-kerne-bam.md)). Im CI
macht dasselbe ein Test mit noodles als Leser (`crates/genome/tests/ausrichten.rs`).

## Abschluss
`finish` schreibt die letzten Blöcke und die BGZF-Endmarke. Ohne sie gilt eine BAM als abgeschnitten.

## Kompressionsstufe
Jeder BGZF-Block wird mit Deflate komprimiert, auf allen Kernen zugleich. Die Stufe ist
`DEFAULT_LEVEL` = 6 wie bei samtools; `create_raw_level` und `create_from_text_level` nehmen eine
andere (0 bis 9). `fq2bam` schreibt mit Stufe 3 (`--bam-stufe`, 6 gibt die bisherige Datei). Die Stufe ändert nur die Bytes der Datei, gelesen
ergeben sich dieselben Sätze in derselben Reihenfolge; der Test `level_changes_bytes_not_records`
prüft das und dass Stufe 6 dieselbe Datei ergibt wie die Vorgabe. Im Zusammenführen von `fq2bam`
ist das Komprimieren der größte Teil der Zeit. Gemessen (Simulation über GRCh38, 300.000 Paare,
4 Kerne, 6 Sortierteile):

| Stufe | Zusammenführen | Größe |
|---|---|---|
| 6 | 1,17 s | 47,3 MB |
| 3 | 0,94 s | 48,9 MB (+3,3 %) |
| 1 | 0,63 s | 66,3 MB (+40 %) |

## Grenzen
Kein CRAM.

## Literatur
Li et al. 2009; SAM-Spezifikation. Siehe [Literatur](../literatur.md).
