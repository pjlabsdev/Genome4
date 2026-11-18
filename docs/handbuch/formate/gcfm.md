# Indexdatei `*.gcfm`

Stand: gebaut (`gc-index`, `fm.rs`), Formatversion 2.

Alle Zahlen little-endian. `u64` = 8 Bytes, `u32` = 4 Bytes.

| Nr. | Inhalt | Typ |
|---|---|---|
| 1 | Kennung `GCFM` | 4 Bytes |
| 2 | Formatversion (2) | u32 |
| 3 | Zeilen der BWT (2L + 1) | u64 |
| 4 | Primärzeile | u64 |
| 5 | Basen je Zählblock `b` (Vielfaches von 32) | u64 |
| 6 | Abstand der Suffix-Array-Stichproben `s` | u64 |
| 7 | `c_before` für A, C, G, T | 4 × u64 |
| 8 | Zahl der Contigs, dann je Contig: Länge des Namens, Name (UTF-8), Länge, Anfang im Vorwärtsstrang | u64, dann u64 + Bytes + u64 + u64 |
| 9 | Zahl der N-Strecken, dann je Strecke: Anfang, Länge | u64, dann 2 × u64 |
| 10 | Länge der Referenz L | u64 |
| 11 | Referenz, 2 Bit je Base: Zahl der Wörter, dann die Wörter | u64, dann u64 je 32 Basen |
| 12 | Zählblöcke: Zahl der Wörter, dann je Block 4 Zählstände und `b`/32 Wörter BWT | u64, dann u64 |
| 13 | Zahl der Stichproben, dann je Stichprobe die Textposition | u64, dann 5 Bytes |
| 14 | CRC-32 über alle Bytes davor | u32 |

Packung: Base i eines Wortes steht in den Bits 2·(i mod 32) und 2·(i mod 32) + 1; A = 0, C = 1,
G = 2, T = 3. Die Zahl der Zählblöcke ist Zeilen / b + 1 (abgerundet), die der Stichproben
Zeilen / s (aufgerundet).

## Prüfung beim Lesen
Kennung, Formatversion, zulässige Abstände, Zeilenzahl, Größe jeder Tabelle und die CRC-32 werden
geprüft. Eine abgeschnittene oder beschädigte Datei wird abgelehnt.

## Änderungen
- Version 1 (nie veröffentlicht): Stichproben mit 32 Bit.
- Version 2: Stichproben mit 40 Bit, damit das ganze Humangenom passt.
