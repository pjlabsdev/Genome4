# Sortier-Teil `*.gcsort`

Stand: gebaut (`gc-bam`, `sort.rs`), Kennung `GCSORT1`.

Zwischendatei von [`genome sort`](../werkzeuge/sort.md). Sie lebt nur während eines Laufs im
Unterordner `genome-sort-<Prozess>-<Zahl>` des Temp-Ordners und ist kein Austauschformat.

Alle Zahlen little-endian.

| Nr. | Inhalt | Typ |
|---|---|---|
| 1 | Kennung `GCSORT1` und Zeilenende | 8 Bytes |
| 2 | Datensätze in sortierter Folge, jeder wie in der BAM (Längenfeld `u32`, dann der Satz), zusammen als ein Deflate-Strom (schnelle Stufe) | variabel |
| 3 | Zahl der Datensätze | u64 |
| 4 | Bytes der Datensätze vor dem Packen | u64 |
| 5 | CRC-32 über diese Bytes | u32 |

Beim Lesen werden 3 bis 5 gegen das Gelesene geprüft. Verfahren:
[Sortieren nach Ort](../verfahren/sortieren.md).
