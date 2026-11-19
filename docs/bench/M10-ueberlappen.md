# M10: Lesen, Trimmen und Schreiben neben dem Ausrichten

Messung 19.11.2026 in der Cloud (4 Kerne), HG002 chr20:10-12 Mb, 421.564 echte Paare, gzip, Index
chr20 kompakt. Vorher: Commit 9ece955 (main), nachher: Zweig m10-tempo.

| Lauf | Ausrichten vorher | nachher | ganzer Lauf nachher | BAM gleich |
|---|---|---|---|---|
| ohne `--trimmen` | 15,8 s | 14,6 s | 17,0 s | ja |
| mit `--trimmen` | 16,5 s | 12,9 s | 20,1 s | ja |

Vorher liefen Lesen (mit Trimmen 3,9 s, davon rund 3 s Trimmen), Ausrichten (11,5 s) und
Schreiben (0,5 s) nacheinander. Nachher wartet das Ausrichten insgesamt nur noch 0,2 s auf das
Lesen. Der Vergleich der beiden Reads eines Paars prüft 8 Basen auf einmal und hört beim ersten
Überschreiten der Grenze auf (Test `within_mismatches_matches_count`); damit braucht
`genome trim` für diese 421.564 Paare 1,9 s statt rund 3,5 s, gleiche Ausgabe Byte für Byte.

Mit mehr Kernen wächst der Gewinn, weil der Anteil, der nacheinander lief, gleich bleibt, während
das Ausrichten schneller wird (auf dem Mac mit 10 Kernen nicht gemessen).

Die Adaptersuche (rund 5 s, einmal je Lauf) läuft jetzt, während der Index lädt. Bei chr20 lädt
der Index in 0,2 s, daher ist sie hier noch sichtbar; beim ganzen Genom (7,75 GB Index) verschwindet
sie darin.
