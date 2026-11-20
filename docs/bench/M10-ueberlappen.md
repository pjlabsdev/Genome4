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

## Adaptersuche auf mehreren Kernen

Die Suche nach bekannten Adaptern rechnet blockweise (8.192 Reads) auf allen Kernen vor, ob ein
Adapter einen Read trifft, und wendet danach die Regeln von fastp der Reihe nach an; die Keime der
zweiten Suche laufen gleichzeitig. Gleiches Ergebnis (Tests, und `genome trim` schreibt für HG002
chr20 dieselben Dateien Byte für Byte). Vorbereitung für beide Seiten: 2,7 s statt 5,2 s auf 4 Kernen.

## Verworfen: Zählblöcke auf Speicherzeilen ausgerichtet

Die Seed-Suche braucht rund 63 % der Rechenzeit beim Ausrichten (Verlängern 23 %, Orte 9 %,
Ketten 4 %; HG002 chr20, Index standard). Ein Zählblock des Index hat 64 Bytes, liegt aber nicht
auf einer Grenze von 64 Bytes und belegt daher oft zwei Speicherzeilen. Versuch: Tabelle an
64 Bytes ausgerichtet. Am ganzen GRCh38-Index (standard), 4 Kerne, abwechselnd gemessen:
Ausrichten vorher 19,7 und 20,0 s, nachher 19,8 und 20,2 s, gleiche BAM. Kein Gewinn, daher nicht
übernommen. Vermutlich holt der Prozessor die Nachbarzeile ohnehin mit.

## Seed-Suche verschränkt, Zählen vom näheren Blockende (Zweig m11-seed-verschraenkt)

Nur die Seed-Suche, 200.000 Reads von HG002 chr20, ganzer GRCh38-Index (standard), ein Kern,
je drei Läufe: vorher 6,37 bis 6,83 s, nachher Read für Read 6,19 bis 6,48 s, verschränkt
(32 Reads) 5,99 bis 6,11 s; gleiche Treffer. Ganzes `fq2bam`, 4 Kerne, abwechselnd: vorher 19,98
und 19,73 s, nachher 19,98 und 19,94 s Ausrichten; gleiche BAM. Der Gewinn der Seed-Suche geht im
Rauschen der Cloud unter. Profil (callgrind, chr20): rund die Hälfte der Befehle der Seed-Suche
steckt in `occ4`, sie ist eher durch Rechnen als durch Speicherzugriffe begrenzt.
