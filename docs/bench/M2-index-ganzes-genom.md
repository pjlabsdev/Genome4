# M2, Teil 2: Index für das ganze Genom (GRCh38)

Stand: 18.11.2026. Der Index für das ganze Humangenom ist gebaut und vollständig geprüft,
ohne Abweichung. Rohdaten des Laufs: [`bench/m2/ergebnisse/grch38_build.json`](../../bench/m2/ergebnisse/grch38_build.json),
Protokoll [`grch38.kompakt.gcfm.log.jsonl`](../../bench/m2/ergebnisse/grch38.kompakt.gcfm.log.jsonl),
Ressourcen [`grch38_build.rusage.json`](../../bench/m2/ergebnisse/grch38_build.rusage.json).

## Aufbau

- Referenz: GRCh38, GIAB-Analyse-Satz ohne ALT-Contigs, 3.099.922.541 Basen in 195 Contigs,
  1.037 N-Strecken. Text des Index (beide Stränge): 6.199.845.082 Zeichen.
- Befehl:
  `genome index build grch38.fa.gz grch38.kompakt.gcfm --variante kompakt --verfahren blockweise --block-mio 96 --pruefen 1000000 --vollpruefung --log-stufe debug`
- Verfahren: [blockweiser Bau](../handbuch/verfahren/blockweiser-bau.md),
  Format: [gcfm](../handbuch/formate/gcfm.md).
- Prüfungen: Stichprobe mit 1.000.000 Abfragen gegen die Referenz und Vollprüfung per
  Rückwärtslauf durch die BWT (siehe [FM-Index](../handbuch/verfahren/fm-index.md)).
- Rechner: Cloud-Container (Linux x86-64), nicht pj's Mac. Die Zeiten sind deshalb nicht mit
  der Basismessung auf dem Mac vergleichbar.

## Ergebnis

Rechner: 4 Kerne (`kerne=4` laut Protokoll). Während der ersten rund 15 Minuten liefen auf
demselben Rechner zeitweise Testläufe mit; die Zeiten des Baus sind deshalb eher zu hoch als zu
niedrig.

| Stufe | Zeit | Kennzahl |
|---|---|---|
| FASTA lesen (gzip) | 16,8 s | 3.099.922.541 Basen, 195 Contigs |
| Text beider Stränge packen | 7,0 s | 6.199.845.082 Zeichen, 1,55 GB |
| Stichprobenränge (Differenzüberdeckung mod 1024, 40 Reste) | 54,0 s | 242.181.480 Plätze, 0,97 GB |
| Fächer zählen und zu Blöcken fassen | 27,9 s | 16.777.216 Fächer, 65 Blöcke zu höchstens 96 Mio. Suffixen |
| Blöcke sortieren, BWT und Stichproben schreiben | 1.001,9 s | 6.199.845.083 Zeilen; 40.905.666 Vergleiche über den Stichprobenrang entschieden |
| **Bau gesamt** | **1.107,5 s (18,5 min)** | |
| Stichprobenprüfung | 97,0 s | 1.000.000 Abfragen, 687.249 mit Treffern, 8.050.828 Fundstellen verglichen, **0 Abweichungen** |
| Vollprüfung | 2.310,0 s (38,5 min) | 6.199.845.082 Zeichen geprüft, **0 Abweichungen** |
| Index schreiben | 3,5 s | 4.843.653.489 Bytes |
| **Lauf gesamt** | **3.518,4 s (58,6 min)** | CPU-Zeit 6.160 s |

- Spitzenbedarf an RAM über den ganzen Lauf (Bau und Prüfungen): **8,68 GB**. Das liegt unter dem
  Ziel von 12 GB aus dem Plan. Welche Stufe die Spitze setzt, ist nicht getrennt gemessen.
- Dateigröße: 4,84 GB, das sind 12,5 Bits je Referenzbase (Variante kompakt), wie gerechnet.

## Aussage der Prüfungen

- **Vollprüfung:** Der Rückwärtslauf durch die BWT gibt Base für Base die Referenz beider Stränge
  wieder, und jede gespeicherte Suffix-Array-Stichprobe steht an der Stelle, die der Rückwärtslauf
  verlangt. Damit sind BWT, Zählstände und Stichproben dieses Index mit der Referenz konsistent.
- **Stichprobenprüfung:** Für 1.000.000 Abfragen stimmen alle gemeldeten Fundstellen mit der
  Referenz überein.

## Grenzen dieser Messung

- Ein Lauf, ein Rechner, keine Wiederholung; die Zeiten haben keine Streuungsangabe.
- Kein Vergleich mit bwa auf dem ganzen Genom. Der Vergleich mit bwa liegt bisher nur für chr20
  vor (siehe [`M2-fm-index.md`](M2-fm-index.md)).
- Die Byte-Gleichheit von blockweisem und direktem Bau ist nur bis chr20 gezeigt; für das ganze
  Genom ist der direkte Bau wegen des Speicherbedarfs nicht gelaufen. Hier trägt die Vollprüfung.
- Nicht auf pj's Mac gemessen (16 GB RAM, langsame externe SSD). Der Lauf dort steht aus.
- Nur die Variante kompakt; die Variante schnell (24 Bits je Base, gerechnet rund 9,3 GB) ist
  für das ganze Genom noch nicht gebaut.
