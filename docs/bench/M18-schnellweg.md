# M18: Schnellweg und Vektorbefehle im lokalen Alignment

Cloud-Rechner, 4 Kerne. Reads: HG002 GIAB, Paare aus der chr20-BAM (421.564 Paare, 148 Basen).
Zahlen: `sekunden_ausrichten` aus `ALIGN.FERTIG`.

## Tempo
Vektoren heißt: lokales Alignment der Partnersuche und Seed-Verlängern mit SSE2/NEON aus
`gc-sys`. Ganzes Genom: je zwei Läufe, abwechselnd.

| Index | vorher (`GENOME_SIMD=0`, `--schnellweg aus`) | Vektoren | Vektoren und Schnellweg |
|---|---|---|---|
| chr20 (nur lokales Alignment in Vektoren) | 13,1 bis 14,6 s | 12,4 bis 12,9 s | 8,7 bis 9,1 s (1,5x) |
| GRCh38 kompakt | 30,9 / 31,4 s | 25,8 / 25,8 s (1,21x) | 21,5 / 21,1 s (1,46x) |

Lokales Alignment allein (148 gegen 1.000 Basen, Zufallsfolgen, `examples/swbench.rs`): 56 µs
mit Feldern, 22 µs mit Vektoren.

## Anteil Schnellweg
| Index | genommen | mehrere Fundstellen | sonst voller Weg |
|---|---|---|---|
| chr20 | 631.430 (75 %) | 152.359 | 59.339 |
| GRCh38 | 502.528 (60 %) | 293.725 | 46.875 |

## Gleichheit
- Vektoren: BAM Zeile für Zeile gleich (chr20 echte Reads, 844.006 Zeilen; chr20 simuliert,
  400.000 Zeilen; GRCh38 echte Reads, 844.609 Zeilen).
- Schnellweg: siehe [Schnellweg](../handbuch/verfahren/schnellweg.md#messung).
- Wahrheit (simuliert chr20, 200.000 Paare, Seed 3): mit und ohne Schnellweg dieselben Zahlen,
  MAPQ ab 20 13 von 386.865 falsch, Read für Read gleich.

## Offen
Ganzes Genom auf dem Mac (HG001/PJ), andere Geräte (MGI), Tempo auf NEON.
