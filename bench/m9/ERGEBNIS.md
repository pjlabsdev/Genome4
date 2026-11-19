# M9: `genome trim` gegen fastp 1.4.0, Read für Read

Lauf 19.11.2026 auf pjs Mac (Apple M5, 16 GB), Zweig `m9-trimmen` Commit bd1a40d, fastp 1.4.0
aus der Genome-3-Umgebung, Aufruf wie in Genome 3 (`--detect_adapter_for_pe --thread 8`).
Skript: `vergleich.sh`, Auswertung `trim-vergleich.py`. Je die ersten 1.000.000 Paare.

## Ergebnis

| Satz | Paare | behalten (beide) | nur bei einem | andere Länge | andere Folge |
|---|---|---|---|---|---|
| HG001, NovaSeq (`@A00744…`, PolyG an) | 1.000.000 | 984.668 | 0 | 0 | 0 |
| PJ, MGI DNB (`@E250…`, PolyG aus) | 1.000.000 | 1.000.000 | 0 | 0 | 0 |

Jeder behaltene Read ist in Folge und Qualität gleich, in R1 und R2. Auch die Zähler stimmen:

| | fastp (Reads) | genome (Paare) |
|---|---|---|
| HG001 Qualität | 29.620 | 14.810 |
| HG001 zu kurz | 932 | 466 |
| HG001 zu viele N | 106 | 53 |
| HG001 Adapter-Dimer | 6 | 3 |
| HG001 Adapterbasen | 1.054.937 | 1.054.937 |
| PJ Adapterbasen | 396.376 | 396.376 |

Beide Werkzeuge erkennen in beiden Sätzen keine Adapterfolge (fastp: `unspecified`); gekürzt wird
über die Überlappung der Paare.

## Laufzeit und Speicher (`/usr/bin/time -l`)

| Satz | fastp (8 Threads) | genome trim |
|---|---|---|
| HG001 | 40,0 s, 1,32 GB | 9,8 s, 0,46 GB |
| PJ | 37,7 s, 1,33 GB | 9,0 s, 0,49 GB |

fastp schreibt zusätzlich HTML und JSON; Eingabe für beide unkomprimiert, Ausgabe unkomprimiert.
