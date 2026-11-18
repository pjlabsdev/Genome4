# Bereiche bereinigen, primär bestimmen, MAPQ

Stand: gebaut und geprüft (`gc-align`, `finalize.rs`). Messung und Vergleich mit bwa:
[M3 Vergleich mit bwa](../../bench/M3-vergleich-bwa.md). Noch nicht gebaut: Zusammenfügen zweier
Bereiche über eine große Lücke (`mem_patch_reg`), CIGAR, Paarung, Ausgabe.

## Doppelte entfernen
Wie `mem_sort_dedup_patch` in bwa, ohne das Zusammenfügen. Verschiedene Seeds derselben Stelle
führen oft zum selben Bereich. Überlappen sich zwei Bereiche in der Referenz und im Read jeweils
zu mehr als 95 % des kürzeren, bleibt der mit dem höheren Score. Danach wird nach fallendem
Score, dann nach Lage geordnet; genau gleiche Bereiche entfallen.

## Primär und untergeordnet
Wie `mem_mark_primary_se`. Die Bereiche werden nach Score geordnet. Bei gleichem Score
entscheidet ein Streuwert aus der laufenden Nummer des Reads und der Position in der Liste
(Thomas Wangs 64-Bit-Hash, wie in bwa). So hängt die Wahl zwischen gleich guten Orten nicht von
der Lage im Genom ab, bleibt aber bei jedem Lauf dieselbe.

Der beste Bereich ist primär. Jeder weitere ist dem ersten besseren untergeordnet, den er im Read
zu mindestens der Hälfte des kürzeren von beiden überlappt; sonst ist er selbst primär (ein
anderer Teil des Reads, etwa bei einem geteilten Read). Der übergeordnete Bereich merkt sich:

- `sub`: den Score des besten ihm untergeordneten Bereichs,
- `sub_n`: wie viele untergeordnete Bereiche höchstens 7 Punkte schlechter sind (das Größere aus
  Treffer plus Fehltreffer und Lücke öffnen plus verlängern).

## MAPQ
Wie `mem_approx_mapq_se`. Mit Score `s`, zweitbestem Score `sub` (mindestens 19, die
Seed-Mindestlänge), Länge `l` des Bereichs (das Größere aus Read- und Referenzlänge):

```text
identität = 1 - (l - s) / 5 / l
t         = (3 / ln l) * identität²          für l ab 50, sonst identität²
MAPQ      = 6,02 * (s - sub) * t²            gerundet
MAPQ     -= 4,343 * ln(sub_n + 1)            gerundet, wenn sub_n > 0
MAPQ      = auf 0 bis 60 begrenzt, dann mal (1 - Anteil des Reads in sehr häufigen Seeds)
```

Ist `sub` mindestens so groß wie `s`, ist die MAPQ 0.

**Die 3 im Zähler.** In bwa ist der Faktor als ln(50) = 3,912 gedacht, wird aber in einer
Ganzzahl gehalten und deshalb zu 3 abgeschnitten. Genome 4.0 bildet das bewusst nach, damit die
MAPQ mit bwa vergleichbar bleibt und nachgelagerte Werkzeuge, die auf bwa-MAPQ geeicht sind,
gleich arbeiten. Gefunden wurde das am Vergleich: Mit 3,912 lag die MAPQ bei 5.762 von 200.000
Reads über der von bwa. Ob eine eigene, besser geeichte MAPQ sinnvoll ist, ist eine spätere
Frage und bräuchte eine Messung gegen die Wahrheit.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Doppelte: gleicher, fast gleicher schlechterer und enthaltener Bereich entfallen, anderer Ort bleibt | bestanden |
| Primär und untergeordnet, `sub`, `sub_n` an gesetzten Fällen; Bereiche in verschiedenen Read-Teilen beide primär | bestanden |
| MAPQ an gesetzten Fällen, darunter zwei von bwa ausgegebene Werte (145 gegen 130 ergibt 32; 130 gegen 120 ergibt 19) | bestanden |
| Reihenfolge gleich guter Orte hängt von der Read-Nummer ab, nicht von der Lage (200 Nummern, beide Orte etwa gleich oft) | bestanden |
| Vergleich mit bwa auf chr20: Ort, Strang, MAPQ, Score je Read | siehe Messung |

## Literatur
Li 2013, siehe [Literatur](../literatur.md).
