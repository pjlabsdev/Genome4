# Reihenfolge bei Gleichstand: warum Genome 4.0 wie bwa sortiert

Stand: gebaut und geprüft (`gc-align`, `ksort.rs`, `kbtree.rs`). Messung:
[M3 Vergleich mit bwa](../../bench/M3-vergleich-bwa.md).

## Das Problem
An mehreren Stellen hängt das Ergebnis von BWA-MEM davon ab, in welcher Reihenfolge gleichrangige
Dinge stehen. Das Verfahren selbst sagt dazu nichts; es ergibt sich aus den Hilfsmitteln, die bwa
benutzt. Wer bwa nachbaut und anders sortiert, bekommt bei einem kleinen Teil der Reads andere
Ketten, andere zweitbeste Bereiche und damit eine andere MAPQ. Auf chr20 waren das 63 von
200.000 Reads.

Genome 4.0 bildet beides nach. Das Ziel ist nicht, dass diese Reihenfolge richtig wäre, sondern
dass der Nachbau zuerst nachweislich dasselbe liefert wie bwa. Erst dann lässt sich jede spätere
Änderung als Verbesserung oder Verschlechterung messen.

## Sortieren (`ksort.rs`)
bwa sortiert mit `ks_introsort` aus klib: Quicksort mit Median aus drei, Teilstücke bis 17
Einträge bleiben liegen und werden am Ende gemeinsam durch Einfügen sortiert. Das Verfahren ist
nicht stabil. Betroffen:

- **Kettenfilter:** Ketten nach fallendem Gewicht. Von den Ketten, die eine behaltene stark
  überlappen, bleibt "die erste" erhalten; bei gleichem Gewicht bestimmt die Sortierung, welche
  das ist, und damit den zweitbesten Bereich für die MAPQ.
- **Bereiche bereinigen:** nach Ende in der Referenz, danach nach Score und Lage.
- **Primär bestimmen:** nach Score und Streuwert.

## Suchbaum beim Verketten (`kbtree.rs`)
Beim Verketten wird für jeden Seed die Kette gesucht, deren Anfang am nächsten vor ihm liegt.
bwa hält die Ketten dafür in einem B-Baum (klib `kbtree.h`, Knoten zu 512 Byte, höchstens 9
Ketten je Knoten). Viele Ketten haben denselben Anfang, weil an einer Stelle Seeds verschiedener
Länge liegen. Welche davon die Suche liefert, hängt vom Aufbau des Baums ab: Ein gleicher
Schlüssel in einem inneren Knoten wird vor denen in den Blättern gefunden, und ein weiterer
gleicher Schlüssel wird hinter dem ersten eingefügt, nicht am Ende. Genome 4.0 baut den Baum mit
denselben Regeln für Suchen, Einfügen und Teilen.

## Früheres Verlängern zählt für alle Ketten
Ein Seed wird nicht verlängert, wenn er schon in einem Bereich liegt und nahe an dessen
Diagonale. bwa prüft das gegen alle bisherigen Bereiche des Reads, auch die aus anderen Ketten;
die erste Fassung von Genome 4.0 prüfte nur die eigene Kette.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Sortieren: 0 bis 5.000 Einträge, Ergebnis gleich einer Bibliothekssortierung bei verschiedenen Schlüsseln; bei gleichen Schlüsseln geordnet und vollständig | bestanden |
| Sortieren: Reihenfolge dreier gleicher Schlüssel wie von Hand am Verfahren nachvollzogen (a, c, b) | bestanden |
| Baum: 3.000 Einträge mit vielen gleichen Schlüsseln, Suche liefert stets den größten Schlüssel nicht über dem gesuchten; Durchlauf geordnet und vollständig | bestanden |
| Baum: Reihenfolge gleicher Schlüssel in einem Blatt wie in klib (a, c, b; Suche liefert a) | bestanden |
| Wirkung: MAPQ-Abweichungen zu bwa auf chr20 von 66 auf 3 | siehe Messung |

## Grenzen
- Geprüft ist die Wirkung am Vergleich mit bwa, nicht der Baum Knoten für Knoten gegen klib.
- Die Knotengröße folgt aus der Größe einer Kette in bwa auf 64-Bit-Systemen (40 Byte). Ein bwa
  auf einem 32-Bit-System hätte einen anderen Baum.

## Quellen
klib (`ksort.h`, `kbtree.h`) von Attractive Chaos, MIT-Lizenz; als Teil von bwa (Li 2013), siehe
[Literatur](../literatur.md).
