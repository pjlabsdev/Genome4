# Globales Alignment, CIGAR und Zusammenfügen von Bereichen

Stand: gebaut und geprüft (`gc-align`, `global.rs`, `finalize.rs`). Messung:
[M3 Vergleich mit bwa](../../bench/M3-vergleich-bwa.md). Noch nicht gebaut: das CIGAR je Read mit
abgeschnittenen Enden, `NM`/`MD` in der Ausgabe.

## Wozu
Das [Verlängern](seed-verlaengern.md) liefert nur, wo ein Alignment beginnt und endet und wie
viele Punkte es hat. Die einzelnen Schritte (Base gegen Base, Insertion, Deletion) stehen im
CIGAR. Dafür wird das Stück Read noch einmal gegen das Stück Referenz ausgerichtet, diesmal mit
festen Enden auf beiden Seiten.

## Verfahren
Wie `ksw_global2` in bwa: dynamische Programmierung nach Gotoh mit affinen Lücken, nur innerhalb
eines Bands der Breite `w` um die Hauptdiagonale. Je Zelle wird neben dem Score festgehalten,
woher er kam (2 Bit für den Score selbst, je 2 Bit für die offene Deletion und Insertion). Vom
Ende aus wird der Weg zurückverfolgt.

- Bei gleich guten Wegen gilt die Reihenfolge von bwa: Base gegen Base vor Deletion vor Insertion.
- Eine Lücke öffnet nur aus einem Basenpaar heraus. `10M3I3D10M` kann so nicht entstehen.
- Speicher: Band mal Länge der Referenz, 1 Byte je Zelle.

**Bandbreite** (`bwa_gen_cigar2`): höchstens die verlangte, höchstens so groß, wie eine Lücke
bei der halben Read-Länge noch Punkte übrig ließe, mindestens der Längenunterschied plus 3.

**Gegenstrang:** Liegt der Bereich im reversen Komplement, werden Read und Referenz vor dem
Alignment umgedreht. Lücken in Wiederholungen stehen dadurch auf dem Vorwärtsstrang so weit
links wie möglich, wie es für VCF üblich ist, und das CIGAR ist schon in Vorwärtsrichtung.

`NM` (Editierabstand) wird aus dem CIGAR gezählt: ungleiche Basen plus Basen in Lücken.

## Zusammenfügen zweier Bereiche
Wie `mem_patch_reg`. Zwei Bereiche eines Reads, die hintereinander auf fast derselben Diagonale
liegen, werden zu einem, wenn

1. beide auf demselben Strang liegen und in Read und Referenz in gleicher Reihenfolge,
2. der Versatz zwischen ihren Diagonalen klein ist: höchstens 200 Basen (400, wenn sie sich
   überlappen) und im Verhältnis zur gemeinsamen Länge unter 5 % (10 %),
3. das globale Alignment über beide mindestens 90 % des Scores erreicht, den die beiden einzeln
   anteilig erwarten lassen.

Der neue Bereich bekommt den Score des gemeinsamen Alignments. Eine Deletion von 3 Basen in einem
Read von 150 wird so zusammengefügt (141 von erwarteten 153 Punkten), eine von 30 Basen nicht
(Versatz 17 % der Länge); solche Reads bleiben zwei Bereiche und werden später als geteilter
Read gemeldet.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Gesetzte Fälle: gleich, Deletion von 3, Insertion von 4, mit erwartetem Score und CIGAR | bestanden |
| 400 zufällige Paare mit Austauschen, Insertionen, Deletionen, N: Score gleich einer unabhängigen Rechnung ohne Band; Punkte des CIGAR ergeben den Score; CIGAR verbraucht Read und Referenz ganz; auch mit engstem Band | bestanden |
| Zusammenfügen: Deletion von 3 wird ein Bereich mit 141 Punkten, Deletion von 30 bleibt zwei | bestanden |
| Vergleich mit bwa auf chr20 | siehe Messung; die erzeugten Reads haben keine Indels, das Zusammenfügen ist dort praktisch nicht gefordert |

## Grenzen
- Noch nicht gegen CIGARs von bwa verglichen. Das folgt mit der SAM-Ausgabe.
- Reads mit Indels sind noch nicht im Messsatz.

## Literatur
Gotoh 1982, Li 2013, siehe [Literatur](../literatur.md).
