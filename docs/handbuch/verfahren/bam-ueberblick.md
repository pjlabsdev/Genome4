# Überblick über eine BAM

Stand: gebaut (`gc-bam`, Modul `qc`).

Alle Zahlen des Überblicks entstehen in einem Durchlauf über die rohen Datensätze einer nach Ort
sortierten BAM. In `fq2bam` ist das der Durchlauf, der ohnehin den Index schreibt; es kommt kein
weiteres Lesen der BAM hinzu.

## Je Datensatz
- Zählung nach Flags wie [flagstat](../werkzeuge/flagstat.md).
- Zusatzfelder: Jedes Tag wird gezählt (Feld mit 65.536 Plätzen, Index = die zwei Buchstaben).
  Dabei wird `NM` gelesen.
- Hauptzeilen: Read-Länge, Summe der Basenqualitäten und Basen mit Qualität ≥ 30.
- Hauptzeilen mit Ort ohne Duplikate: MAPQ, mittlere MAPQ je Contig; aus Paaren mit beiden
  Partnern auf demselben Contig und TLEN > 0 die Insertgröße (so zählt jedes Paar einmal).
- Sätze mit Ort, ohne Nebenzeilen, Duplikate und nicht bestandene Prüfung (ergänzende Alignments
  zählen mit, wie bei `samtools coverage`): ausgerichtete Basen, eingefügte, fehlende und weich
  abgeschnittene Basen, Fehlpaarungen = NM minus eingefügte und fehlende Basen, und die
  Abdeckung.

## Abdeckung je Stelle ohne Feld über das Genom
Ein Feld mit einer Zahl je Stelle bräuchte für chr1 allein rund 1 GB. Stattdessen trägt jeder
Datensatz für jedes ausgerichtete Stück (M, =, X) am Anfang +1 und am Ende −1 in einen Ring von
4 Mio. Plätzen ein. Weil die Datensätze nach ihrem Anfang sortiert kommen, ändert sich an
Stellen vor dem Anfang des aktuellen Datensatzes nichts mehr: Sie werden ausgezählt, indem die
Einträge aufsummiert werden. Jede fertige Stelle geht mit ihrer Abdeckung in

- die Zähler des Contigs für ≥ 1, 10, 20, 30x,
- die Summe ihres Fachs von 1 Mio. Basen (für die Grafik entlang des Genoms),
- auf den Autosomen in die Verteilung (Abdeckung bis 1.000, darüber im letzten Fach).

Strecken ohne offene Einträge (Abdeckung 0, etwa N-Strecken oder Contigs ohne Reads) werden in
einem Schritt gezählt. Kommt ein Datensatz mit kleinerer Contig-Nummer oder kleinerem Anfang als
schon ausgezählt, bricht der Durchlauf ab: Die BAM ist nicht sortiert.

Die mittlere Abdeckung eines Contigs ist die Summe der ausgerichteten Basen geteilt durch seine
Länge, N-Strecken eingerechnet, wie bei `samtools coverage`.

## Seite und Datei
Die Zahlen gehen in die Datei [`genome.bamqc.v1`](../formate/genome-bamqc-v1.md). Die HTML-Seite
wird nur aus dieser Datei gesetzt: Grafiken als SVG im Text, keine Bibliotheken, kein Netz. Jeder
Text steht in beiden Sprachen in der Seite; der Knopf blendet eine aus. Zahlen stehen ebenfalls
zweimal, mit deutscher und englischer Schreibweise (1.234,5 und 1,234.5).
