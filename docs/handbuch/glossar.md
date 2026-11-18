# Glossar

| Begriff | Bedeutung |
|---|---|
| Ablösetor | Bedingungen, die ein eigenes Werkzeug erfüllen muss, bevor es ein fremdes ersetzt ([Validierung](validierung.md)) |
| BGZF | blockweise gzip-Kompression der BAM; jeder Block ist für sich lesbar |
| BWT | Burrows-Wheeler-Transformation: zu jedem Suffix in sortierter Reihenfolge das Zeichen davor |
| Contig | zusammenhängende Sequenz der Referenz, z. B. ein Chromosom |
| Differenzüberdeckung | kleine Menge von Resten modulo v, die jeden Abstand als Differenz zweier Elemente enthält ([blockweiser Bau](verfahren/blockweiser-bau.md)) |
| Endzeichen | gedachtes Zeichen hinter dem Text, kleiner als jede Base |
| FM-Index | BWT mit Zähltabelle und Stichproben des Suffix-Arrays; findet exakte Treffer in einer Zeit, die nur von der Musterlänge abhängt ([Verfahren](verfahren/fm-index.md)) |
| Fundstelle, roh | jede Zeile des Suffix-Array-Bereichs eines Musters, auch über Contig-Grenzen und in N-Strecken |
| Lauf-ID | Kennung eines Laufs im Protokoll und im Ergebnis |
| LMS-Suffix | Suffix vom Typ S, dessen Vorgänger vom Typ L ist (Baustein von SA-IS) |
| N-Strecke | Folge unbekannter Basen in der Referenz; im Index mit festen Pseudozufallsbasen gefüllt und vermerkt |
| Primärzeile | Zeile der BWT, deren Suffix an Textposition 0 beginnt |
| Schattenbetrieb | altes und neues Werkzeug laufen auf derselben Probe, Abweichungen werden erklärt |
| Stichprobe (Suffix-Array) | jede k-te Zeile des Suffix-Arrays wird gespeichert, die übrigen werden über die BWT erreicht |
| Suffix-Array | Anfangspositionen aller Suffixe eines Textes in sortierter Reihenfolge |
| Treffer, gefiltert | Fundstelle, die ganz in einem Contig liegt und keine N-Strecke berührt |
| Zählblock | Abschnitt der BWT mit vorangestelltem Zählstand der vier Basen |
