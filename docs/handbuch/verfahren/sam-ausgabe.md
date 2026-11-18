# Ausgabe als SAM (einzelne Reads)

Stand: gebaut und geprüft (`gc-align`, `sam.rs`; Prüfbefehl `genome index sam`). Messung:
[M3 SAM gegen bwa](../../bench/M3-sam.md). Noch nicht gebaut: Paare, BAM, mehrere Kerne,
Lesegruppen (`@RG`), ALT-Contigs.

## Was je Read geschrieben wird
Wie `mem_reg2sam` und `mem_aln2sam` in bwa.

- **Primäre Zeile:** der beste Bereich des Reads, sofern sein Score mindestens 30 ist (`-T`).
- **Ergänzende Zeilen (Flag 2048):** jeder weitere Bereich mit Score ab 30, der einen anderen
  Teil des Reads abdeckt (geteilter Read). Ihre MAPQ ist höchstens die der primären Zeile. Die
  abgeschnittenen Read-Enden stehen dort als harter Schnitt (`H`), Basen und Qualitäten nur für
  den eigenen Teil. Alle Zeilen eines Reads verweisen mit `SA` aufeinander.
- **Untergeordnete Bereiche** (derselbe Teil des Reads passt auch woanders) bekommen keine eigene
  Zeile. Sind sie fast so gut wie ihr übergeordneter (über 80 % des Scores) und höchstens 5,
  stehen sie in dessen `XA` als `Contig,±Position,CIGAR,NM;`. Bei mehr als 5 entfällt `XA`.
- **Ohne Ort:** Flag 4, keine Position, Basen und Qualitäten wie gelesen.

Basen und Qualitäten stehen in Richtung des Vorwärtsstrangs; bei Flag 16 sind die Basen das
reverse Komplement des Reads und die Qualitäten umgedreht.

## Angaben je Zeile

| Angabe | Bedeutung |
|---|---|
| `NM` | Editierabstand zur Referenz |
| `MD` | Referenz an den abweichenden Stellen |
| `AS` | Score des Alignments |
| `XS` | Score des besten untergeordneten Bereichs, 0 wenn es keinen gibt |
| `SA` | die anderen Zeilen desselben Reads: `Contig,Position,Strang,CIGAR,MAPQ,NM;` |
| `XA` | fast gleich gute andere Orte, siehe oben |

Kopf: je Contig eine `@SQ`-Zeile, dazu `@PG` mit Version und Aufruf.

## Zwei Feinheiten
- **Die Schwelle für `XA`** steht in bwa in einer Zahl einfacher Genauigkeit; 0,8 ist dort
  geringfügig größer als 0,8. Ein Bereich mit genau 80 % des Scores kommt deshalb nicht in `XA`.
  Genome 4.0 rechnet genauso. Gefunden am Vergleich: ohne das wich `XA` bei 322 Reads ab.
- **N-Strecken, bewusste Abweichung von bwa.** bwa füllt N-Strecken der Referenz mit
  Zufallsbasen und findet dort gelegentlich einen Seed von 19 oder 20 Basen. Der taucht dann als
  `XS:i:19` oder `XS:i:20` auf. Genome 4.0 verwirft Fundstellen, die eine N-Strecke berühren;
  `XS` ist dort 0. Das betrifft auf chr20 rund 25 von 200.000 Reads und in 3 Fällen die MAPQ
  um einen Punkt. Das Verhalten von bwa ist hier ein Zufallsprodukt, deshalb wird es nicht
  nachgebaut.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Einfacher Read auf beiden Strängen: alle Felder, Qualitäten bei Flag 16 umgedreht | bestanden |
| Read ohne Ort: Flag 4, `*`-Felder | bestanden |
| Geteilter Read: primäre und ergänzende Zeile, `90M60S` und `90H60M`, `SA` gegenseitig, MAPQ der ergänzenden nicht höher | bestanden |
| Zwei gleich gute Orte: MAPQ 0, `XS` gleich `AS`, der andere Ort in `XA` | bestanden |
| Ganze Zeilen gegen bwa, chr20, 200.000 Reads: alle Felder gleich | 199.964 ohne Indels, 199.956 mit Indels; siehe Messung |

## Grenzen
- Die Testsätze enthalten keine geteilten Reads. Ergänzende Zeilen und `SA` sind nur durch den
  eigenen Test belegt, nicht durch den Vergleich mit bwa.
- Ein Kern, kein BAM, keine Paare.

## Literatur
Li 2013; Li et al. 2009; SAM-Spezifikation, siehe [Literatur](../literatur.md).
