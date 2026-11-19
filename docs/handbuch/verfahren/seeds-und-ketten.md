# Seeds sammeln, orten und zu Ketten verbinden

Stand: gebaut und geprüft (`gc-align`, `seed.rs` und `chain.rs`). Messung:
[M3 Ketten](../../bench/M3-ketten.md). Noch nicht gebaut: Ketten filtern, Alignment
(Smith-Waterman), Paarung, MAPQ.

## Seeds sammeln: drei Runden
Wie in BWA-MEM (`mem_collect_intv`), mit denselben Vorgaben:

1. **Längste exakte Treffer** ab 19 Basen, siehe [Seed-Suche](seed-suche.md).
2. **Nachsuche:** Ein langer, seltener Treffer kann kürzere, häufigere Treffer verdecken, die zum
   richtigen Ort führen. Jeder Treffer ab 28 Basen (19 × 1,5) mit höchstens 10 Fundstellen wird
   deshalb in seiner Mitte neu angesetzt; zugelassen sind dabei nur Treffer, die häufiger
   vorkommen als er selbst.
3. **Seltene Seeds von links nach rechts:** ab jeder Startposition nach rechts verlängern, bis das
   Muster länger als 19 Basen ist und seltener als 20-mal vorkommt; dann hinter dem Treffer
   weitermachen. Diese Seeds sind nicht nach links verlängert.

Die Treffer aller Runden werden nach Anfang und Ende im Read geordnet. Doppelte bleiben stehen, wie in bwa: In Wiederholungen entscheiden sie mit, welche Ketten entstehen ([Messung](../../bench/M5-echte-reads.md)).

## Orten
Jede Zeile eines Treffers wird über die Suffix-Array-Stichproben in eine Textposition
umgerechnet. Kommt ein Treffer öfter als 500-mal vor, wird gleichmäßig ausgedünnt: nur jede
(Anzahl / 500)-te Zeile, höchstens 500. Fundstellen über eine Contig- oder Stranggrenze oder in
einer gefüllten N-Strecke fallen weg. (bwa verwirft Fundstellen über Contig-Grenzen ebenfalls,
kennt aber keinen Filter für N-Strecken; dort können sich die Seeds unterscheiden.)

## Verketten
Eine Kette ist eine Folge von Seeds auf demselben Strang desselben Contigs, die im Read und in
der Referenz in gleicher Reihenfolge und ungefähr auf derselben Diagonale liegen
(`mem_chain`, `test_and_merge`). Jeder Seed wird mit der Kette verglichen, deren erster Seed in
der Referenz am nächsten vor ihm beginnt:

- Liegt er ganz innerhalb der Kette (im Read und in der Referenz), gehört er dazu und wird nicht
  eigens aufgenommen.
- Sonst wird er angehängt, wenn er in der Referenz nicht vor dem letzten Seed liegt, die
  Abstände zum letzten Seed im Read und in der Referenz sich um höchstens 100 Basen
  unterscheiden (Bandbreite) und beide Lücken kleiner als 10.000 Basen sind.
- Andernfalls beginnt er eine neue Kette. Über die Stranggrenze wird nie verkettet.

**Gewicht** einer Kette (`mem_chain_weight`): die von den Seeds überdeckten Basen, Überlappungen
einfach gezählt, im Read und in der Referenz getrennt gerechnet; das kleinere von beiden.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Dritte Runde gegen eine Vergleichssuche durch Abzählen im Text (200 Reads, wechselnde Mindestlänge und Häufigkeitsgrenze) | gleich |
| Sammeln und Orten an einem Genom mit drei Kopien eines Stücks: Nachsuche findet die Kopien, jeder Seed gleicht dem Text, Ausdünnen begrenzt die Zahl | bestanden |
| Verketten an gesetzten Fällen: eine Diagonale, Überlappung, zu weit weg, anderer Strang, anderer Contig, enthaltener Seed, rückläufige Referenz | bestanden |
| chr20, 200.000 erzeugte Reads: Anteil mit einer Kette am wahren Ort | 199.909 (99,95 %), siehe Messung |
| Die 74 Reads ohne Seed am wahren Ort einzeln erklärt (`bench/m3/explain_missing.py`) | 40 verdeckt durch längeren Treffer anderswo, 34 durch Ausdünnen, 0 unerklärt |

## Grenzen
- Die Übereinstimmung mit bwa ist hier noch nicht gemessen; Seeds und Ketten sind in bwa nicht
  von außen sichtbar. Der Vergleich folgt an den fertigen Alignments.
- Die Reihenfolge, in der Seeds an Ketten angehängt werden, folgt bwa nach bestem Verständnis des
  veröffentlichten Quelltexts; Abweichungen im Einzelfall sind möglich, bis der Vergleich steht.
- Erzeugte Reads mit reinen Austauschfehlern sind leichter als echte.

## Literatur
Li 2013, siehe [Literatur](../literatur.md).
