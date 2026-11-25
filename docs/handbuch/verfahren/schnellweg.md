# Schnellweg für leichte Reads

Die meisten Reads eines menschlichen Genoms passen in voller Länge ohne Lücke an genau eine Stelle,
mit höchstens einer Abweichung. Für sie liefert das volle Verfahren (Seed-Suche mit Nachsuche,
Ketten, Verlängern) genau einen Bereich über den ganzen Read. Der Schnellweg findet diesen Bereich
direkt und überspringt die Seed-Suche. `fq2bam --schnellweg an` schaltet ihn ein; Vorgabe ist
`aus`. Quelltext: `crates/gc-align/src/fast.rs`.

## Ablauf
1. Den Read in Stücke von mindestens 21 Basen teilen (7 bei 148 Basen) und jedes Stück im
   FM-Index suchen, auf beiden Strängen.
2. Hat ein Stück mehr als eine Fundstelle, geht der Read den vollen Weg. Höchstens ein Stück darf
   gar keine Fundstelle haben.
3. Das erste Stück mit einer Fundstelle gibt den Anfang. Der ganze Read wird dort Base für Base
   mit der Referenz verglichen: höchstens eine Abweichung, kein N, innerhalb eines Contigs.
4. Jedes Stück ohne die Abweichung muss genau eine Fundstelle haben, das Stück mit der Abweichung
   keine. Sonst voller Weg.
5. Der Bereich bekommt die Punkte, die das volle Verfahren vergibt: über den ganzen Read, und als
   besten lokalen Wert den ohne das Ende hinter der Abweichung, wenn diese weniger als fünf Basen
   vor einem Ende liegt.

## Warum die MAPQ gleich bleibt
bwa senkt die MAPQ bei 148 Basen erst unter 60, wenn eine zweite Kopie höchstens etwa 16 Punkte
schlechter ist, also höchstens vier Ereignisse (Abweichung oder Lücke) vom Read entfernt. Eine
solche Kopie passt in mindestens drei der sieben Stücke genau. Höchstens eines davon ist das Stück
mit der Abweichung; ein anderes hat dann zwei Fundstellen, und der Read geht den vollen Weg.
Dasselbe gilt bis sechs Ereignisse, also für alle Kopien, die bwa als `XA` meldet (ab 0,8 des
besten Scores).

Eine erste Fassung mit drei Stücken war schneller, gab aber auf chr20 bei 16 Reads mit naher
zweiter Kopie MAPQ 60 statt 40 bis 57. Sie wurde verworfen.

## Was sich ändert
- `XS` bei Reads, für die das volle Verfahren nur zufällige kurze Treffer anderswo findet
  (Score unter 0,8 des besten, meist 19 bis 50). Der Schnellweg schreibt dort `XS:i:0`.
- In sehr seltenen Fällen bei Gleichstand ein anderer von zwei gleich guten Orten eines geretteten
  Partners (1 Read von 843.128 auf dem ganzen Genom, gleicher Score).

## Prüfmodus
`--schnellweg pruefen` rechnet für jeden Read beide Wege, gibt den vollen aus und zählt in
`ALIGN.FERTIG` unter `schnellweg`, wie oft der Schnellweg denselben Bereich geliefert hätte
(`gleich_wie_voller_weg`, `anders_als_voller_weg`) und warum Reads den vollen Weg gingen
(`mehrere_fundstellen`, `zu_wenig_stuecke`, `stuecke_uneinig`, `pruefung_nicht_bestanden`).

## Messung
Siehe [M18 Schnellweg](../../bench/M18-schnellweg.md).

| Prüfung | Ergebnis |
|---|---|
| HG002 chr20, 843.128 echte Reads, Index chr20: Zeilen gegen `--schnellweg aus` | MAPQ, Flags, Orte, CIGAR gleich; `XS` anders in 9 % |
| dieselben Reads, Index GRCh38: Zeilen gegen `aus` | MAPQ gleich; `XS` anders in 44 %, 1 CIGAR anders (Gleichstand) |
| 400.000 simulierte Reads chr20 gegen die Wahrheit (`genome-validate`) | jeder Read gleich richtig oder falsch wie ohne Schnellweg |
| Tests: leichte Reads gleich wie der volle Weg, zweite Kopie und N gehen den vollen Weg | bestanden |

## Grenzen
- Gemessen bisher auf chr20-Reads einer Probe (HG002, Illumina 148 Basen). Das ganze Genom und
  andere Geräte stehen aus; bis dahin bleibt die Vorgabe `aus`.
- Reads unter 60 Basen gehen immer den vollen Weg.
