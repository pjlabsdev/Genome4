# M3: Vergleich mit bwa (Ort, Strang, MAPQ, Score)

Stand: 18.11.2026. Verfahren: [Bereiche bewerten](../handbuch/verfahren/bereiche-bewerten.md) und
die Seiten davor. Verglichen wird je Read der primäre Bereich von Genome 4.0 mit dem primären
Alignment von bwa mem. CIGAR und Ausgabe als BAM gibt es noch nicht; sie sind nicht Teil dieses
Vergleichs.

## Aufbau
- Referenz: chr20 aus GRCh38. Reads: 200.000 erzeugte zu 150 Basen, 1 % Austauschfehler, beide
  Stränge (`bench/m3/simulate_reads.py`, Saat 1), einzeln ausgerichtet (kein Paar).
- Genome 4.0: `genome index ketten chr20.kompakt.gcfm --fastq sim150.fq --tabelle unser.tsv`.
- bwa: das in pybwa 2.3.2 eingebaute bwa mem mit Vorgaben, über `bench/m3/compare_bwa_align.py`
  (Zuordnung der Reads nach Reihenfolge).
  Das ist nicht das Programm `bwa 0.7.19` von pj's Mac; welche bwa-Version pybwa einbettet, ist
  hier nicht geprüft. Der Vergleich mit genau 0.7.19 steht aus.
- Verglichen: Contig, Position des linken Endes, Strang, MAPQ, und der Alignment-Score (bwa `AS`
  gegen unseren besten lokalen Score).

## Ergebnis

| | Reads | Anteil |
|---|---|---|
| Ort und Strang gleich | 199.993 | 99,9965 % |
| davon auch Score gleich | 199.993 | |
| davon auch MAPQ gleich | 199.990 | 99,995 % aller Reads |
| Ort gleich, MAPQ verschieden | 3 | |
| anderer Ort, bei beiden MAPQ 0 | 7 | |
| anderer Ort, sonst | 0 | |
| nur bei einem der beiden ohne Alignment | 0 | |

Verteilung unserer MAPQ: 186.390 Reads mit 60, 3.030 mit 30 bis 59, 2.547 mit 1 bis 29, 8.033
mit 0.

## Verlauf

| Stand | Ort gleich | Score gleich | MAPQ gleich |
|---|---|---|---|
| erster Vergleich | 199.992 | 199.990 | 199.926 |
| Zuordnung nach Reihenfolge statt Name (Fehler im Vergleichsskript) | 199.992 | 199.992 | 199.926 |
| Seed gilt als enthalten auch in Bereichen anderer Ketten | 199.993 | 199.993 | 199.934 |
| Sortieren wie bwa (`ks_introsort`) | 199.993 | 199.993 | 199.989 |
| Suchbaum beim Verketten wie bwa (`kbtree`) | 199.993 | 199.993 | 199.990 |

Das Zusammenfügen von Bereichen (`mem_patch_reg`) ist gebaut, ändert auf diesen Reads aber
nichts: Sie enthalten keine Indels. Die Ursachen sind auf
[Reihenfolge bei Gleichstand](../handbuch/verfahren/reihenfolge-wie-bwa.md) beschrieben.
Gefunden wurden sie, nachdem der Nachbau Zeile für Zeile gegen den Quelltext von bwa gelesen und
für einzelne Reads mit der ausführlichen Ausgabe von bwa (über pybwa, Stufe 4) verglichen wurde.

## Die verbleibenden 10 Reads
Alle zehn liegen in chr20:26,4 bis 28,6 Mb, dem Satelliten-DNA-Bereich am Zentromer.

- **7 Reads, anderer Ort, beide MAPQ 0, gleicher Score:** gleich gute Orte. Warum die Wahl hier
  anders fällt als bei bwa, ist nicht geklärt. Bei den übrigen rund 6.900 Reads mit gleich guten
  Orten fällt sie gleich aus.
- **3 Reads, gleiche Lage, andere MAPQ** (0 gegen 1, 41 gegen 48, 23 gegen 21):
  - `chr20:26596727:-:1` ist geklärt: Der zweitbeste Bereich grenzt an eine N-Strecke. bwa und
    Genome 4.0 füllen N-Strecken mit verschiedenen Zufallsbasen; bei uns passt zufällig eine Base
    mehr (135 statt 134 Punkte). Beide Programme rechnen hier mit erfundenen Basen. Das ist ein
    Ansatz für eine Verbesserung: N-Strecken beim Verlängern als N werten.
  - Die beiden anderen sind nicht geklärt. Einer hat sehr häufige Seeds (19 % des Reads), bei
    denen nur jede soundsovielte Fundstelle benutzt wird.

## Zwei Berichtigungen gegenüber der ersten Fassung
- Die "2 Reads mit anderem Score" waren ein Fehler im Vergleichsskript: 33 Read-Namen kommen im
  Testsatz doppelt vor (gleicher Ort, gleiche Fehlerzahl), und das Skript ordnete nach Namen zu.
  Jetzt ordnet es nach Reihenfolge zu. Echte Score-Abweichungen gibt es nicht.
- Der Read `chr20:63035123:+:3` (bwa meldet `57M36D93M` mit 98 Punkten, Genome 4.0 meldete den
  wahren Ort mit 135) war nicht fehlendes Zusammenfügen, wie vermutet. bwa verlängert den Seed am
  wahren Ort nicht, weil er in einem Bereich einer anderen Kette liegt. Genome 4.0 bildet das
  jetzt nach und meldet dasselbe wie bwa. Das Ergebnis ist hier schlechter als vorher, aber
  gleich bwa; die Abweichung bleibt als Ansatz für eine spätere Verbesserung vermerkt.

## Der Fund zur MAPQ
Im ersten Lauf wich die MAPQ bei 5.762 Reads ab, immer nach oben. Ursache: bwa hält einen Faktor
seiner MAPQ-Formel, gedacht als ln(50) = 3,912, in einer Ganzzahl und rechnet mit 3. Nach dem
Nachbilden waren es 66, heute 3. Einzelheiten in
[Bereiche bewerten](../handbuch/verfahren/bereiche-bewerten.md).

## Gegen die Wahrheit
Aus demselben Lauf (siehe [M3 Bereiche](M3-bereiche.md)): Der beste Bereich liegt bei 99,67 % der
Reads am wahren Ort. Da Ort und Strang bei 199.993 von 200.000 Reads mit bwa übereinstimmen,
liegt bwa auf diesen Daten praktisch gleich. Eine Auswertung "falsch bei MAPQ ab 30" für beide
Programme fehlt noch.

## Was der Vergleich nicht zeigt
- Nur chr20, nur erzeugte Reads mit Austauschfehlern, keine Indels, keine echten Reads, keine
  Paare. Indels und geteilte Reads sind genau die Fälle, in denen das Zusammenfügen und der
  CIGAR zählen; sie sind hier nicht geprüft.
- Kein Tempovergleich: bwa brauchte hier 59 s, Genome 4.0 78 s, beide auf einem Kern; Genome 4.0
  rechnet aber noch keinen CIGAR und schreibt keine Ausgabe, und nichts ist optimiert.

Rohdaten: `bench/m3/ergebnisse/vergleich_bwa.json`, `mapq.json`.
