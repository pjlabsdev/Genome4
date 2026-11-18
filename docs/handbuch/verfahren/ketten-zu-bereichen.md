# Ketten filtern und zu Alignment-Bereichen verlängern

Stand: gebaut und geprüft (`gc-align`, `chain.rs` `filter_chains`, `region.rs`). Messung:
[M3 Bereiche](../../bench/M3-bereiche.md). Noch nicht gebaut: Bereiche zusammenführen und
bewerten (primär, sekundär, MAPQ), CIGAR, Paarung, Ausgabe.

## Ketten filtern
Wie `mem_chain_flt` in bwa. Die Ketten werden nach [Gewicht](seeds-und-ketten.md) geordnet, die
schwerste zuerst. Eine Kette entfällt, wenn

- sie im Read mit einer bereits behaltenen Kette um mindestens die Hälfte der kürzeren von beiden
  überlappt, und
- ihr Gewicht unter der Hälfte des Gewichts jener Kette liegt, und
- der Abstand der Gewichte mindestens 38 beträgt (zweimal die Seed-Mindestlänge).

Je behaltener Kette bleibt die erste, die sie stark überlappt, in jedem Fall erhalten. Sie wird
später gebraucht, um zu beurteilen, wie eindeutig der Ort ist.

## Von der Kette zum Bereich
Wie `mem_chain2aln` in bwa. Für jede Kette:

1. **Referenzfenster:** das größte Stück Referenz, das eine Verlängerung von irgendeinem Seed
   aus erreichen kann (Read-Rest plus die längste lohnende Lücke), begrenzt auf den Contig und
   den Strang der Kette. Ein Bereich reicht deshalb nie über eine Contig-Grenze.
2. **Seeds der Länge nach, längster zuerst.** Liegt ein Seed schon in einem früheren Bereich
   derselben Kette und nahe an dessen Diagonale, wird er übersprungen, außer ein bereits
   verlängerter, etwa gleich langer Seed überlappt ihn auf einer anderen Diagonale.
3. **Nach links verlängern** mit dem [Rechenkern](seed-verlaengern.md), Ankerscore ist die
   Seed-Länge. Lag das Maximum nahe am Rand des Bandes (mehr als drei Viertel der Bandbreite von
   der Diagonale), wird einmal mit doppeltem Band wiederholt.
4. **Read-Ende:** Ist der Score bis zum Ende des Reads höchstens 5 Punkte schlechter als der
   beste lokale Score, wird der Read bis zum Ende ausgerichtet; sonst wird das Ende
   abgeschnitten (Soft-Clip in der späteren Ausgabe).
5. **Nach rechts verlängern**, ebenso, Ankerscore ist der Score nach Schritt 3.

Ein Bereich enthält Anfang und Ende im Read und im Indextext, den besten lokalen Score, den
Score des gemeldeten Bereichs (`truesc`, kleiner, wenn bis zum Read-Ende verlängert wurde), die
verwendete Bandbreite und die Seed-Basen, die er überdeckt.

## Validierung

| Prüfung | Ergebnis |
|---|---|
| Kettenfilter an gesetzten Fällen (leichte überlappende Kette entfällt, erster Überlapper bleibt, Kette daneben bleibt) | bestanden |
| Bereiche mit bekanntem Ergebnis: exakter Read 100, ein Fehltreffer 95, Deletion von 3 Basen 88, Insertion von 2 Basen 92, Gegenstrang, jeweils mit genauen Grenzen | bestanden |
| Read-Enden: Fehltreffer an der vorletzten Base bleibt ausgerichtet (lokal 98, gemeldet 95), zehn fremde Basen am Anfang oder Ende werden abgeschnitten | bestanden |
| Read über eine Contig-Grenze: kein Bereich reicht darüber | bestanden |
| chr20, 200.000 erzeugte Reads: bester Bereich am wahren Ort | 99,67 %, siehe Messung |

## Grenzen
- Die Reihenfolge bei gleichem Gewicht ist die von bwa, siehe
  [Reihenfolge bei Gleichstand](reihenfolge-wie-bwa.md). Der Vergleich mit bwa steht in der
  [Messung](../../bench/M3-vergleich-bwa.md).
- Nicht gebaut ist `mem_flt_chained_seeds`: bwa prüft bei langen Reads (ab etwa 800 Basen) jeden
  Seed einer Kette mit einem kurzen Smith-Waterman und verwirft schwache. Bei Reads von 150
  Basen läuft dieser Schritt in bwa nicht; für lange Reads fehlt er hier noch.
- Ein Seed, der in einem Bereich einer anderen Kette liegt, wird wie in bwa nicht verlängert.
  Das kann das bessere Alignment kosten: Beim Read `chr20:63035123:+:3` meldet bwa (und jetzt
  auch Genome 4.0) ein Alignment mit einer Deletion von 36 Basen und 98 Punkten, obwohl der Read
  am wahren Ort ohne Lücke 135 Punkte hat. Ein Ansatz für eine spätere, gemessene Verbesserung.

## Literatur
Li 2013, siehe [Literatur](../literatur.md).
