# M3, Teil 3: Alignment-Bereiche

Stand: 18.11.2026. Verfahren: [Ketten zu Bereichen](../handbuch/verfahren/ketten-zu-bereichen.md),
[Seed verlängern](../handbuch/verfahren/seed-verlaengern.md).

## Aufbau
Wie [M3 Ketten](M3-ketten.md): chr20, 200.000 erzeugte Reads zu 150 Basen, 1 % Austauschfehler,
wahrer Ort im Namen, Vorgaben von BWA-MEM. Befehl:
`genome index ketten chr20.kompakt.gcfm --fastq sim150.fq --wahrheit-im-namen`. Ein Kern, ein Lauf.

## Ergebnis

| | Anzahl | Anteil |
|---|---|---|
| Ketten vor dem Filter | 4.848.788 | |
| Ketten nach dem Filter | 2.206.967 | |
| Bereiche | 2.212.596 | |
| Reads: bester Bereich am wahren Ort, eindeutig | 192.470 | 96,24 % |
| davon mit genau dem wahren Anfang | 192.464 | |
| Reads: bester Bereich am wahren Ort, gleich gut wie ein anderer Ort | 6.864 | 3,43 % |
| Reads: bester Bereich woanders | 666 | 0,33 % |
| davon Score unter dem des wahren Orts | 32 | 0,016 % |
| Zeit | 72,2 s | 2.771 Reads je Sekunde |

"Bester Bereich" ist der mit dem höchsten lokalen Score. Der Score des wahren Orts ist als
untere Schranke gerechnet: der ganze Read ohne Lücken, also 150 minus 5 je Fehler.

## Einordnung
- Bei 634 der 666 Reads mit dem besten Bereich an anderer Stelle ist der gefundene Score
  mindestens so hoch wie diese Schranke. Der Read passt dort also nicht schlechter als am wahren
  Ort; das ist die Folge von Wiederholungen und Lesefehlern, nicht ein verpasstes Alignment. Ein
  fertiges Programm meldet solche Reads mit niedriger MAPQ.
- Bei 32 Reads (0,016 %) wurde der wahre Ort verpasst, obwohl er besser gepasst hätte. Das passt
  zur Größenordnung der 74 Reads ohne Seed am wahren Ort aus [M3 Ketten](M3-ketten.md); einzeln
  zugeordnet sind die 32 nicht.
- Gleichstand heißt hier nur gleicher Score an verschiedenen Orten. MAPQ gibt es noch nicht.
- Wie bwa 0.7.19 bei denselben Reads abschneidet, ist noch nicht gemessen. Erst dieser Vergleich
  zeigt, ob 0,33 % und 0,016 % dem Verfahren entsprechen oder auf Abweichungen hinweisen.
- Tempo: nicht optimiert, ein Kern, Smith-Waterman Zelle für Zelle.

Rohdaten: `bench/m3/ergebnisse/bereiche.json`.
