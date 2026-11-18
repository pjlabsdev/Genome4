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
| davon Score über dem des wahren Orts | 538 | 0,269 % |
| davon Score gleich dem des wahren Orts | 96 | 0,048 % |
| davon Score unter dem des wahren Orts | 32 | 0,016 % |
| Zeit | 72,2 s | 2.771 Reads je Sekunde |

"Bester Bereich" ist der mit dem höchsten lokalen Score. Der Score des wahren Orts ist für die
666 Reads direkt gemessen: der beste zusammenhängende Abschnitt auf der wahren Diagonale ohne
Lücken, mit denselben Punkten (1 je Treffer, 4 Abzug je Austausch). Die erzeugten Reads enthalten
nur Austausche, deshalb ist das dort das richtige Alignment. Ein Alignment mit Lücken könnte am
wahren Ort höchstens noch mehr Punkte bringen; der Wert bleibt also streng genommen eine untere
Schranke, und die drei Zahlen sind entsprechend zu lesen: 538 und 96 sind Obergrenzen, 32 ist eine
Untergrenze. Die erste Fassung dieser Seite rechnete mit einer gröberen Schranke (ganzer Read ohne
Lücken) und schloss daraus zu viel; die Zahl 32 hat sich mit der genaueren Messung nicht geändert.

## Einordnung
- Bei 538 Reads passt der gefundene Ort nach Punkten besser als der wahre Ort ohne Lücken. Das
  ist bei Wiederholungen mit Lesefehlern zu erwarten, kein Programm kann den wahren Ort dann am
  Score erkennen. Ob jeder dieser Fälle so liegt, ist nicht einzeln geprüft.
- Bei 96 Reads sind beide Orte gleich gut, aber am wahren Ort steht kein gleich guter Bereich in
  unserer Liste. Das ist ein verpasster Gleichstand. Er verfälscht nicht den Ort (beide sind
  gleich berechtigt), kann aber die MAPQ zu hoch ausfallen lassen.
- Bei mindestens 32 Reads (0,016 %) wurde der wahre Ort verpasst, obwohl er besser gepasst hätte.
  Das passt zur Größenordnung der 74 Reads ohne Seed am wahren Ort aus [M3 Ketten](M3-ketten.md);
  einzeln zugeordnet sind die 32 nicht.
- Gleichstand heißt hier nur gleicher Score an verschiedenen Orten. MAPQ gibt es noch nicht.
- Wie bwa 0.7.19 bei denselben Reads abschneidet, ist noch nicht gemessen. Erst dieser Vergleich
  zeigt, ob 0,33 % und 0,016 % dem Verfahren entsprechen oder auf Abweichungen hinweisen.
- Tempo: nicht optimiert, ein Kern, Smith-Waterman Zelle für Zelle.

Rohdaten: `bench/m3/ergebnisse/bereiche.json`, Aufteilung der 666 Reads in `bench/m3/ergebnisse/bereiche-wahrer-ort.json`.
