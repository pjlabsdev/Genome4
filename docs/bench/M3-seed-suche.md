# M3, Teil 1: Seed-Suche (SMEM)

Stand: 18.11.2026. Verfahren: [Seed-Suche](../handbuch/verfahren/seed-suche.md).

## Aufbau
- Referenz: chr20 aus GRCh38 (64.444.167 Basen), Index in beiden Varianten, blockweise gebaut.
- Reads: 200.000 Stück zu 150 Basen, erzeugt mit `bench/m3/simulate_reads.py` (Saat 1): zufällige
  Orte ohne N, beide Stränge, 1 % Substitutionsfehler. Kein realistisches Fehlermodell, keine
  Indels; echte Reads folgen mit dem Alignment.
- Befehl: `genome index smem chr20.<variante>.gcfm --fastq sim150.fq [--pruefen]`, Mindestlänge 19.
- Rechner: Cloud-Container, Intel Xeon 2,3 GHz, ein Kern genutzt. Ein Lauf je Zeile, keine
  Wiederholung.

## Ergebnis

| Variante | Seeds | Reads ohne Seed | rohe Fundstellen | Zeit | Reads je Sekunde |
|---|---|---|---|---|---|
| kompakt | 471.822 | 0 | 4.098.313 | 13,6 s | 14.736 |
| schnell | 471.822 | 0 | 4.098.313 | 14,5 s | 13.779 |

Im Mittel 2,36 Seeds je Read. Die Zeit enthält das Lesen der FASTQ-Datei und das Laden des
Index nicht getrennt; sie ist eine erste Orientierung, kein Leistungsnachweis.

Mit `--pruefen` wurde jeder der 471.822 Treffer in beiden Varianten nachgeprüft (Bereich gleich
der einfachen Rückwärtssuche, links und rechts nicht verlängerbar, keiner in einem anderen
enthalten): **0 Abweichungen**.

## Einordnung
- Beide Varianten liefern dieselben Seeds, wie es sein muss.
- Die Variante schnell ist hier nicht schneller. Die Seed-Suche braucht nur die Zählblöcke, nicht
  die Suffix-Array-Stichproben, und der erwartete Vorteil der kleineren Zählblöcke zeigt sich in
  diesem einen Lauf nicht. Ob das an chr20 (passt weitgehend in den Cache) oder am Programm liegt,
  ist nicht untersucht.
- Hochrechnung auf einen 30x-Lauf wäre unseriös: Beim ganzen Genom ist der Index viel größer als
  der Cache, und die Suche ist noch nicht optimiert.

Rohdaten: `bench/m3/ergebnisse/`.
