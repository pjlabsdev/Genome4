# M2, Teil 2: Index für das ganze Genom (GRCh38)

Stand: Lauf vom 18.11.2026 ist noch nicht abgeschlossen. Diese Seite enthält bisher nur den
Aufbau der Messung; Zeiten, Speicherbedarf, Dateigröße und Prüfergebnis werden eingetragen,
sobald Bau und Vollprüfung durch sind. Bis dahin gilt keine Zahl als gemessen.

## Aufbau

- Referenz: GRCh38, GIAB-Analyse-Satz ohne ALT-Contigs, 3.099.922.541 Basen in 195 Contigs,
  1.037 N-Strecken. Text des Index (beide Stränge): 6.199.845.082 Zeichen.
- Befehl:
  `genome index build grch38.fa.gz grch38.kompakt.gcfm --variante kompakt --verfahren blockweise --block-mio 96 --pruefen 1000000 --vollpruefung --log-stufe debug`
- Verfahren: [blockweiser Bau](../handbuch/verfahren/blockweiser-bau.md),
  Format: [gcfm](../handbuch/formate/gcfm.md).
- Prüfungen: Stichprobe mit 1.000.000 Abfragen gegen die Referenz und Vollprüfung per
  Rückwärtslauf durch die BWT (siehe [FM-Index](../handbuch/verfahren/fm-index.md)).
- Rechner: Cloud-Container (Linux x86-64), nicht pj's Mac. Die Zeiten sind deshalb nicht mit
  der Basismessung auf dem Mac vergleichbar.

## Ergebnis

Noch offen.
