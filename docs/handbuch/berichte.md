# Berichte

Stand: geplant (Plan Abschnitt 14). Der Rahmen `gc-report` entsteht mit `fq2bam` in Meilenstein M5.
Gebaut ist der erste eigenständige Bericht: der Überblick über eine BAM
([`bam-bericht`](werkzeuge/bam-bericht.md)), den `fq2bam` neben jede BAM schreibt.

## Grundsatz

- Jedes Werkzeug (HLA, BRCA, KIR, Lp(a), Pharmakogenetik und alle anderen) erstellt einen
  eigenständigen Bericht, auch ohne die App.
- Der Gesamtbericht rechnet nichts. Er setzt die Berichte der Werkzeuge zusammen; jede Zahl ist
  zum Werkzeugbericht rückverfolgbar.

## Drei Ebenen je Werkzeug

1. Ergebnis `<werkzeug>.json` ([Ergebnis-Hülle](formate/ergebnis-huelle.md)): die Fakten.
2. Berichtsbaustein `<werkzeug>.report.json` (Schema `genome.report.v1`): die Darstellung.
3. Bericht als PDF und als Ansicht in der App, gesetzt aus dem Baustein.

## Teile eines Berichtsbausteins

Kopf (Werkzeug, Version, Probe, Datenstand, Lauf-ID), Kernaussagen, Ergebnisse, Qualität des
Laufs, Methode, Validierungsstand, Grenzen, Literatur.

Die Abschnitte Methode, Validierungsstand, Grenzen und Literatur eines Berichts kommen aus der
Handbuchseite des Werkzeugs. Handbuch und Bericht sagen deshalb dasselbe.

Einzelheiten: [Plan Abschnitt 14](../plan/Genome-4.0-Plan.md).
