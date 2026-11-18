# Validierung: wie eine Aussage belegt wird

Stand: Methode festgelegt (Plan Abschnitt 2 und 5), angewandt in M1 und M2.

## Drei Ebenen

1. **Tests im Quelltext** (`cargo test`): prüfen jede Funktion gegen eine unabhängige, einfache
   Umsetzung (zum Beispiel naive Suche statt FM-Index) und gegen Grenzfälle. Sie laufen in der CI
   auf allen Systemen.
2. **Selbstprüfung im Lauf**: Werkzeuge prüfen ihr eigenes Ergebnis, wo das möglich ist. Beispiel:
   `genome index build --vollpruefung` liest die fertige BWT rückwärts und vergleicht sie Base für
   Base mit der Referenz.
3. **Messung gegen fremde Wahrheit**: Vergleich mit dem bisherigen Standardwerkzeug und, wo es
   sie gibt, mit einer unabhängigen Wahrheit (GIAB, Assemblies, Laborergebnisse). Abgelegt unter
   [Messungen](messungen.md) mit Befehlen zum Nachrechnen.

## Regeln für Messungen

- **Version einfrieren**, dann erst auf ungesehene Proben. Fehler fließen nur in die nächste Version.
- **Entwicklungs- und Testproben getrennt**, als Listen im Repository.
- **Gleiche Daten, gleiche Hardware, gleiche Kernzahl** für das Original und das eigene Werkzeug;
  Versionen des Originals werden festgehalten.
- **Unabhängige Prüfer**: Vergleichswerkzeuge (samtools, bwa, hap.py, Truvari) bleiben fremd,
  weil ein Prüfer nicht vom Geprüften abhängen soll. Sie laufen nur in der Messung, nie im Programm.
- **Jede Abweichung wird erklärt** oder als offen ausgewiesen. Beispiel: die 65 abweichenden
  Trefferzahlen gegenüber bwa in der [M2-Messung](../bench/M2-fm-index.md).
- **Kein fairer Vergleich, keine Tempoaussage.** Zeiten aus ungleichen Bedingungen werden
  ausdrücklich so gekennzeichnet.

## Das Ablösetor

Ein eigenes Werkzeug ersetzt ein fremdes erst, wenn alle Punkte aus Plan Abschnitt 2 belegt sind:
Genauigkeit (nicht unterlegen, auf einer Hauptkennzahl besser), ungesehene Proben, Tempo und
Speicher, Robustheit, ehrliche Ausgabe, Reproduzierbarkeit, bestandener Schattenbetrieb.

## Was eine Messung enthalten muss

1. Frage und Ergebnis in einer Tabelle.
2. Daten (Quelle, Version, Umfang) und Rechner.
3. Versionen aller beteiligten Programme.
4. Befehle zum Nachrechnen und die Rohdaten im Repository (`bench/<meilenstein>/ergebnisse/`).
5. Grenzen: was nicht geprüft wurde.
