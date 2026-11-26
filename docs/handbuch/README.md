# Handbuch Genome 4.0

Dieses Handbuch beschreibt Genome 4.0 so, dass jedes Ergebnis wissenschaftlich nachvollziehbar ist:
was ein Werkzeug tut, nach welchem Verfahren, wie es geprüft wurde, wo seine Grenzen liegen und
woher die Methode stammt. Es ist modular: eine Datei je Werkzeug, je Verfahren und je Dateiformat.

## Stand der Seiten

Jede Seite trägt oben einen Stand:

- **gebaut und geprüft**: im Programm vorhanden, mit Messung belegt (Link zur Messung auf der Seite).
- **gebaut**: im Programm vorhanden, Messung steht aus.
- **geplant**: noch nicht gebaut. Die Seite beschreibt nur, was der Plan festlegt.

Auf keiner Seite steht eine Genauigkeits- oder Tempoaussage ohne Messung.

## Inhalt

### Grundlagen
- [Überblick](ueberblick.md): was Genome 4.0 ist, Grundsätze, Stand aller Bausteine
- [Bauen und Starten](bauen-und-starten.md)
- [Architektur](architektur.md): Crates, Schichten, Regeln für macOS, Linux und Windows
- [Validierung](validierung.md): wie eine Aussage belegt wird
- [Protokollierung](protokollierung.md): Lauf-Protokoll lesen, Fehler finden
- [Ereignisse und Fehlercodes](ereignisse-und-fehlercodes.md): jede Kennung, die im Protokoll vorkommt
- [Berichte](berichte.md): ein Bericht je Werkzeug, daraus der Gesamtbericht
- [Glossar](glossar.md)
- [Literatur](literatur.md)

### Werkzeuge
- [Übersicht aller Werkzeuge](werkzeuge/README.md)
- [`genome fastq-check`](werkzeuge/fastq-check.md)
- [`genome bam-roundtrip`](werkzeuge/bam-roundtrip.md)
- [`genome index`](werkzeuge/index.md)
- [`genome sort`](werkzeuge/sort.md)
- [`genome markdup`](werkzeuge/markdup.md)
- [`genome bam-index`](werkzeuge/bam-index.md)
- [`genome flagstat`](werkzeuge/flagstat.md)
- [`genome fq2bam`](werkzeuge/fq2bam.md)
- [Vorlage für neue Werkzeugseiten](werkzeuge/vorlage.md)

### Verfahren
- [FASTQ-Prüfung und Herkunft der Lesungen](verfahren/fastq-pruefung.md)
- [BAM lesen und schreiben](verfahren/bam-schreiben.md)
- [Ausrichten auf mehreren Kernen](verfahren/mehrere-kerne.md)
- [Suffix-Array nach SA-IS](verfahren/suffix-array-sais.md)
- [FM-Index](verfahren/fm-index.md)
- [Blockweiser Index-Bau für das ganze Genom](verfahren/blockweiser-bau.md)
- [Seed-Suche: längste exakte Treffer (SMEM)](verfahren/seed-suche.md)
- [Seeds sammeln, orten und zu Ketten verbinden](verfahren/seeds-und-ketten.md)
- [Seed verlängern: Smith-Waterman mit Band](verfahren/seed-verlaengern.md)
- [Ketten filtern und zu Alignment-Bereichen verlängern](verfahren/ketten-zu-bereichen.md)
- [Bereiche bereinigen, primär bestimmen, MAPQ](verfahren/bereiche-bewerten.md)
- [Globales Alignment, CIGAR und Zusammenfügen von Bereichen](verfahren/globales-alignment.md)
- [Reihenfolge bei Gleichstand: warum Genome 4.0 wie bwa sortiert](verfahren/reihenfolge-wie-bwa.md)
- [Ausgabe als SAM (einzelne Reads)](verfahren/sam-ausgabe.md)
- [Read-Paare: Insertgröße, Partnersuche, Paarbewertung](verfahren/read-paare.md)
- [Schnellweg für leichte Reads](verfahren/schnellweg.md)
- [Grafikeinheit für die Partnersuche](verfahren/grafikeinheit.md)
- [Sortieren nach Ort mit festem Speicher](verfahren/sortieren.md)
- [Duplikate markieren](verfahren/duplikate-markieren.md)
- [BAM-Index](verfahren/bam-index.md)
- [Überblick über eine BAM: Zahlen und Abdeckung in einem Durchlauf](verfahren/bam-ueberblick.md)
- [Panel-Alignment: Reads gegen viele ähnliche Sequenzen einer Genfamilie](verfahren/panel-alignment.md)

### Formate
- [Ergebnis-Hülle `genome.<werkzeug>.vN`](formate/ergebnis-huelle.md)
- [Lauf-Protokoll `genome.log.v1`](formate/genome-log-v1.md)
- [Indexdatei `*.gcfm`](formate/gcfm.md)
- [Sortier-Teil `*.gcsort`](formate/gcsort.md)
- [Überblick über eine BAM `genome.bamqc.v1`](formate/genome-bamqc-v1.md)

### Messungen
- [Verzeichnis der Messungen](messungen.md)

### Plan
- [Genome-4.0-Plan](../plan/Genome-4.0-Plan.md)

## Regeln für dieses Handbuch

1. Jede Änderung am Programm ändert im selben Pull Request die betroffene Seite.
2. Jedes neue Werkzeug bekommt eine Seite nach der [Vorlage](werkzeuge/vorlage.md), jedes neue
   Verfahren eine Seite unter `verfahren/`, jedes neue Dateiformat eine unter `formate/`.
3. Jede Kennung im Protokoll steht in [Ereignisse und Fehlercodes](ereignisse-und-fehlercodes.md).
4. Zahlen kommen nur aus Messungen, die sich mit den genannten Befehlen wiederholen lassen.
5. Ein Test (`crates/genome/tests/handbuch.rs`) prüft bei jedem Lauf: alle Verweise im Handbuch
   zeigen auf vorhandene Dateien, jeder Unterbefehl hat eine Werkzeugseite, jede Kennung aus dem
   Quelltext ist dokumentiert.
