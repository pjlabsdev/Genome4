# Überblick

Stand: 18.11.2026.

Genome 4.0 ist die Neufassung von Genome vollständig in Rust: alle Analysewerkzeuge vom FASTQ bis
zum Bericht und die App selbst. Es läuft auf macOS, Linux und Windows. Zur Laufzeit ruft es keine
fremden Programme auf (kein conda, Docker, Java, Python, Nextflow).

## Grundsätze

1. **Eigene Werkzeuge, belegt besser.** Ein eigenes Werkzeug ersetzt ein fremdes erst, wenn es auf
   ungesehenen Proben mit vorher eingefrorener Version mindestens gleich genau ist. Siehe
   [Validierung](validierung.md).
2. **Nachvollziehbar.** Jeder Lauf schreibt ein lückenloses Protokoll
   ([Protokollierung](protokollierung.md)); jedes Ergebnis trägt Schema, Werkzeug und Version
   ([Ergebnis-Hülle](formate/ergebnis-huelle.md)).
3. **Modular.** Ein Werkzeug je Aufgabe, jedes mit eigenem Bericht; der Gesamtbericht setzt nur
   zusammen ([Berichte](berichte.md)).
4. **Ehrlich.** Was ein Werkzeug nicht kann oder nicht entscheiden konnte, steht in der Ausgabe.
5. **Überall gleich.** Gleiches Rechenergebnis auf allen drei Systemen ([Architektur](architektur.md)).

## Stand der Bausteine

| Baustein | Aufgabe | Stand | Seite |
|---|---|---|---|
| `gc-io` FASTQ | lesen, vollständig prüfen, Paare abgleichen | gebaut | [Verfahren](verfahren/fastq-pruefung.md) |
| `gc-io` FASTA | lesen (auch gzip) | gebaut | [Werkzeug index](werkzeuge/index.md) |
| `gc-io` BAM | lesen und schreiben, BGZF über alle Kerne | gebaut und geprüft | [Verfahren](verfahren/bam-schreiben.md) |
| `gc-schema` | Ergebnis-Hülle | gebaut | [Format](formate/ergebnis-huelle.md) |
| `gc-log` | Lauf-Protokoll | gebaut | [Protokollierung](protokollierung.md) |
| `gc-index` SA-IS | Suffix-Array | gebaut und geprüft | [Verfahren](verfahren/suffix-array-sais.md) |
| `gc-index` FM-Index | Suche und Ortung | gebaut und geprüft (chr20) | [Verfahren](verfahren/fm-index.md) |
| `gc-index` blockweise | Index für das ganze Genom | gebaut, Messung siehe Seite | [Verfahren](verfahren/blockweiser-bau.md) |
| `genome fq2bam` | FASTQ zur fertigen BAM | geplant (M3 bis M6) | [Werkzeug](werkzeuge/fq2bam.md) |
| `gc-sys` | Systemnahes je Betriebssystem | geplant | [Architektur](architektur.md) |
| `gc-report` | Berichtsbaustein | geplant (M5) | [Berichte](berichte.md) |
| alle weiteren Werkzeuge | siehe Plan Abschnitt 6 | geplant | [Übersicht](werkzeuge/README.md) |

## Wo was steht

- Was gebaut werden soll und in welcher Reihenfolge: [Plan](../plan/Genome-4.0-Plan.md).
- Was gebaut ist und wie es arbeitet: dieses Handbuch.
- Womit es belegt ist: [Messungen](messungen.md).
