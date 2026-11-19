# `genome flagstat`

Stand: gebaut und geprüft ([Messung](../../bench/M5-flagstat.md)). Erster Teil der
Qualitätsauswertung von [`fq2bam`](fq2bam.md), einzeln aufrufbar.

## Zweck
Zählt die Datensätze einer BAM nach ihren Flags: wie viele einen Ort haben, wie viele als
richtiges Paar liegen, wie viele Duplikate sind. Das ist die erste Kontrolle nach dem Ausrichten.
Ersetzt `samtools flagstat`.

## Aufruf
```bash
genome flagstat datei.bam [--kerne 0]
```

Die BAM muss nicht sortiert sein.

## Ausgaben
JSON, Schema `genome.flagstat.v1`, mit zwei gleich gebauten Gruppen: `bestanden` und
`nicht_bestanden` (Flag 512, Qualitätsprüfung des Geräts). Je Gruppe:

| Feld | Zeile bei samtools | gezählt wird |
|---|---|---|
| `gesamt` | total | alle Datensätze |
| `hauptzeilen` | primary | weder Flag 256 noch 2048 |
| `sekundaer`, `ergaenzend` | secondary, supplementary | Flag 256; sonst Flag 2048 |
| `duplikate`, `duplikate_hauptzeilen` | duplicates, primary duplicates | Flag 1024 |
| `mit_ort`, `hauptzeilen_mit_ort` | mapped, primary mapped | Flag 4 nicht gesetzt |
| `aus_paaren`, `read1`, `read2` | paired in sequencing, read1, read2 | Hauptzeilen mit Flag 1 |
| `richtiges_paar` | properly paired | Flag 2 und selbst mit Ort |
| `beide_mit_ort` | with itself and mate mapped | weder Flag 4 noch 8 |
| `nur_selbst_mit_ort` | singletons | Flag 8, selbst mit Ort |
| `partner_auf_anderem_contig` | with mate mapped to a different chr | beide mit Ort, Contigs verschieden |
| `partner_auf_anderem_contig_mapq5` | dasselbe mit mapQ >= 5 | zusätzlich MAPQ mindestens 5 |

Dazu die Anteile `anteil_mit_ort`, `anteil_hauptzeilen_mit_ort`, `anteil_richtiges_paar` und
`anteil_nur_selbst_mit_ort` als Zahl zwischen 0 und 1 (`null`, wenn der Nenner 0 ist), und
`sekunden`. Die Zeilen ab `aus_paaren` zählen nur Hauptzeilen, wie bei samtools.

Das Werkzeug liest nur und schreibt kein Lauf-Protokoll.

## Validierung
Gegen `samtools flagstat -O json` 1.24 an drei Dateien, je 32 Zahlen: chr20 mit markierten
Duplikaten (260.536 Datensätze), 7 Contigs (300.000) und eine BAM mit 400.000 Datensätzen und
zufälligen Flags, Contigs und MAPQ, in der jede Zeile und beide Gruppen besetzt sind. Keine
Abweichung.

## Grenzen
- Nur die Zählung nach Flags. Insertgrößen, Abdeckung, Fehlerraten und Basenqualitäten gehören
  zur Qualitätsauswertung und sind offen.
- Laufzeit am ganzen Genom nicht gemessen.
