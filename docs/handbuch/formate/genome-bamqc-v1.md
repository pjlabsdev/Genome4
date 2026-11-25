# Überblick über eine BAM `genome.bamqc.v1`

Stand: gebaut (`genome bam-bericht`, `genome fq2bam`).

Eine JSON-Datei `<name>.bamqc.json` neben der BAM. Aus ihr wird die HTML-Seite gesetzt; jede Zahl
der Seite steht hier. Verteilungen stehen dünn besetzt als `[[wert, anzahl], …]`.

| Feld | Inhalt |
|---|---|
| `schema`, `werkzeug`, `version` | `genome.bamqc.v1`, `genome`, Programmversion |
| `datei` | `bam`, `bytes`, `bai`, `bai_bytes` |
| `kopf` | `hd` (Felder der `@HD`-Zeile), `contigs`, `laenge` (Basen aller `@SQ`), `gruppen` (je Art: `gruppe`, `contigs`, `laenge`), `lesegruppen` und `programme` (Felder der `@RG`- und `@PG`-Zeilen), `kommentare` |
| `saetze` | `gesamt`; `bestanden` und `nicht_bestanden` mit den Feldern von [flagstat](../werkzeuge/flagstat.md); `hauptzeilen_ohne_ort` |
| `tags` | je Zusatzfeld `tag` und `saetze` (Datensätze, die es tragen), absteigend |
| `mapq` | `verteilung`, `reads`, `anteil_ab_20`, `anteil_0`, `anteil_60` (Hauptzeilen mit Ort ohne Duplikate) |
| `insertgroesse` | `kennzahlen` (`n`, `mittel`, `streuung`, `median`, `p05`, `p95` bis `max`), `ueber_max`, `max` (2.000), `verteilung` |
| `read_laenge` | `kennzahlen`, `verteilung` (Hauptzeilen) |
| `qualitaet` | `basen`, `mittel`, `anteil_q30` (Hauptzeilen) |
| `fehler` | `ausgerichtete_basen` (M, =, X), `fehlpaarungen` (NM minus I und D), `rate`, `nm_summe`, `rate_wie_samtools` (NM je M+I), `saetze_mit_nm`, `insertierte_basen`, `deletierte_basen`, `soft_clip_basen` |
| `abdeckung.autosomen` | `laenge`, `mittel`, `median`, `p05`, `p95`, `anteil_ab` (`1`, `10`, `20`, `30`), `verteilung` (Stellen je Abdeckung bis `max` = 1.000) |
| `abdeckung` | `fach_basen` (1.000.000), `schwellen`, `geschlecht` (`x_zu_autosomen`, `y_zu_autosomen`, `deutung`: `XY`, `XX`, `unklar`, `unbekannt`), `gekuerzte_stuecke` |
| `abdeckung.contigs[]` | `name`, `gruppe`, `laenge`, `reads` (Hauptzeilen mit Ort), `duplikate`, `basen`, `mittel`, `stellen_ab` (Stellen mit ≥ 1, 10, 20, 30x), `mapq_mittel`; für Chromosomen 1–22, X, Y, M zusätzlich `faecher` (Summe der Abdeckung je 1 Mio. Basen) |
| `lauf` | nur mit Protokoll von `fq2bam`, sonst `null`: `version`, `lauf`, `beginn`, `ende`, `befehl`, `parameter`, `system`, `ergebnis` (das ganze Ergebnis von `FQ2BAM.FERTIG`), `adapter`, `polyg`, `speicherplan`, `ausrichten` |

Arten der Contigs (`gruppe`): `autosom`, `geschlecht`, `mito`, `decoy` (Name endet auf
`_decoy`), `alt` (`_alt`), `hla` (`HLA-…`), `ebv`, `unplatziert` (`chrUn…`, `…_random`),
`andere`.

Neue Felder können hinzukommen; vorhandene ändern ihre Bedeutung nicht.
