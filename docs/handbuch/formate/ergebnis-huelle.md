# Ergebnis-Hülle `genome.<werkzeug>.vN`

Stand: gebaut (`gc-schema`).

Jedes Werkzeug gibt sein Ergebnis als JSON in derselben Hülle aus:

```json
{
  "schema": "genome.index-build.v1",
  "tool": "genome",
  "version": "0.1.0",
  "result": { }
}
```

| Feld | Bedeutung |
|---|---|
| `schema` | Name und Version des Ergebnisformats. Die Version steigt, wenn ein Feld entfällt oder seine Bedeutung ändert; neue Felder erhöhen sie nicht |
| `tool` | Programm |
| `version` | Programmversion |
| `result` | Inhalt nach dem Schema, beschrieben auf der Seite des Werkzeugs |

## Vorhandene Schemata

| Schema | Werkzeug |
|---|---|
| `genome.fastq-check.v1` | [fastq-check](../werkzeuge/fastq-check.md) |
| `genome.bam-roundtrip.v1` | [bam-roundtrip](../werkzeuge/bam-roundtrip.md) |
| `genome.index-build.v1` | [index build](../werkzeuge/index.md) |
| `genome.index-locate.v1` | [index locate](../werkzeuge/index.md) |

Geplant: Commit-Kennung des Programms und Prüfsummen der Datenpakete in der Hülle (Manifest).
