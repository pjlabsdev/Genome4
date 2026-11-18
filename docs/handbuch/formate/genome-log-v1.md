# Lauf-Protokoll `genome.log.v1`

Stand: gebaut (`gc-log`).

Eine Textdatei, UTF-8, ein JSON-Objekt je Zeile, Zeilenende `\n`. Jede Zeile wird sofort
geschrieben.

| Feld | Typ | Bedeutung |
|---|---|---|
| `schema` | Text | immer `genome.log.v1` |
| `time` | Text | Zeitpunkt in UTC, RFC 3339 mit Millisekunden |
| `t_ms` | Zahl | Millisekunden seit Beginn des Laufs (monoton) |
| `run` | Text | Lauf-ID, gleich in allen Zeilen eines Laufs und im Ergebnis |
| `tool` | Text | Werkzeug, z. B. `index.build` |
| `version` | Text | Programmversion |
| `level` | Text | `error`, `warn`, `info`, `debug`, `trace` |
| `code` | Text | feste Kennung, siehe [Ereignisse und Fehlercodes](../ereignisse-und-fehlercodes.md) |
| `msg` | Text | lesbarer Text; bei Stufen der Name der Stufe |
| `fields` | Objekt | Zahlen, Pfade und Namen zum Ereignis; fehlt, wenn leer |

Fehler tragen in `fields` immer `wo` (Ort), `ursache` und `was_tun`. Entscheidungen tragen `wahl`
und `grund`.

Beispielzeile:

```json
{"code":"STUFE.ENDE","fields":{"sekunden":21.75,"zeilen":128888335},"level":"info","msg":"index.suffix_array","run":"1a1eef2071a-af0","schema":"genome.log.v1","t_ms":30576,"time":"2026-11-18T13:54:03.412Z","tool":"index.build","version":"0.1.0"}
```

Ein vollständiges Beispiel liegt unter `bench/m2/ergebnisse/chr20-bau-kompakt.log.jsonl`.

## Regeln
- Keine Sequenzen, keine Genotypen.
- Ein Lauf ohne `LAUF.FERTIG` und ohne `LAUF.FEHLGESCHLAGEN` wurde von außen beendet.
- Neue Felder können hinzukommen; vorhandene ändern ihre Bedeutung nicht.
