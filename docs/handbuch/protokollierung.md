# Protokollierung

Stand: Grundlage gebaut (`gc-log`), genutzt von `genome index build`. Weitere Teile nach Plan
Abschnitt 15 sind unten als geplant markiert.

## Was jeder Lauf schreibt

Jeder Lauf erzeugt einen Ereignisstrom in zwei Formen:

- lesbare Zeilen auf der Fehlerausgabe,
- dieselben Ereignisse als JSON-Zeilen in der Protokolldatei
  ([Format `genome.log.v1`](formate/genome-log-v1.md)). Bei `genome index build` heißt sie
  `<index>.log.jsonl`.

Jede Zeile wird sofort auf die Platte geschrieben. Nach einem Absturz oder Plattenausfall ist
alles bis zur letzten Zeile lesbar.

## Aufbau eines Laufs

| Ereignis | Inhalt |
|---|---|
| `LAUF.START` | vollständige Befehlszeile, alle Parameter einschließlich Standardwerten, System, Architektur, Kerne, Arbeitsordner |
| `LAUF.EINGABE` | jede Eingabedatei mit Rolle, Pfad, Größe, Änderungszeit |
| Entscheidungen | jede automatische Wahl mit Grund, z. B. `INDEX.VERFAHREN` |
| `STUFE.START` / `STUFE.ENDE` | jede Stufe mit Dauer und ihren Zahlen; `STUFE.ENDE` und `STUFE.LEBT` tragen `speicher_spitze_bytes`, den bis dahin größten belegten Arbeitsspeicher des Prozesses (Linux, macOS; fehlt, wo das System ihn nicht nennt). So zeigt das Protokoll, welche Stufe den Speicher wachsen lässt, etwa Laden des Index, Laden des Strobe-Index oder Ausrichten. |
| `STUFE.LEBT` | Lebenszeichen langer Stufen mit Fortschritt |
| `LAUF.FERTIG` oder `LAUF.FEHLGESCHLAGEN` | Zusammenfassung bzw. Fehlertext |

Alle Kennungen: [Ereignisse und Fehlercodes](ereignisse-und-fehlercodes.md).

## Einen Fehler finden

1. Letzte Zeile der Protokolldatei lesen. Steht dort `LAUF.FERTIG`, war der Lauf erfolgreich.
2. Steht dort `LAUF.FEHLGESCHLAGEN`: die Zeile davor mit Stufe `error` nennt den Fehlercode, was
   passiert ist, wo (`wo`), die wahrscheinliche Ursache (`ursache`) und was zu tun ist (`was_tun`).
3. Fehlt beides: der Lauf wurde von außen beendet (Platte weg, Rechner aus, Prozess beendet).
   Die letzte `STUFE.START`- oder `STUFE.LEBT`-Zeile zeigt, wo er war.
4. `STUFE.ABGEBROCHEN` nennt die Stufe, in der ein Fehler auftrat.
5. `PROGRAMM.ABSTURZ` ist ein Fehler im Programm selbst, mit Quelltextstelle.

Beispiel:

```text
[    0.087s] ERROR INDEX.CONTIG_FEHLT: Contig chrZZ steht nicht in der Referenz
             ursache=Name falsch geschrieben oder andere Referenz (z. B. 20 statt chr20)
             was_tun=Namen in der FASTA-Datei prüfen (Zeilen, die mit > beginnen)
             wo={"contigs_in_datei":1,"datei":"chr20.fa"}
[    0.089s] ERROR STUFE.ABGEBROCHEN: index.fasta_lesen
[    0.089s] ERROR LAUF.FEHLGESCHLAGEN: Lauf fehlgeschlagen
```

## Stufen der Ausführlichkeit

`error`, `warn`, `info` (Standard), `debug`, `trace`. Bei `genome index build` über
`--log-stufe`. Zeilen unterhalb der gewählten Stufe kosten keine Laufzeit.

## Schutz der Daten

Protokolle enthalten keine Sequenzen und keine Genotypen, nur Zahlen, Pfade und Namen.

## Geplant (Plan Abschnitt 15)

- Prüfsummen der Eingabedateien im Laufkopf, Spitzen-RAM und freier Plattenplatz (braucht `gc-sys`).
- Bilanz je Stufe (Lesungen hinein = heraus + begründet verworfen) ab `fq2bam`.
- `--trace-read <Name>`: der ganze Entscheidungsweg einer einzelnen Lesung.
- Protokollansicht und Diagnosepaket in der App.
- Für jeden Fehlercode ein Test, der ihn auslöst.
