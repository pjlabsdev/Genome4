# Ereignisse und Fehlercodes

Stand: vollständig für den gebauten Teil. Ein Test prüft, dass jede Kennung aus dem Quelltext hier
steht.

Kennungen sind fest: sie ändern ihre Bedeutung nicht. Neue kommen hinzu.

## Lauf und Stufen (alle Werkzeuge)

| Kennung | Stufe | Bedeutung | Felder |
|---|---|---|---|
| `LAUF.START` | info | Lauf beginnt | `befehl`, `parameter`, `system`, `architektur`, `kerne`, `arbeitsordner` |
| `LAUF.EINGABE` | info | eine Eingabedatei | `rolle`, `pfad`, `bytes`, `geaendert` |
| `LAUF.FERTIG` | info | Lauf erfolgreich beendet | Zusammenfassung des Werkzeugs |
| `LAUF.FEHLGESCHLAGEN` | error | Lauf mit Fehler beendet | `fehler` |
| `STUFE.START` | info | Stufe beginnt (Name in `msg`) | je Stufe |
| `STUFE.LEBT` | info | Lebenszeichen einer laufenden Stufe | Fortschritt, `sekunden` |
| `STUFE.ENDE` | info | Stufe fertig | Zahlen der Stufe, `sekunden` |
| `STUFE.ABGEBROCHEN` | error | Stufe wurde ohne Abschluss verlassen | `sekunden` |
| `PROGRAMM.ABSTURZ` | error | Fehler im Programm selbst | `wo.quelltext`, `ursache`, `was_tun` |

## Index (`genome index build`)

| Kennung | Stufe | Bedeutung | Was tun |
|---|---|---|---|
| `INDEX.BAU.PARAMETER` | info | Variante, Verfahren, Größe der Referenz, Contigs, N-Strecken | |
| `INDEX.VERFAHREN` | info | Entscheidung zwischen direktem und blockweisem Bau, mit Grund | |
| `ALIGN.INSERTGROESSE` | info | Paare: je Block die geschätzte Insertgröße für die vier Richtungen (Zahl der Paare, Perzentile, Mittel, Streuung, Grenzen für "richtiges Paar") und welche Richtungen benutzt werden | |
| `ALIGN.FERTIG` | info | Ausrichten beendet: Reads, geschriebene Zeilen, Blöcke, Kerne, Laufzeit, Reads je Sekunde | |
| `INDEX.STICHPROBE.RUNDE` | debug | eine Verfeinerungsrunde der Stichproben-Ränge | |
| `INDEX.PRUEFUNG.STAND` | debug | Fortschritt der Selbstprüfung | |
| `ALIGN.LESEGRUPPE_UNGUELTIG` | error | die Angabe bei `--lesegruppe` ist keine gültige `@RG`-Zeile (kein `@RG` am Anfang, echte Tabulatoren, kein `ID`, Feld ohne `XX:Wert`, Zeilenende) | Form `'@RG\tID:kennung\tSM:probe'`, Tabulatoren als `\t` geschrieben |
| `INDEX.CONTIG_FEHLT` | error | ein mit `--contig` gewähltes Contig steht nicht in der FASTA | Namen in der FASTA prüfen (z. B. `20` statt `chr20`) |
| `INDEX.BAU_FEHLGESCHLAGEN` | error | Bau nicht möglich: leere Referenz, zu groß für das gewählte Verfahren, oder interner Fehler | steht im Feld `was_tun` |
| `INDEX.PRUEFUNG_FEHLGESCHLAGEN` | error | Stichproben- oder Vollprüfung hat eine Abweichung gefunden | Index nicht verwenden, Protokoll aufbewahren und melden |

Stufen des Index-Baus (Name in `msg` von `STUFE.*`): `index.fasta_lesen`, `index.text`,
`index.suffix_array`, `index.suffix_array.pruefung`, `index.bwt`, `index.selbstpruefung` (direkter
Bau); `index.stichprobe`, `index.faecher`, `index.bloecke`, `index.stichprobenpruefung`
(blockweiser Bau); `index.vollpruefung`, `index.schreiben`.

## Nur in Tests

`TEST.WAHL`, `TEST.DEBUG`, `TEST.FEHLER` kommen nur in den Tests von `gc-log` vor.
