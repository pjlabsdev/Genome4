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
| `STUFE.LEBT` | info | Lebenszeichen einer laufenden Stufe | Fortschritt, `sekunden`, `speicher_spitze_bytes` |
| `STUFE.ENDE` | info | Stufe fertig | Zahlen der Stufe, `sekunden`, `speicher_spitze_bytes` |
| `STUFE.ABGEBROCHEN` | error | Stufe wurde ohne Abschluss verlassen | `sekunden` |
| `PROGRAMM.ABSTURZ` | error | Fehler im Programm selbst | `wo.quelltext`, `ursache`, `was_tun` |

## Index (`genome index build`)

| Kennung | Stufe | Bedeutung | Was tun |
|---|---|---|---|
| `INDEX.BAU.PARAMETER` | info | Variante, Verfahren, Größe der Referenz, Contigs, N-Strecken | |
| `INDEX.VERFAHREN` | info | Entscheidung zwischen direktem und blockweisem Bau, mit Grund | |
| `ALIGN.INSERTGROESSE` | info | Paare: je Block die geschätzte Insertgröße für die vier Richtungen (Zahl der Paare, Perzentile, Mittel, Streuung, Grenzen für "richtiges Paar") und welche Richtungen benutzt werden | |
| `ALIGN.FERTIG` | info | Ausrichten beendet: Reads, geschriebene Zeilen, Blöcke, Kerne, Laufzeit, Reads je Sekunde; Sekunden für Lesen (mit Trimmen), Ausrichten und Schreiben, die gleichzeitig laufen, und Wartezeit des Ausrichtens auf das Lesen; unter `schnellweg` die Zähler des [Schnellwegs](verfahren/schnellweg.md), unter `strobeweg` die des [Strobemer-Wegs](verfahren/strobeweg.md) | |
| `GPU.BEREIT` | info | `--grafikeinheit an`: Die Grafikeinheit rechnet die Partnersuche; Name und Schnittstelle unter `geraet`, siehe [Grafikeinheit](verfahren/grafikeinheit.md) | |
| `GPU.FEHLT` | warn | `--grafikeinheit an`, aber keine nutzbare Grafikeinheit; der Prozessor rechnet, die Ausgabe ist dieselbe. Grund unter `ursache` | Treiber prüfen (Metal auf macOS, Vulkan unter Linux, DX12 unter Windows) oder ohne die Option laufen lassen |
| `GPU.FEHLER` | warn | Die Grafikeinheit ist während des Laufs ausgefallen; der Prozessor rechnet diesen Block, die Ausgabe ist dieselbe | Bei Wiederholung ohne `--grafikeinheit an` laufen lassen und `ursache` melden |
| `FQ2BAM.SPEICHERPLAN` | info | fq2bam hat das Ausdünnen der Suffix-Array-Stichproben nach dem Arbeitsspeicher gewählt (`--sa-ausduennen 0`): Arbeitsspeicher, Größe von Index und Stichproben, Sortierpuffer, geschätzter Strobe-Index, gewählter Faktor und erwarteter Bedarf | |
| `SORT.TEIL` | info | Sortieren: der Speicher war voll, ein sortierter Teil liegt auf der Platte (Zahl der Teile, Bytes auf der Platte, Datensätze bisher) | |
| `SORT.FERTIG` | info | Sortieren beendet: Datensätze, Bytes, Teile, Laufzeit, Datensätze je Sekunde | |
| `MARKDUP.FERTIG` | info | Duplikate markiert: gelesen, bewertet, Duplikate mit und ohne Partner, Anteil, geschätzte Bibliotheksgröße, Spitze des Fensters, Laufzeit | |
| `BAMINDEX.FERTIG` | info | BAM-Index geschrieben: Datensätze, davon ohne Contig, Contigs mit Daten, Fächer, Größe des Index, Laufzeit | |
| `TRIM.ADAPTER` | info | `fq2bam --trimmen` und `trim`: je Seite (`read1`, `read2`) die erkannte Adapterfolge oder `keiner`, mit Grund, Name des bekannten Adapters, Weg der Erkennung und Zahl der angesehenen Reads ([Trimmen](verfahren/trimmen.md)) | |
| `TRIM.POLYG` | info | `fq2bam --trimmen` und `trim`: PolyG kürzen `an` oder `aus`, mit dem ersten Read-Namen, nach dem fastp entscheidet | |
| `TRIM.FERTIG` | info | `fq2bam --trimmen` und `trim`: Paare gesamt, behalten, verworfen nach Grund (Qualität, zu kurz, zu viele N, Adapter-Dimer), PolyG und Adapter gekürzt; behalten plus verworfen ergibt die Paare | |
| `BERICHT.FERTIG` | info | Überblick über eine BAM geschrieben: Pfade der HTML-Seite und der JSON-Datei, Schema, ob ein Protokoll von `fq2bam` eingeflossen ist | |
| `FQ2BAM.FERTIG` | info | `fq2bam` beendet: der ganze Bericht (Reads, Sortieren, Duplikate, Zählung nach Flags, Index, Zeiten je Schritt) | |
| `INDEX.STICHPROBE.RUNDE` | debug | eine Verfeinerungsrunde der Stichproben-Ränge | |
| `INDEX.PRUEFUNG.STAND` | debug | Fortschritt der Selbstprüfung | |
| `ALIGN.LESEGRUPPE_UNGUELTIG` | error | die Angabe bei `--lesegruppe` ist keine gültige `@RG`-Zeile (kein `@RG` am Anfang, echte Tabulatoren, kein `ID`, Feld ohne `XX:Wert`, Zeilenende) | Form `'@RG\tID:kennung\tSM:probe'`, Tabulatoren als `\t` geschrieben |
| `SORT.FEHLGESCHLAGEN` | error | Sortieren abgebrochen (Eingabe nicht lesbar oder abgeschnitten, Platte voll, Teil beschädigt); es gibt keine Ausgabedatei | Eingabe mit `samtools quickcheck` prüfen, freien Platz im Temp-Ordner und am Ziel prüfen, neu starten |
| `MARKDUP.FEHLGESCHLAGEN` | error | Markieren abgebrochen (Eingabe nicht nach Ort sortiert, Partner-Angaben `MC` oder `ms` fehlen, Eingabe abgeschnitten, Platte voll); es gibt keine Ausgabedatei | erst `genome sort`; Partner-Angaben mit `genome index sam --partner-angaben` oder `samtools fixmate -m` |
| `BAMINDEX.FEHLGESCHLAGEN` | error | BAM-Index nicht geschrieben (BAM nicht nach Ort sortiert, abgeschnitten, Position über 512 Mbp) | erst `genome sort`; die BAM mit `samtools quickcheck` prüfen |
| `BERICHT.FEHLGESCHLAGEN` | error (`bam-bericht`), warn (`fq2bam`) | Überblick nicht geschrieben (BAM nicht lesbar oder nicht nach Ort sortiert, Platte voll). In `fq2bam` bleibt die BAM gültig | BAM mit `samtools quickcheck` prüfen, Platz prüfen, mit `genome bam-bericht` nachholen |
| `FQ2BAM.FEHLGESCHLAGEN` | error | `fq2bam` abgebrochen (FASTQ fehlerhaft oder ungleich lang, Lesegruppe ungültig, Index nicht lesbar, Platte voll); es gibt keine Ausgabe-BAM | Ursache im Feld `ursache`; FASTQ mit `genome fastq-check` prüfen, Platz im Temp-Ordner und am Ziel prüfen |
| `INDEX.CONTIG_FEHLT` | error | ein mit `--contig` gewähltes Contig steht nicht in der FASTA | Namen in der FASTA prüfen (z. B. `20` statt `chr20`) |
| `INDEX.BAU_FEHLGESCHLAGEN` | error | Bau nicht möglich: leere Referenz, zu groß für das gewählte Verfahren, oder interner Fehler | steht im Feld `was_tun` |
| `INDEX.PRUEFUNG_FEHLGESCHLAGEN` | error | Stichproben- oder Vollprüfung hat eine Abweichung gefunden | Index nicht verwenden, Protokoll aufbewahren und melden |

Stufen des Index-Baus (Name in `msg` von `STUFE.*`): `index.fasta_lesen`, `index.text`,
`index.suffix_array`, `index.suffix_array.pruefung`, `index.bwt`, `index.selbstpruefung` (direkter
Bau); `index.stichprobe`, `index.faecher`, `index.bloecke`, `index.stichprobenpruefung`
(blockweiser Bau); `index.vollpruefung`, `index.schreiben`.

Stufen des Sortierens: `sort.lesen` (lesen, im Speicher sortieren, Teile schreiben),
`sort.schreiben` (zusammenführen und BAM schreiben). Stufe des Markierens: `markdup.durchlauf`. Stufe des BAM-Index: `bamindex.lesen`. Stufen von `fq2bam`: `fq2bam.trimmen_vorbereiten` (nur mit `--trimmen`), `fq2bam.index_laden`, `fq2bam.strobe_index` (nur mit `--strobeweg an`),
`fq2bam.ausrichten`, `fq2bam.zusammenfuehren`, danach `bamindex.lesen` (mit dem Überblick). Stufe von `bam-bericht`: `bericht.lesen`. Stufen von `trim`:
`fq2bam.trimmen_vorbereiten`, `trim.schreiben`.

## Nur in Tests

`TEST.WAHL`, `TEST.DEBUG`, `TEST.FEHLER` kommen nur in den Tests von `gc-log` vor.
