# FASTQ-Prüfung und Herkunft der Lesungen

Stand: gebaut (`gc-io`, `fastq.rs` und `readname.rs`).

## Lesen
- gzip wird am Dateianfang (`1f 8b`) erkannt, mehrteilige gzip-Dateien werden ganz gelesen.
- Zeilenenden `\n` und `\r\n` gelten gleich.
- Die Sequenz wird in Großbuchstaben gewandelt; `.` gilt als N.

## Was geprüft wird

Jeder Datensatz, vollständig. Ein Fehler bricht mit Datei und Zeilennummer ab.

| Prüfung | Meldung |
|---|---|
| Kopfzeile beginnt mit `@` | Kopfzeile beginnt nicht mit '@' |
| Name nicht leer | leerer Lesungsname |
| Datensatz vollständig | Datei endet nach der Kopfzeile (abgeschnitten?) / Qualitätszeile fehlt (abgeschnitten?) |
| Sequenz nur A, C, G, T, N | unerlaubtes Zeichen in der Sequenz |
| Trennzeile `+` vorhanden | Trennzeile '+' fehlt |
| Qualität so lang wie die Sequenz | Qualität hat q Zeichen, Sequenz s |
| Qualitätszeichen im Bereich 33 bis 126 | ungültiger Qualitätswert |

Gepaart zusätzlich:

| Prüfung | Fehler |
|---|---|
| beide Dateien gleich viele Lesungen | die längere und die kürzere Datei werden genannt |
| Namen des Paars stimmen überein (ohne `/1`, `/2` und ohne Kommentar) | Nummer des Paars und beide Namen |

## Herkunft aus dem Namen

Gebraucht für Read Groups und später für optische Duplikate.

| Plattform | Schema | gelesen |
|---|---|---|
| Illumina (CASAVA 1.8+) | `INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y`, optional ein achtes Feld (UMI) | Flowcell, Lane, Kachel, X, Y |
| MGI / DNBSEQ | `<FLOWCELL>L<lane>C<spalte>R<reihe><nummer>`, z. B. `V300012345L1C001R0020000123` | Flowcell, Lane, Feld aus Spalte und Reihe; keine Koordinaten |

Unbekannte Schemata liefern keine Herkunft; das ist kein Fehler.

## Kennzahlen
Lesungen, Basen, kürzeste und längste Lesung, Anteil Basen mit Qualität ab 30 (Phred+33),
GC-Anteil, N-Anteil.

## Grenzen
- Qualitätskodierung Phred+64 (sehr alte Illumina-Daten) wird nicht erkannt.
- Verschachtelte FASTQ (beide Lesungen in einer Datei) werden noch nicht als Paar gelesen.
