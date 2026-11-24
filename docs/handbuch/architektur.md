# Architektur

Stand: gebaut (Grundgerüst), Erweiterung nach Plan Abschnitt 4 und 12.

## Crates

| Crate | Aufgabe | hängt ab von |
|---|---|---|
| `gc-io` | FASTQ, FASTA, BAM, BGZF lesen und schreiben | noodles, flate2 (zlib-rs) |
| `gc-index` | Referenz als 2-Bit-Folge, Suffix-Array, FM-Index | `gc-log`, flate2 (nur CRC-32) |
| `gc-bam` | BAM-Datensätze nach Ort sortieren, Duplikate markieren, Index `.bai` bauen, nach Flags zählen, Kopf anpassen | flate2 |
| `gc-align` | Seeds sammeln und orten, Ketten; später Alignment | `gc-index` |
| `gc-panel` | Panel-Aligner: Reads und Allele gegen viele ähnliche Sequenzen einer Genfamilie, aus hla-wgs ([Verfahren](verfahren/panel-alignment.md)) | `gc-sys`, ahash |
| `gc-log` | Ereignisstrom und Lauf-Protokoll | serde_json |
| `gc-schema` | Ergebnis-Hülle | serde, serde_json |
| `gc-sys` | Systemnahes je Prozessor und Betriebssystem: Vorab-Laden in den Cache (`prefetch`) für x86_64 und aarch64, Code mit AVX2 übersetzen, wenn der Prozessor es kann (`with_avx2`) | keine |
| `genome` | das Programm mit allen Unterbefehlen | alle oben, clap, anyhow |

Geplant: `gc-ref`, `gc-stats`, `gc-bench`, `gc-report` und die Werkzeug-
Crates aus Plan Abschnitt 4.1.

## Schichten

Die Schichtung wird durch die Abhängigkeiten der Crates erzwungen: eine Crate kann nur nutzen,
was in ihrer `Cargo.toml` steht. Fachcode (`gc-index`, `gc-io`) kennt weder die Kommandozeile noch
die Oberfläche. Das Programm `genome` verbindet nur.

## Regeln für macOS, Linux und Windows

- Nur reines Rust oder plattformneutrale Crates. Systemnahes (Speicherabbild, Dateisperren,
  RAM-Größe, Schlaf verhindern, Befehle für einen bestimmten Prozessor) kommt in die Crate
  `gc-sys`, nie in den Fachcode. Was dort steht, ändert kein Ergebnis, nur das Tempo oder den
  Umgang mit dem System; auf einem Prozessor ohne passenden Befehl geschieht nichts.
- Pfade immer als `PathBuf`.
- Dateien werden mit Zeilenende `\n` geschrieben; beim Lesen wird `\r\n` vertragen (FASTQ, FASTA).
- Binärformate sind ausdrücklich little-endian ([`*.gcfm`](formate/gcfm.md)).
- Das Rechenergebnis ist auf allen Systemen gleich. Komprimierte Bytes einer BAM können je nach
  Prozessor abweichen, der Inhalt nicht.
- Ergebnisse hängen nicht von der Kernzahl ab. Beispiel: der blockweise Index-Bau liefert mit
  1 bis 4 Kernen dieselbe Datei ([Test](verfahren/blockweiser-bau.md)).

## Unsicherer Code

Im Fachcode kommt `unsafe` genau an einer Stelle vor: beim blockweisen Index-Bau schreiben mehrere
Kerne in getrennte, vorab zugeteilte Plätze desselben Blocks (`gc-index/src/blockwise.rs`,
`fill_block`). Die Begründung steht dort am Code. Dazu kommen die systemnahen Befehle in `gc-sys`
(`prefetch`, `with_avx2`), jeweils mit Begründung am Code.

## Atomares Schreiben

Ergebnisdateien werden erst unter einem Nachbarnamen (`.teil`) geschrieben und dann umbenannt.
Ein Abbruch hinterlässt nie eine halbe Datei unter dem endgültigen Namen.
