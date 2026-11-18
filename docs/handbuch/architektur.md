# Architektur

Stand: gebaut (Grundgerüst), Erweiterung nach Plan Abschnitt 4 und 12.

## Crates

| Crate | Aufgabe | hängt ab von |
|---|---|---|
| `gc-io` | FASTQ, FASTA, BAM, BGZF lesen und schreiben | noodles, flate2 (zlib-rs) |
| `gc-index` | Referenz als 2-Bit-Folge, Suffix-Array, FM-Index | `gc-log`, flate2 (nur CRC-32) |
| `gc-log` | Ereignisstrom und Lauf-Protokoll | serde_json |
| `gc-schema` | Ergebnis-Hülle | serde, serde_json |
| `genome` | das Programm mit allen Unterbefehlen | alle oben, clap, anyhow |

Geplant: `gc-sys` (Systemnahes), `gc-ref`, `gc-stats`, `gc-bench`, `gc-report` und die Werkzeug-
Crates aus Plan Abschnitt 4.1.

## Schichten

Die Schichtung wird durch die Abhängigkeiten der Crates erzwungen: eine Crate kann nur nutzen,
was in ihrer `Cargo.toml` steht. Fachcode (`gc-index`, `gc-io`) kennt weder die Kommandozeile noch
die Oberfläche. Das Programm `genome` verbindet nur.

## Regeln für macOS, Linux und Windows

- Nur reines Rust oder plattformneutrale Crates. Systemnahes (Speicherabbild, Dateisperren,
  RAM-Größe, Schlaf verhindern) kommt in die Crate `gc-sys`, nie in den Fachcode.
- Pfade immer als `PathBuf`.
- Dateien werden mit Zeilenende `\n` geschrieben; beim Lesen wird `\r\n` vertragen (FASTQ, FASTA).
- Binärformate sind ausdrücklich little-endian ([`*.gcfm`](formate/gcfm.md)).
- Das Rechenergebnis ist auf allen Systemen gleich. Komprimierte Bytes einer BAM können je nach
  Prozessor abweichen, der Inhalt nicht.
- Ergebnisse hängen nicht von der Kernzahl ab. Beispiel: der blockweise Index-Bau liefert mit
  1 bis 4 Kernen dieselbe Datei ([Test](verfahren/blockweiser-bau.md)).

## Unsicherer Code

`unsafe` kommt genau an einer Stelle vor: beim blockweisen Index-Bau schreiben mehrere Kerne in
getrennte, vorab zugeteilte Plätze desselben Blocks (`gc-index/src/blockwise.rs`, `fill_block`).
Die Begründung steht dort am Code.

## Atomares Schreiben

Ergebnisdateien werden erst unter einem Nachbarnamen (`.teil`) geschrieben und dann umbenannt.
Ein Abbruch hinterlässt nie eine halbe Datei unter dem endgültigen Namen.
