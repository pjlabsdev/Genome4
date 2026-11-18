# Arbeitsregeln für Genome 4.0

- Alles in Rust. Keine fremden Programme zur Laufzeit (kein conda, Docker, Java, Python, Nextflow).
- Muss auf macOS, Linux und Windows laufen. CI auf allen drei Systemen muss grün sein.
  Systemnahes nur in einer eigenen Crate (`gc-sys`), nie im Fachcode. Pfade immer als `PathBuf`.
- Kommunikation mit pj auf Deutsch, Du-Form.
- Commits ohne `Co-Authored-By` oder andere KI-Urheberhinweise. Push nur auf pj's Wort.
- Absicht vor dem Bauen: größere Schritte erst nach pj's Freigabe.
- Jede Genauigkeitsaussage braucht einen Benchmark gegen eine Wahrheit (GIAB, HPRC, …)
  mit vorher eingefrorener Version. Kein „besser“ ohne Messung.
- Ausgaben als JSON mit Schema `genome.<werkzeug>.vN` über `gc-schema`.
- Texte in Oberfläche und Berichten: keine Gedankenstriche (— / –).
- Protokollierung über `gc-log` (Plan Abschnitt 15): jeder Lauf mit Laufkopf, Stufen, Entscheidungen
  und Abschluss; Fehler mit festem Code, Ort, Ursache und Handlungsempfehlung. Keine Sequenzen im Protokoll.
- Jedes Werkzeug liefert seinen eigenen Bericht; der Gesamtbericht setzt nur zusammen (Plan Abschnitt 14).
- Handbuch unter `docs/handbuch` gehört zum Code: jedes Werkzeug, Verfahren, Dateiformat und jeder
  Protokollcode hat dort seine Seite, mit Quellen und Messung. Änderung am Programm und Handbuchseite
  kommen im selben Commit; der Test `crates/genome/tests/handbuch.rs` prüft den Zusammenhang.
