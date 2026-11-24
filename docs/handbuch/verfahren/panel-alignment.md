# Panel-Alignment: Reads gegen viele ähnliche Sequenzen einer Genfamilie

Stand: gebaut und geprüft (`gc-panel`). Übernommen aus hla-wgs 1.2.0, wo dieser Aligner seit 1.1
minimap2 ersetzt. Genutzt von hla-wgs; später auch für weitere Gen-Werkzeuge (KIR, PGx-Gene).

## Aufgabe
Ein Gen-Werkzeug richtet Reads nicht gegen das ganze Genom aus, sondern gegen ein Panel: alle
Allele eines Gens, seine Paraloge und Pseudogene. Diese Sequenzen unterscheiden sich oft nur an
wenigen Stellen. Gefragt ist für jeden Read die beste Stelle in jeder Sequenzgruppe (Gen) und ein
genaues Alignment, damit später ein Allelpaar gewählt und neue Varianten erkannt werden können.
Derselbe Aligner richtet auch Allele gegeneinander aus, etwa wenn ein kurzes Allel die Enden eines
längeren erhält.

## Verfahren
1. **Anker:** exakte k-mere (21 für kurze Reads) aus dem Index aller Panel-Sequenzen. Mit
   Minimierern wird nur jedes w-te k-mer gesucht; ist die beste Kette dann kürzer als ein Anteil
   der Read-Länge, wird mit allen k-meren nachgesucht. k-mere, die in einer Sequenz zu oft
   vorkommen (Wiederholungen), geben dort keine Anker.
2. **Ketten:** Anker werden wie in minimap2 zu Ketten verbunden (Lücke, Diagonalabstand,
   begrenzter Rückblick). Je Sequenzgruppe werden höchstens `per_group` Ketten ausgerichtet, damit
   ein Gen mit vielen Allelen andere nicht verdrängt.
3. **Alignment:** zwischen den Ankern global, über sie hinaus als Verlängerung mit z-drop, mit
   zweistufigen affinen Lückenkosten und linksbündigen Lücken. Große Lücken werden an eindeutigen
   13-meren geteilt; gleiche Rechenaufgaben werden über ihren Inhalt wiederverwendet.
4. **Rettung des Partner-Reads:** ein gebändertes lokales Alignment in einem Fenster neben dem
   Partner, nur an Diagonalen mit genügend gemeinsamen 11-meren.

Die Punktwerte folgen den Voreinstellungen von minimap2 (`sr`, `asm5`, `asm10`, `hifi`), damit
Schwellen, die auf minimap2-Werten gebaut wurden, ihre Bedeutung behalten.

## Gleiches Ergebnis auf jedem Rechner
Die beiden innersten Schleifen (lokales Alignment im Band, Verlängerung) laufen über
Antidiagonalen, deren Zellen nicht voneinander abhängen; der Übersetzer macht daraus
Vektorbefehle. Auf x86_64 werden sie über `gc_sys::with_avx2` mit AVX2 übersetzt, wenn der
Prozessor es kann. Es läuft derselbe Code, das Ergebnis ist gleich; nur das Tempo ändert sich.
Die zeilenweise Fassung bleibt zum Vergleich (`HLA_EXT_AD=0`, `HLA_RESCUE_AD=0`).

## Prüfung
- Selbstprüfungen in `crates/gc-panel/src/props.rs`: jedes Alignment auf Zufallsfolgen gegen ein
  vollständiges Alignment über die ganze Matrix, jeder gemeldete Pfad neu bewertet, gleiche
  Ergebnisse von Antidiagonal- und Zeilenfassung, Zuordnung unabhängig von der Kernzahl.
- In hla-wgs gegen minimap2 auf 80 HPRC- und 59 HGSVC3-Genomen mit assemblierten Haplotypen
  (Wahrheit): zwei Felder 2485 / 2489 gegen 2486 / 2489, gesicherte vier Felder exakt 1779 / 1798
  gegen 1778 / 1801; HGSVC3 zwei Felder 1833 / 1833 bei beiden; Zeit je Probe 22,2 s gegen 31,5 s
  (4 Kerne). Tabelle und Einzelheiten: hla-wgs, `docs/aligner.md`.
- Der Umzug nach Genome 4 ändert nichts: hla-wgs mit `gc-panel` schreibt auf 8 HPRC-Proben
  byte-gleiche Ergebnisdateien wie mit dem eingebauten Aligner.

## Schalter
Aus hla-wgs übernommen, damit dessen Läufe vergleichbar bleiben: `HLA_RESCUE_V`, `HLA_RESCUE_BAND`,
`HLA_EXT_AD`, `HLA_RESCUE_AD`, `HLA_ALNSTAT`, `HLA_ALNBIG`, `HLA_PROPS`. Bedeutung: hla-wgs,
`docs/manual/switches.md`.

## Quellen
- Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18),
  3094–3100. https://doi.org/10.1093/bioinformatics/bty191
- Suzuki H., Kasahara M. (2018). Introducing difference recurrence relations for faster semi-global
  alignment of long sequences. BMC Bioinformatics, 19(Suppl 1), 45.
  https://doi.org/10.1186/s12859-018-2014-8
- hla-wgs: https://github.com/pjlabsdev/hla-wgs, DOI https://doi.org/10.5281/zenodo.24367469
