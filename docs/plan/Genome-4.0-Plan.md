# Genome 4.0: eigene Rust-Werkzeuge für jeden Analyseschritt

Stand 18.11.2026, von pj gestartet (Ergänzungen: ganze App in Rust, neues UX; erster Schritt FASTQ → BAM, Abschnitt 0; läuft künftig auch auf Linux und Windows, Abschnitt 12.3). Noch nichts gebaut, nichts geändert.

## Kurzfassung

- Genome 3.0 ruft heute rund 30 fremde Werkzeuge über conda, Docker, Java und Nextflow auf. Genau dort entstehen die meisten Ausfälle (VEP-Segfaults, Manta nur über Docker, GATK-Versionsabfrage ohne Ende, CNVkit-SIGBUS von der SSD, 27,6 GB VEP-Cache).
- Genome 4.0 ersetzt sie durch **eine eigene Rust-Werkzeugfamilie `genome-core`**: signierte Programme, die direkt in der App liegen, ohne conda, ohne Docker, ohne Java, ohne Nextflow. Jedes Programm liefert JSON nach festem Schema, wie hla-wgs und brca12-wgs.
- Maßstab ist das HLA-Vorgehen: Wahrheitssätze aus Assemblies und Laborergebnissen, Version vor jeder Testrunde einfrieren, ungesehene Proben, Vergleich mit dem Original auf denselben Lesungen. Ein Werkzeug ersetzt das alte erst, wenn es dieses Tor besteht.
- **Ehrlich:** In acht Bereichen ist "besser als das Original" realistisch und messbar (Annotation, QC, PGx, KIR, Repeats, Lp(a), Haplogruppen, Spezialgene). Bei SV/CNV ist es mit Pangenom-Genotypisierung erreichbar, aber aufwendig. Beim Variant Calling ist GATK-Niveau gut erreichbar, DeepVariant-Niveau nur mit eigenem neuronalen Modell und viel Arbeit. Beim Aligner lohnt ein Neuschreiben nicht; wir binden minimap2/BWA-MEM2 als Bibliothek in ein eigenes Programm ein (wie bei HLA) und forschen erst später an einem Pangenom-Aligner.
- "Komplett unabhängig" heißt: keine fremden Programme mehr installieren. Referenzdaten (GRCh38, ClinVar, gnomAD, PharmVar, IMGT/HLA, PGS Catalog) bleiben fremde Daten. Die spiegeln wir mit Prüfsummen (wie genome-plus).
- **Alles in Rust, auch die App:** Genome 4.0 ist kein Swift-Programm mehr, sondern eine native Rust-App (GPU-gerenderte Oberfläche, Berichte mit Typst, Lokalisierung mit Fluent). Analysekern und Oberfläche teilen dieselben Crates. Abschnitt 12.
- **Neues UX:** eine zoombare Genom-Karte als Startfläche (vom Chromosom bis zur einzelnen Lesung), Befunde mit sichtbarer Beweiskette, Ergebnisse, die schon während der Analyse erscheinen, und ein Bericht, der in der App genau so aussieht wie im PDF. Abschnitt 13.
- **Erster Schritt (Festlegung pj): FASTQ → BAM komplett in Rust**, ein Programm `genome fq2bam` vom FASTQ bis zur indizierten, duplikatmarkierten BAM mit QC. Ausgearbeitet in Abschnitt 0.
- Danach: schnelle sichere Gewinne (Annotation, QC-Rest, Haplogruppen, HLA-Einbau), Spezialgene mit wissenschaftlichem Vorsprung (PGx, KIR, Repeats, Lp(a)), Variant Calling und SV/CNV, zuletzt Forschung (neuronaler Caller, Pangenom-Alignment).

---

## 0. Erster Schritt: FASTQ → BAM in Rust

Festlegung pj (18.11.): Genome 4.0 beginnt damit, wie aus einer FASTQ eine BAM wird. Ziel ist ein einziges Programm `genome fq2bam`, das ohne fremde Werkzeuge aus gepaarten FASTQ-Dateien eine sortierte, duplikatmarkierte, indizierte BAM (oder CRAM) mit vollständiger Qualitätsauswertung macht, gleich genau wie BWA-MEM2 + samtools und auf dem Mac schneller.

### 0.1 Heute und morgen

| Schritt | heute (3.x) | `genome fq2bam` |
|---|---|---|
| FASTQ prüfen, trimmen | fastp, FastQC (eigene Läufe, Zwischen-FASTQ auf Platte) | im Lesestrom, keine Zwischendatei |
| Alignment | BWA bzw. minimap2, SAM als Text | BWA-MEM-Verfahren in Rust, Datensätze direkt im Speicher |
| Sortieren | samtools sort (eigener Lauf) | im selben Prozess, Zwischenstände komprimiert auf interner SSD |
| Duplikate | samtools fixmate + markdup (zwei Läufe) | beim Zusammenführen der sortierten Teile |
| Index, Kennzahlen | samtools index, mosdepth, flagstat (drei Läufe) | im selben Schreibdurchgang |
| Ergebnis | BAM + verstreute Logdateien | BAM/CRAM + Index + `fq2bam.json` (Schema `genome.fq2bam.v1`) mit Kennzahlen, Versionen, Prüfsummen |

Heute liest und schreibt die Kette die Daten fünf- bis sechsmal komplett. `fq2bam` liest das FASTQ einmal und schreibt die BAM einmal.

### 0.2 Die Stufen im Programm

Alle Stufen laufen gleichzeitig als Datenstrom über alle Kerne:

1. **Einlesen:** FASTQ und FASTQ.gz (gepaart, verschachtelt, mehrere Lanes), auch unalignierte oder bereits alignierte BAM/CRAM als Eingang (zum Neu-Alignieren, z. B. pj's PJ.bam). Read Groups aus dem Kopf der Lesungen, für Illumina- und MGI/DNBSEQ-Schreibweise.
2. **Prüfen, bevor gerechnet wird:** Paare passen zusammen, Qualitätskodierung, Längen, abgeschnittene Dateien. Fehler brechen sofort mit klarer Meldung ab, nicht nach Stunden.
3. **Trimmen:** Adapter über die Überlappung der Paare (fastp-Prinzip) plus bekannte Illumina- und MGI-Adapter; Poly-G-Enden (zweifarbige Geräte). Kein hartes Qualitätstrimmen als Standard, weil das Alignment ohnehin weich clippt; die Wirkung beider Varianten wird gemessen und im Bericht begründet.
4. **Referenz und Index:** GRCh38 Analyse-Set mit Decoys (hs38d1, wie pj's Daten), mit oder ohne ALT. Suffix-Array (SA-IS) und FM-Index in Rust gebaut, als fertiges Datenpaket aus genome-reflib geladen und per Speicherabbild genutzt. Zwei Varianten: **schnell** (BWA-MEM2-artige Belegungstabellen, mehr RAM) und **kompakt** (unter 6 GB, für Macs mit 16 GB); Genome wählt nach Speicher.
5. **Alignment (BWA-MEM-Verfahren):** SMEM-Seeds mit Reseeding, Chaining und Chain-Filter, gebänderte Smith-Waterman-Erweiterung mit SIMD (NEON auf ARM, AVX2/AVX-512 auf x86, zur Laufzeit gewählt) über viele Lesungen gleichzeitig, Insert-Größen je Block geschätzt, Mate-Rescue, Paarbewertung, MAPQ, supplementäre Alignments, alternative Treffer (XA), ALT-bewusste Behandlung. **Ergebnis unabhängig von der Kernzahl** (feste Blockgröße, wie BWA `-K`), damit Läufe reproduzierbar sind.
6. **Sortieren:** Datensätze gehen ohne SAM-Text in Positionsfächer im Speicher; ist das RAM-Budget voll, werden sortierte, komprimierte Teile auf die interne SSD geschrieben und am Ende zusammengeführt.
7. **Duplikate:** beim Zusammenführen markiert (Position, Orientierung, Mate, Bibliothek), optische Duplikate über Kachel-Koordinaten aus dem Lesungsnamen (Illumina und MGI).
8. **Schreiben:** BAM mit paralleler BGZF-Kompression (zlib-rs, reines Rust) oder CRAM 3.1; BAI/CSI/CRAI im selben Durchgang; MD/NM-Tags; `@PG` mit exakter Version und allen Parametern; Prüfsumme.
9. **Kennzahlen unterwegs:** Anteil aligniert und richtig gepaart, Duplikatrate, Insert-Größen, MAPQ-Verteilung, Abdeckung je Chromosom und Region (mosdepth-gleich), GC-Verzerrung, Geschlecht aus X/Y-Abdeckung, Fehlerrate je Zyklus.

Fortschritt geht als Ereignisse an die App (Stufe, Lesungen pro Sekunde, Restzeit), Abbruch ist jederzeit sauber, ein Neustart setzt an fertigen Teilen wieder an.

### 0.3 Validierung

Die Frage ist nicht "gleich wie BWA-MEM2", sondern "landet jede Lesung dort, wo sie herkommt, und werden die Varianten danach mindestens so gut gefunden".

1. **Simulierte Lesungen mit bekannter Herkunft:** aus der diploiden HG002-Q100-Assembly und HPRC-Assemblies, mit Fehlerprofilen, die aus echten Illumina- und MGI-Läufen gelernt sind. Jede Lesung kennt ihre wahre Position, also wird die Fehlplatzierungsrate exakt messbar, getrennt nach GIAB-Schichten (Segmentduplikationen, Homopolymere, MHC, KIR, CYP2D6, LPA). Zusätzlich die **MAPQ-Kalibrierung**: stimmt "MAPQ 30" mit einer tatsächlichen Fehlerrate von 1:1000 überein. Dieselben Lesungen gehen durch BWA-MEM2 und minimap2.
2. **Echte Daten gegen BWA-MEM2 + samtools:** GIAB HG001 bis HG007 (Illumina 30x), HG002 auf MGI/DNBSEQ, 1000G-30x-Proben, pj's eigene Lesungen. Gemessen: Anteil gleicher Primärpositionen, jede Abweichungsgruppe nach Region erklärt, Duplikat-Übereinstimmung, flagstat-Zahlen.
3. **Wirkung auf die Varianten (der eigentliche Richter):** auf unserer BAM und auf der BWA-MEM2-BAM laufen dieselben eingefrorenen Caller (DeepVariant und GATK HC als neutrale Richter). hap.py gegen GIAB v4.2.1, HG002 Q100 und CMRG, alle Schichten.
4. **Reproduzierbarkeit:** byteidentische BAM bei 4, 8 und 16 Kernen und bei Wiederholung.
5. **Robustheit:** abgeschnittene Dateien, vertauschte Paare, gemischte Leselängen (2×100 bis 2×150), 10x bis 60x, BAM als Eingang.
6. **Testdesign wie bei HLA:** Version vor der Testrunde einfrieren; Entwicklung auf chr20 und HG001, Test auf den übrigen Proben.

### 0.4 Messziele

| Kennzahl | Ziel |
|---|---|
| Fehlplatzierung (simuliert, MAPQ ≥ 20) | höchstens so hoch wie BWA-MEM2, in jeder Schicht |
| MAPQ-Kalibrierung | genauer als BWA-MEM2 (gemessener Fehler je MAPQ-Stufe näher am Sollwert) |
| Gleiche Primärposition wie BWA-MEM2 (echt) | über 99,5 %, alle Abweichungen erklärt |
| SNV-/Indel-F1 danach (DeepVariant und HC als Richter) | nicht schlechter als mit BWA-MEM2-BAM (95-%-Intervall), in schweren Regionen besser angestrebt |
| Duplikate | Übereinstimmung mit samtools markdup über 99,9 % |
| Wanduhr FASTQ bis indizierte BAM, 30x, pj's Mac | höchstens die Zeit, die heute BWA-MEM2 allein braucht (Sortieren, Duplikate, Index, QC kosten nichts extra); Streckziel 1,5-mal schneller als die heutige Gesamtkette |
| Spitzen-RAM | unter 16 GB mit kompaktem Index |
| Zwischenspeicher auf Platte | höchstens 1,2-mal die fertige BAM |

Die Basiswerte der heutigen Kette messe ich zuerst auf deinem Mac; erst dann stehen hinter den Zielen echte Zahlen.

### 0.5 Meilensteine

| M | Inhalt | Prüfung (geht weiter nur, wenn bestanden) | Dauer |
|---|---|---|---|
| M1 | Workspace, gc-io: FASTQ/BAM lesen und schreiben, BGZF parallel; Basismessung der heutigen Kette | Rundreise identisch zu samtools | 1 Woche |
| M2 | SA-IS, FM-Index, beide Varianten | exakte Trefferzahlen wie BWA bei 1 Mio. Abfragen | 1 bis 2 Wochen |
| M3 | Seeds, Chaining, SIMD-Erweiterung, Einzel-Lesungen auf chr20 | über 99 % gleiche Positionen wie BWA-MEM2, Tempo auf Kurs (sonst Halt und Neubewertung) | 2 bis 3 Wochen |
| M4 | Paarung, MAPQ, Supplementär, ALT, ganzes Genom | simulierte Fehlplatzierung höchstens BWA-MEM2 | 2 Wochen |
| M5 | Sortieren, Duplikate, Index, QC verschmolzen; Trimmen; Einbau in Genome 3.x im Schattenbetrieb neben dem heutigen BWA-Weg | Duplikate und Kennzahlen gleich; App-Lauf auf pj's Probe | 2 Wochen |
| M6 | Eingefrorene Version, volle Testrunde (simuliert, GIAB, MGI, Varianten-Richter), Forschungsbericht, Repo mit DOI | alle Messziele aus 0.4 | 1 bis 2 Wochen |

Zusammen etwa 9 bis 12 Wochen. Größtes Risiko: die vielen kleinen Heuristiken von BWA-MEM (MAPQ, Rescue, Clipping-Strafen). Ziel ist darum nicht bitgleich mit BWA-MEM2, sondern gemessen gleich gut gegen die Wahrheit.

### 0.6 In der App

- In Genome 3.x kommt `fq2bam` als Alternative neben den heutigen Weg, zuerst im Schattenbetrieb (beide laufen, Abweichungen sichtbar). Nach M6 wird es Standard, BWA, samtools, fastp und FastQC fallen für diesen Schritt weg.
- In der Rust-App ist es der erste Knoten im Fluss-Graphen (Abschnitt 13.4): Lesungen pro Sekunde, Restzeit, und Kennzahlen wie Duplikatrate und Geschlecht erscheinen schon während des Laufs.

### 0.7 Veröffentlichung

Eigenes Repo (Name offen, z. B. `fq2bam-rs`), MIT/Apache-2.0, Zenodo-DOI, Forschungsbericht auf pjlabs.dev mit dem Simulations- und GIAB-Vergleich. Das BWA-MEM-Verfahren wird in README und Bericht klar zitiert (Li 2013; Vasimuddin et al. 2019), weil es eine Neuumsetzung ist, keine neue Methode.

---

## 1. Ausgangslage: Werkzeuge in Genome 3.0

Quelle: README 3.0, gvali-Prüfung vom 16.11., brca12-Benchmark, HLA-Arbeit. Die genaue Versionsliste aus dem Werkzeugkatalog der App prüfe ich in Welle 0 auf dem Mac nach (heute war der Ordner aus diesem Thread nicht lesbar).

| Bereich | Werkzeug heute | Laufzeitbasis | bekannte Probleme |
|---|---|---|---|
| FASTQ-QC, Trimmen | fastp, FastQC | conda, C++/Java | defekte fastp-Version als "Abgebrochen" (behoben) |
| Alignment | BWA, minimap2 | conda | lange Laufzeit, Wiederholungslogik |
| BAM/CRAM, Indizes | samtools, htslib, bgzip, tabix | conda | samtools nur in SSD-conda |
| Variant Calling | bcftools (App), GATK HaplotypeCaller (brca12) | conda, Java | GATK-JVM hängt bei Versionsabfrage |
| SV | Delly 1.3.3, Manta 1.6 | nativ / Docker | Manta nur über Docker |
| CNV | CNVkit (GATK gCNV im Gespräch) | conda, Python | SIGBUS beim Laden von der SSD |
| Annotation | Ensembl VEP 116, dbSNP, ClinVar, gnomAD | conda, Perl | Segfault (Exit 139), Cache 27,6 GB |
| Zusatzbefunde | ACMG-SF-Liste (eigen, Swift) | App | |
| HLA | HLA*LA, T1K | conda | HLA*LA bis 14 min, Phantom-DRB4; **eigener Ersatz fertig (pjlabshla v16)** |
| KIR | T1K (KIR-Modus), KIR-HLA-Liganden (eigen) | conda | KIR-Typisierung aus WGS generell schwach |
| Pharmakogenetik | Aldy + CPIC/FDA-Regeln (eigen) | conda, Python/ILP | stille Phasen, Lizenz von Aldy prüfen |
| Lp(a) KIV-2 | KILDA | Adapter | |
| Repeat-Expansionen | ExpansionHunter | conda | X/Y-Behandlung (behoben) |
| Telomere | TelomereHunter2 | conda, Python | |
| Haplogruppen | Haplogrep (mt), Y-Ableitung | conda/Java | |
| Mikroben | Kraken2, Bracken (EBV) | conda | große Datenbank |
| Probenqualität | mosdepth, VerifyBamID2, Somalier | conda | |
| Polygene Scores | eigen (Swift, PGS Catalog) | App | 128-MiB-Grenze bei großen Scores |
| Liftover | Ketten-Liftover | App/Tool | |
| SNP-Werkzeuge | SNP-Suche, Microarray-Export (eigen), MANE-Proteine über VEP | App + VEP | |
| BRCA1/2 | brca12-wgs (Nextflow, BWA, GATK, Manta, RTG) | Nextflow, Java, Docker | schwerste Abhängigkeitskette |
| Benchmarks | RTG vcfeval, hap.py | Java/Python | nur zur Validierung |

---

## 2. Was "besser" heißt: das Ablösetor

Ein eigenes Werkzeug ersetzt ein fremdes erst, wenn alle Punkte belegt sind:

1. **Genauigkeit:** auf allen Wahrheitssätzen des Bereichs mindestens gleich gut wie das Original (Nicht-Unterlegenheit mit 95-%-Konfidenzintervall), auf mindestens einer Hauptkennzahl messbar besser. Gemessen auf denselben Lesungen und derselben Hardware.
2. **Ungesehene Proben:** die Zahlen stammen aus Runden mit vorher eingefrorener Version. Entwicklungsproben und Testproben sind getrennt und stehen in Dateien im Repo.
3. **Tempo und Speicher:** schneller als das Original auf pj's Mac und höchstens 16 GB RAM (damit Genome auf jedem Apple-Silicon-Mac läuft).
4. **Robustheit:** 30x, 20x, 15x, 10x (Downsampling); Illumina und MGI/DNBSEQ (pj's eigene Daten sind MGI); BAM und CRAM; GRCh38 mit und ohne ALT; männlich und weiblich.
5. **Ehrlichkeit im Ergebnis:** jede Aussage trägt eine Sicherheitsangabe und die Alternativen, wie beim HLA-Typisierer. Was das Werkzeug nicht kann, steht in der Ausgabe.
6. **Reproduzierbarkeit:** byteidentische Ausgabe bei gleichem Input, Version und Datenstand; Prüfsummen im Manifest.
7. **Schattenbetrieb bestanden:** in der App liefen altes und neues Werkzeug parallel auf echten Proben, die Abweichungen sind erklärt.

Kennzahlen je Bereich: Recall, Precision, F1 (Varianten, nach GIAB-Schichten), Genotyp-Konkordanz, Kopienzahl-Treffer, Allel-Treffer je Feld (HLA/KIR/PGx), Mendel-Fehlerrate in Trios, Laufzeit (Wanduhr und CPU), Spitzen-RAM.

---

## 3. Ehrliche Einordnung je Bereich

| Bereich | besser als Original realistisch? | warum |
|---|---|---|
| Annotation (VEP) | **ja, deutlich** | Regelwerk mit klarer Wahrheit (HGVS, Konsequenz); Gewinn bei Tempo, Größe, Stabilität |
| QC (fastp, FastQC, mosdepth, VerifyBamID2, Somalier) | **ja** | einfache Algorithmen; Gewinn durch einen einzigen Lesedurchgang statt fünf |
| HLA | **ja, bewiesen** | 2-Feld 2937/2941 vs. T1K 2921/2941, DRB3/4/5-Kopien fast fehlerfrei |
| KIR | **ja, großer Abstand möglich** | heutige WGS-Lösungen schwach, Wahrheit aus HPRC-Assemblies verfügbar |
| Pharmakogenetik (CYP2D6 u. a.) | **ja** | Hybrid-/Kopienzahl-Allele sind die Schwäche; Wahrheit GeT-RM + Assemblies |
| Lp(a) KIV-2, SMN1/2 u. a. Spezialgene | **ja** | k-mer-/Paralog-Modelle, Wahrheit aus Assemblies |
| Repeat-Expansionen | **ja** | Wahrheit aus Assemblies und Coriell-Proben mit bekannten Expansionen |
| Haplogruppen mt/Y, Telomere, EBV, Liftover | **ja, aber kleiner Gewinn** | Original ist schon gut; Gewinn ist Integration und Tempo |
| Polygene Scores | **ja** | Engine ist eigen; Gewinn durch Ahnenkorrektur und vollständige Abdeckung aus WGS |
| BAM-Verarbeitung (Sortieren, Duplikate, Index) | **ja beim Tempo**, gleich bei Genauigkeit | Verschmelzen der Schritte spart Plattenzugriffe |
| SV | **teilweise** | Entdeckung neuer SV aus Kurzlesungen bleibt schwer; Genotypisierung bekannter SV aus dem Pangenom schlägt Delly/Manta klar |
| CNV | **teilweise** | Einzelproben-CNV ist begrenzt; gute GC/Mappability-Korrektur und Referenzpanel bringen mehr als CNVkit im WGS-Modus |
| Small-Variant-Calling | **GATK-Niveau ja, DeepVariant-Niveau schwer** | DeepVariant ist der Stand der Technik; erreichbar nur mit eigenem Modell, Training auf GIAB/Q100, eigener Inferenz auf Apple-Chip |
| Aligner | **nein, Neuschreiben lohnt nicht** | BWA-MEM2/minimap2 sind ausgereift; eigener Aligner wäre Jahre Arbeit für kaum messbaren Gewinn. Sinnvoll: als Bibliothek einbinden, später pangenom-bewusst |

---

## 4. Architektur: `genome-core`

### 4.1 Aufbau

Ein Rust-Workspace mit gemeinsamen Bausteinen und einzelnen Programmen:

| Crate | Aufgabe |
|---|---|
| `gc-io` | BAM/CRAM/FASTQ/VCF/BCF/bgzip/tabix/CSI lesen und schreiben, auf `noodles` (reines Rust, MIT); Regionen-Abfrage, paralleles Dekodieren |
| `gc-ref` | Referenzgenom, Ketten, MANE/Ensembl-Transkripte, Datenbank-Spiegel mit Prüfsummen |
| `gc-stats` | gemeinsame Statistik: EM, HMM, Wilson-Intervalle, Kalibrierung |
| `gc-bench` | Benchmark-Rahmen (siehe Abschnitt 5) |
| `gc-schema` | JSON-Schemata `genome.<werkzeug>.vN` und Manifest |
| Programme | `gqc`, `galign`, `gbam`, `gcall`, `gsv`, `gcnv`, `gannot`, `ghla` (= pjlabshla), `gkir`, `gpgx`, `gloci` (Lp(a), SMN, …), `gstr`, `ghap`, `gtel`, `gmicro`, `gpgs`, `glift` |

Grundsätze:
- **Ein Programm je Aufgabe, gemeinsame Bausteine.** Jedes Programm ist allein lauffähig (für die Veröffentlichung) und nutzt dieselben Crates (für die App).
- **Aufruf aus der App:** Die Rust-App linkt die Crates direkt (gleiche Typen, kein JSON-Umweg für Zustände), startet schwere Läufe aber als **Arbeitsprozess desselben Programms** (`genome worker <werkzeug>`). Begründung: ein Absturz oder Speicherüberlauf reißt die Oberfläche nicht mit, Abbruch ist sauber, und dieselben Werkzeuge bleiben als eigenständige Programme veröffentlichbar. Fortschritt kommt als typisierte Ereignisse über eine Pipe (behebt die "stillen Phasen", die Hermes bei Aldy fand).
- **macOS, Linux und Windows gleichwertig** (Abschnitt 12.3): SIMD für ARM (NEON) und x86 (AVX2/AVX-512) mit Auswahl zur Laufzeit, alle Kerne, GPU über wgpu (Metal, Vulkan, DirectX 12) für neuronale Teile. Kein Code darf nur auf einem System laufen.
- **Ausgabe:** JSON nach Schema plus Standardformate (VCF, BAM/CRAM, TSV), damit Ergebnisse auch außerhalb von Genome nutzbar sind.
- **Determinismus:** feste Zufallssaaten, stabile Sortierung (beim HLA-Port nötig gewesen, siehe typer_v16b).

### 4.2 Verteilung

- Je System ein signiertes Paket (macOS: notarisiertes `Genome.app`/DMG; Windows: signiertes MSI; Linux: AppImage und .deb/.rpm) mit einem Rust-Programm, das Oberfläche und alle Werkzeuge enthält; die veröffentlichten Einzelprogramme sind dieselben Crates mit eigener Kommandozeile. Keine Installation, kein conda, kein Docker. Das löst das Kernproblem aus dem GPT-6-Auftrag (zuverlässige Werkzeugverwaltung) an der Wurzel.
- Große Daten (Referenz, Annotation, Pangenom-Panels, Modelle) lädt Genome wie heute, aber als **versionierte, geprüfte Pakete** aus genome-plus/genome-reflib (Prüfsumme, atomarer Austausch, Wiederaufnahme).
- Zielgröße Annotationsdaten: unter 6 GB statt 27,6 GB VEP-Cache.

### 4.3 Arbeitsdaten

Zwischenergebnisse auf die interne SSD, nicht auf die ExtremeSSD, solange die sich unter Last trennt. Endergebnisse mit Prüfsumme auf den gewählten Speicherort.

---

## 5. Validierungsplattform `genome-bench`

Ein eigenes Repo mit Skripten, die jede Zahl im Bericht neu erzeugen. Aufgebaut auf dem, was beim HLA-Typisierer funktioniert hat.

### 5.1 Wahrheitssätze

| Wahrheit | wofür | Zugang |
|---|---|---|
| GIAB v4.2.1 HG001 bis HG007, GRCh38, mit GIAB-Schichten (schwierige Regionen, Homopolymere, Segmentduplikationen) | SNV/Indel | öffentlich (NIST FTP/S3) |
| GIAB HG002 T2T-Q100-Benchmark (kleine Varianten und SV) | SNV/Indel/SV in schweren Regionen | öffentlich |
| GIAB CMRG (273 medizinisch relevante, schwierige Gene) | klinisch wichtige Gene, PGx-Gene | öffentlich |
| Platinum Pedigree CEPH-1463 (vererbungsgeprüfte Wahrheit) | Varianten, SV, Repeats, Mendel | öffentlich |
| HPRC Release 2 (rund 230 Assemblies) + 1000G-30x-CRAMs | HLA, KIR, PGx, Lp(a), SMN, Repeats, SV: Wahrheit aus den Haplotypen ablesen | öffentlich, im Container erreichbar (bewährt) |
| HGSVC3 (65 Proben, Langlesungen) | SV | öffentlich |
| 1000G NYGC 30x, 3.202 Proben, 602 Trios | Populationsbreite, Mendel-Prüfung, Haplogruppen | öffentlich |
| 1000G-Labor-HLA (Sanger SBT) | HLA | bereits genutzt |
| GeT-RM-Konsens (CDC, Coriell) | PGx-Sternallele, viele Proben in 1000G | öffentlich |
| ClinVar mit HGVS, MANE Select | Annotation | öffentlich |
| pj's eigene Proben (PJ.bam, MGI) und HG1 | Praxistest, Schattenbetrieb | lokal auf dem Mac |

### 5.2 Regeln

- Entwicklungssatz, Testsatz und Reserve als Listen im Repo, vor dem ersten Lauf festgelegt.
- Jede Runde: Version einfrieren (Commit-Hash), dann erst auf ungesehene Proben. Fehler fließen nur in die nächste Version ein.
- Gleiche Lesungen, gleiche Hardware, gleiche Kernzahl für Original und eigenes Werkzeug. Originalversionen festgehalten.
- Vergleichswerkzeuge: hap.py/vcfeval für SNV/Indel, Truvari für SV, eigene Abgleiche für Allele. Diese Benchmarks bleiben fremd, weil ein Prüfer unabhängig vom Geprüften sein soll.
- 95-%-Wilson-Intervalle; McNemar-Test für gepaarte Vergleiche.
- Jede Runde endet mit einem Bericht wie `hla-typer/Ergebnis-Testsatz.md`.

### 5.3 Rechenorte

Container (4 Kerne, 15 GB) für Breitentests mit öffentlichen Daten; pj's Mac für Tempo-Messungen auf Zielhardware und für lokale Proben. EBI-FTP und NCBI sind aus dem Container gesperrt; diese Daten holt der Mac.

---

## 6. Die Werkzeuge im Einzelnen

Für jedes: was ersetzt wird, Ansatz, Validierung, Messziele, Risiko, Veröffentlichung. Basiswerte (Laufzeit, Speicher, Genauigkeit der Originale) messe ich in Welle 0 auf pj's Mac und HG002 30x; die Ziele sind deshalb relativ zum Original angegeben.

### 6.1 `gc-io` und `gbam`: BAM/CRAM-Verarbeitung
- **Ersetzt:** samtools (view, sort, index, fixmate, markdup, depth, faidx), htslib, bgzip, tabix, Teile von bcftools (norm, view, index), Picard-Funktionen.
- **Ansatz:** noodles als Grundlage. Kern ist ein **verschmolzener Durchgang**: Aligner-Ausgabe direkt sortieren, Duplikate markieren (samblaster-Logik, mate-aware), indizieren, Abdeckung und QC-Zahlen mitschreiben. Heute sind das vier bis fünf getrennte Läufe mit je einem kompletten Schreib-/Lesezyklus.
- **Validierung:** byteweiser Vergleich der Datensätze mit samtools-Ausgabe (gleiche Flags, gleiche Duplikate) auf 20 Proben; CRAM-Rundreise gegen samtools; Fuzzing mit kaputten Dateien.
- **Ziele:** Datensätze identisch zu samtools markdup; Wanduhr der Kette Alignment-Ende bis indizierte BAM halbiert; RAM unter 8 GB.
- **Risiko:** gering. CRAM 3.1-Kodierer sind Fleißarbeit.
- **Veröffentlichung:** als Crate in genome-core, kein eigenes Paper.

### 6.2 `gqc`: Qualität in einem Durchgang
- **Ersetzt:** fastp, FastQC, mosdepth, VerifyBamID2, Somalier, Sequenzierqualität.
- **Ansatz:** FASTQ-Prüfung und Adapter-Trimmen im Lesestrom vor dem Aligner (kein Zwischen-FASTQ); BAM-Kennzahlen (Abdeckung je Region, Insert-Größe, Duplikate, GC-Verlauf, Basenqualität, MAPQ) im Durchgang von `gbam`. Kontamination nach dem VerifyBamID2-Modell (Allelanteile an 1000G-SNPs mit Ahnen-PCs); Probenidentität und Geschlecht über einen festen SNP-Satz (Somalier-Prinzip); Verwandtschaft in Trios.
- **Validierung:** Abdeckung identisch zu mosdepth; Kontamination mit künstlich gemischten 1000G-Proben (0,5 bis 10 %) gegen VerifyBamID2; Identität und Verwandtschaft an den 602 1000G-Trios.
- **Ziele:** gleiche Zahlen wie die Originale (Abweichung unter 0,1 %); Kontaminationsschätzung mit kleinerem Fehler als VerifyBamID2 bei unter 2 %; QC kostet keine eigene Laufzeit mehr, weil sie im Durchgang mitläuft.
- **Risiko:** gering.
- **Veröffentlichung:** `gqc` als eigenes Repo mit DOI (kleiner Methodenbericht).

### 6.3 `galign`: Alignment
- **Ersetzt:** den Aufruf von BWA/minimap2 als fremdes Programm.
- **Ehrliche Einschätzung:** Ein neuer Alignment-Algorithmus, der BWA-MEM2 an Genauigkeit deutlich schlägt, ist nicht realistisch; das Feld hat 15 Jahre Optimierung hinter sich. Realistisch ist eine **reine Rust-Umsetzung des BWA-MEM-Verfahrens**, die gleich genau ist, auf Apple Silicon schneller läuft (BWA-MEM2 ist auf x86/AVX-512 optimiert, nicht auf NEON) und mit Trimmen, Sortieren, Duplikaten und QC in einem Prozess verschmilzt.
- **Ansatz in zwei Stufen:**
  1. **Stufe 1 (erster Schritt, Abschnitt 0):** BWA-MEM-Verfahren in Rust (FM-Index, SMEM-Seeds, Chaining, SIMD-Smith-Waterman, Paarung, MAPQ). minimap2 und BWA-MEM2 dienen nur noch als Vergleich im Benchmark.
  2. **Stufe 2 (Welle 5, Forschung):** pangenom-bewusstes Alignment gegen GRCh38 plus HPRC-Haplotypen in schwierigen Regionen (MHC, KIR, CYP2D6, LPA, SMN, Segmentduplikationen), mit Rückprojektion auf GRCh38. Vorbild: vg Giraffe + DeepVariant zeigte weniger Fehler in schweren Regionen. Hier liegt der echte wissenschaftliche Gewinn, nicht in der Grundgeschwindigkeit.
- **Validierung:** Konkordanz der Primärplatzierung mit BWA-MEM2 (Ziel über 99,9 %), vor allem aber Wirkung auf das Variant Calling (F1 nach GIAB-Schichten). Stufe 2: F1 in CMRG-Genen und Q100-schweren Regionen.
- **Ziele:** Stufe 1 gleiche Genauigkeit, Kette FASTQ bis fertige BAM mindestens 30 % schneller als heute. Stufe 2 messbar weniger Fehler in schweren Regionen bei höchstens 1,5-facher Laufzeit.
- **Risiko:** Stufe 1 gering, Stufe 2 hoch.
- **Veröffentlichung:** Stufe 2 wäre ein eigenes Paper.

### 6.4 `gcall`: kleine Varianten (SNV/Indel)
- **Ersetzt:** bcftools mpileup/call (App), GATK HaplotypeCaller (brca12-wgs).
- **Stand der Technik:** DeepVariant (CNN auf Pileup-Bildern) ist bei Kurzlesungen führend, GATK HC liegt deutlich dahinter bei Indels, bcftools noch weiter. DeepVariant läuft auf dem Mac nur auf der CPU und ist dort langsam.
- **Ansatz in drei Stufen:**
  1. **Haplotyp-Caller:** aktive Regionen, lokale De-Bruijn-Assemblierung, Pair-HMM (SIMD), Genotyp-Likelihoods, gVCF mit Referenz-Konfidenz. Ziel HC-Niveau. Löst sofort Java aus brca12-wgs.
  2. **Lernender Filter:** Gradient-Boosting auf Merkmalen je Kandidat (Allelbalance, Strang, Basenqualität, Kontext, Mappability), trainiert auf GIAB HG001/HG005-7, getestet auf HG002-4 (nie auf Trainingschromosomen).
  3. **Neuronales Modell:** eigenes Netz auf Pileup-Tensoren, trainiert auf GIAB + Q100 + Platinum Pedigree, Inferenz in reinem Rust (Burn oder Candle) über wgpu auf der GPU jedes Systems, CPU als Rückfall. Hier ist DeepVariant-Niveau das Ziel; dazu Merkmale aus Stufe 2 von `galign` (Pangenom-Kontext).
  - **Abkürzung, falls gewünscht:** Die veröffentlichten DeepVariant-Modelle (BSD-3) über eine eigene Rust/Core-ML-Inferenz betreiben. Das wäre schnell und genau, aber algorithmisch nicht eigen. Entscheidung für pj (Abschnitt 11).
- **Validierung:** hap.py/vcfeval gegen GIAB v4.2.1 und Q100, alle GIAB-Schichten, CMRG; Mendel an Platinum Pedigree und 1000G-Trios; 30x bis 10x; Illumina und MGI. Training und Test strikt getrennt (Proben und Chromosomen).
- **Ziele:** Stufe 1: SNV- und Indel-F1 mindestens gleich GATK HC auf allen Schichten, Laufzeit ein Bruchteil von HC. Stufe 3: Abstand zu DeepVariant beim SNV-F1 unter 0,05 Prozentpunkte, bei Indels unter 0,2 Prozentpunkte, dabei schneller als DeepVariant auf dem Mac.
- **Risiko:** Stufe 1 mittel, Stufe 3 hoch (Trainingsdaten, Rechenzeit für Training, Gefahr der Überanpassung an GIAB).
- **Veröffentlichung:** eigenes Repo mit DOI und Forschungsbericht; Stufe 3 als Preprint.

### 6.5 `gsv`: strukturelle Varianten
- **Ersetzt:** Delly 1.3.3, Manta 1.6 (Docker).
- **Ansatz:** zwei Wege kombiniert.
  1. **Genotypisierung bekannter SV** aus HPRC/HGSVC3 (k-mer- und Lesungsabgleich an Allel-Pfaden, Prinzip wie PanGenie/Paragraph). Ein großer Teil der SV eines Menschen ist in diesen Katalogen; hier schlagen Pangenom-Verfahren Kurzlesungs-Entdeckung deutlich.
  2. **Entdeckung neuer SV** aus Split-Reads, Discordant Pairs, Abdeckung und lokaler Assemblierung (Manta-Prinzip), mit Breakpoint-Verfeinerung.
- **Validierung:** Truvari gegen HG002 Q100-SV, HGSVC3, Platinum Pedigree; getrennt nach Typ und Größe; Vergleich mit Delly und Manta auf denselben Lesungen.
- **Ziele:** F1 für Deletionen über 50 bp deutlich über Delly und Manta; Insertionen (Schwäche der Kurzlesungs-Caller) über Genotypisierung mindestens doppelte Sensitivität; Laufzeit unter der von Manta.
- **Risiko:** mittel bis hoch.
- **Veröffentlichung:** eigenes Repo mit DOI, Forschungsbericht.

### 6.6 `gcnv`: Kopienzahl
- **Ersetzt:** CNVkit (und geplantes GATK gCNV).
- **Ansatz:** Lesetiefe in Fenstern, Korrektur für GC und Mappability, Normalisierung gegen ein Referenzpanel aus 1000G-30x (gleiche Plattform wählbar), HMM-Segmentierung mit Allelanteilen an heterozygoten SNPs (unterscheidet Deletion von Kopie-neutralem LOH). Paraloge Spezialgene gehen an `gloci`.
- **Validierung:** bekannte CNV aus HGSVC3 und Platinum Pedigree; Mendel in Trios; simulierte CNV.
- **Ziele:** höhere Sensitivität als CNVkit bei CNV über 10 kb bei gleicher Falsch-Positiv-Rate; Ergebnis mit Sicherheitsangabe je Segment.
- **Risiko:** mittel.
- **Veröffentlichung:** gemeinsam mit `gsv` (ein Repo "sv-cnv-wgs").

### 6.7 `gannot`: Annotation
- **Ersetzt:** Ensembl VEP (inkl. MANE-Proteinsequenzen für AlphaFold), dbSNP-rsID, ClinVar-Abgleich, gnomAD-Frequenzen.
- **Ansatz:** eigene Konsequenzvorhersage auf MANE Select/Plus Clinical, RefSeq und Ensembl-Transkripten: HGVS c./p./g. nach HGVS-Regeln (3'-Verschiebung, Duplikationen, Spleißregionen, Start-/Stop-Verlust, Frameshift mit neuer Stoppposition, Selenocystein, Mitochondrien-Code). Annotationsdaten als kompakte Bitpack-Tabellen (Prinzip echtvar), Abfrage mit Speicherabbild. Normalisierung (Links-Ausrichtung, Mehrallel-Zerlegung) eingebaut. Vorberechnete Scores (SpliceAI, AlphaMissense, REVEL, CADD) nur, wo die Lizenz es erlaubt.
- **Validierung:** alle ClinVar-Varianten mit Experten-HGVS (rund 2 Mio.): HGVS-Übereinstimmung mit ClinVar und mit VEP; Konsequenz-Konkordanz mit VEP; jeder Unterschied wird einzeln aufgelöst (wer hat recht, belegt mit Transkript). Zusätzlich Mutalyzer-Prüfung einer Stichprobe.
- **Ziele:** HGVS-Übereinstimmung mit ClinVar mindestens so hoch wie VEP; 5 Mio. Varianten (ein WGS) in unter 2 Minuten; Daten unter 6 GB; kein Absturz auf 1000G-Gesamt-VCF.
- **Risiko:** gering bis mittel (viele Sonderfälle, aber klare Wahrheit).
- **Gewinn für pj sofort:** beseitigt Segfaults, 27,6 GB Cache und Perl.
- **Veröffentlichung:** eigenes Repo mit DOI, Forschungsbericht mit dem ClinVar-Abgleich.

### 6.8 `ghla`: HLA (fertig, Einbau offen)
- **Ersetzt:** HLA*LA und T1K.
- **Stand:** pjlabshla v16 in Rust, 48 s, 1,3 GB; 184 ungesehene HPRC-Proben 2-Feld 2937/2941 (T1K 2921/2941), DRB3/4/5-Kopien fast fehlerfrei; Laborwahrheit 1000G läuft.
- **Offen:** Panel-Mapping und Stufe A in Rust (ein Programm, JSON), alle Gene inkl. DRA, E, F, G, H, K, V, GRCh37/no-ALT, Trio-Test, 4-Feld-Kalibrierung, Einbau in die App, Repo "hla-wgs" mit DOI.
- **Rolle im Plan:** Pilot für Einbau, Signierung, Schema und Schattenbetrieb. Was hier gelernt wird, gilt für alle folgenden Werkzeuge.

### 6.9 `gkir`: KIR und KIR-HLA-Liganden
- **Ersetzt:** T1K im KIR-Modus.
- **Ansatz:** gleiche Bauweise wie `ghla`: Lesungen der KIR-Region und ungemappte Lesungen sammeln, k-mer-EM über IPD-KIR, exakte Lesungs-Likelihood auf Kandidaten. Zusätzlich Gen-Kopienzahl und Haplotyp-Struktur (cA01/cB01 usw.), weil KIR-Gene fehlen oder doppelt vorkommen.
- **Validierung:** KIR-Wahrheit aus HPRC-Assemblies (Annotation mit IPD-KIR, wie bei HLA); Trios; Vergleich mit T1K.
- **Ziele:** Gen-Präsenz und Kopienzahl über 99 %, 3-stellige Allele klar über T1K.
- **Risiko:** mittel; der Wissenschaftsgewinn ist hier am größten, weil KIR aus WGS bisher kaum zuverlässig geht.
- **Veröffentlichung:** eigenes Repo "kir-wgs" oder gemeinsam mit hla-wgs; Paper-tauglich.

### 6.10 `gpgx`: Pharmakogenetik
- **Ersetzt:** Aldy (die CPIC/FDA-Regeln sind bereits eigen).
- **Ansatz:** Sternallele für alle CPIC-Gene aus PharmVar-Definitionen; für CYP2D6, CYP2A6, CYP2B6, UGT2B17, GSTM1, CYP4F2: Kopienzahl und Hybride (CYP2D6/2D7: *13, *36, *68) über paralog-spezifische Positionen (Cyrius-Prinzip) plus Lesungsphasierung; Unsicherheit als Liste möglicher Diplotypen; Phänotyp nach CPIC.
- **Validierung:** GeT-RM-Konsens (viele Proben in 1000G-30x), Wahrheit aus HPRC-Assemblies für CYP2D6-Struktur; Vergleich mit Aldy und PharmCAT auf denselben Proben.
- **Ziele:** CYP2D6-Diplotyp-Treffer über Aldy und Cyrius auf GeT-RM; alle anderen Gene 100 % auf GeT-RM-Konsens; unter 2 Minuten.
- **Risiko:** mittel.
- **Veröffentlichung:** eigenes Repo "pgx-wgs" mit DOI und Bericht.

### 6.11 `gloci`: Spezialgene mit Paralogen
- **Ersetzt:** KILDA (Lp(a) KIV-2); **neu** in 4.0: SMN1/SMN2 (SMA-Trägerstatus), GBA/GBAP1, CYP21A2, HBA1/HBA2, RHD, OPN1LW/MW, NCF1, PMS2.
- **Ansatz:** je Locus ein Modell aus paralog-spezifischen k-meren und Positionen, Kopienzahl relativ zu Kontrollregionen, Phasierung über Langlesungs-Haplotypen als Referenz (Prinzip Paraphase, aber für Kurzlesungen).
- **Validierung:** Wahrheit aus HPRC-Assemblies (KIV-2-Zahl, SMN-Kopien), Coriell-Proben mit bekannten SMA-Trägern, 1000G-Trios.
- **Ziele:** KIV-2-Zahl mit kleinerem Fehler als KILDA; SMN1-Kopienzahl über 99 %.
- **Risiko:** mittel; jedes Gen ist ein eigenes kleines Projekt.
- **Veröffentlichung:** Repo "loci-wgs", je Gen ein Abschnitt im Bericht.

### 6.12 `gstr`: Repeat-Expansionen
- **Ersetzt:** ExpansionHunter.
- **Ansatz:** Katalog der krankheitsrelevanten Repeats (STRchive), Lesungen an Graph aus Flanken und Repeat-Einheit ausrichten, Längenschätzung aus überspannenden, flankierenden und vollständig im Repeat liegenden Lesungen (auch aus ungemappten Paaren), Unsicherheitsintervall, Interruptionen (z. B. FMR1 AGG).
- **Validierung:** Repeat-Längen aus HPRC-Assemblies (exakte Wahrheit), Platinum Pedigree, Coriell-Proben mit bekannten Expansionen (HTT, FMR1, C9orf72, DMPK u. a.).
- **Ziele:** Längenfehler im Normalbereich kleiner als ExpansionHunter; Expansionen über Leselänge mit mindestens gleicher Sensitivität.
- **Risiko:** mittel.
- **Veröffentlichung:** eigenes Repo mit DOI.

### 6.13 `ghap`: Haplogruppen mtDNA und Y
- **Ersetzt:** Haplogrep (mt), bisherige Y-Ableitung.
- **Ansatz:** eigener mtDNA-Caller mit Heteroplasmie und NUMT-Filter, Klassifikation gegen PhyloTree/RSRS mit dem Haplogrep-Bewertungsprinzip; Y gegen ISOGG/YFull-Baum mit Bewertung fehlender und widersprüchlicher Marker.
- **Validierung:** 1000G-30x (veröffentlichte mt- und Y-Haplogruppen) mit Haplogrep3 als Vergleich; Heteroplasmie mit Mischproben.
- **Ziele:** Haplogruppen-Übereinstimmung mit Haplogrep3 über 99,5 %, Y bis auf die tiefste vom Baum gestützte Ebene; unter 10 s.
- **Risiko:** gering.

### 6.14 `gtel`: Telomere
- **Ersetzt:** TelomereHunter2.
- **Ansatz:** Zählung telomerischer Lesungen (TTAGGG und Varianten) mit GC-Korrektur und Normierung auf das Genom; Varianten-Repeats getrennt.
- **Validierung:** Korrelation mit TelomereHunter2 und TelSeq auf 1000G; bekannter Altersgang; wo möglich Telomerlängen aus T2T-artigen Assemblies.
- **Ziele:** gleiche Rangfolge wie TelomereHunter2 (Spearman über 0,98), Laufzeit läuft im Durchgang mit.
- **Risiko:** gering, Gewinn klein. Ehrlich: Telomerlänge aus Kurzlesungen bleibt eine grobe Schätzung.

### 6.15 `gmicro`: EBV und andere Viren
- **Ersetzt:** Kraken2, Bracken.
- **Ansatz:** kleiner, gezielter k-mer-Index (EBV, HHV-6 inkl. integrierter Form, HPV, HBV, HCV, HIV, Mykoplasmen als Kontaminationsanzeiger) statt einer Gesamtdatenbank; nur ungemappte Lesungen.
- **Validierung:** 1000G-Proben sind EBV-transformierte Zelllinien, also durchweg positiv; Proben aus Blut als Gegenprobe; künstlich beigemischte Viruslesungen.
- **Ziele:** gleiche EBV-Befunde wie Kraken2, Datenbank unter 500 MB, unter 1 Minute.
- **Risiko:** gering.

### 6.16 `gpgs`: polygene Scores
- **Ersetzt:** die Swift-Rechnung (behebt die 128-MiB-Grenze) und schließt zu pgsc_calc auf.
- **Ansatz:** Streaming-Abgleich harmonisierter PGS-Catalog-Dateien gegen gVCF mit Referenz-Konfidenz (bei WGS ist fast jede Position direkt genotypisiert, Imputation entfällt); Strang- und Allel-Abgleich; Ahnenkorrektur durch Projektion auf PCA der 1000G/HGDP-Referenz und Normierung wie pgsc_calc; Abdeckungsangabe je Score.
- **Validierung:** gleiche Rohscores wie pgsc_calc auf 1000G (Korrelation über 0,999); Prozentränge gegen pgsc_calc-Ahnenmodus.
- **Ziele:** auch Scores mit über 6 Mio. Varianten in unter 1 Minute; jede fehlende Variante begründet.
- **Risiko:** gering.

### 6.17 `glift`: Liftover
- **Ersetzt:** heutigen Ketten-Liftover.
- **Ansatz:** UCSC-Ketten plus Referenzallel-Prüfung und Behandlung vertauschter Allele (Prinzip bcftools +liftover).
- **Validierung:** dbSNP-Positionen in beiden Builds; Vergleich mit bcftools +liftover.
- **Ziele:** über 99,9 % Übereinstimmung, verworfene Varianten mit Grund.
- **Risiko:** gering.

### 6.18 SNP-Werkzeuge
SNP-Suche und Microarray-Export sind schon eigen; sie wandern auf `gc-io` (schneller in CRAM) und nutzen Genotypen aus `gcall` statt eigener Pileups. Proteinsequenzen für AlphaFold kommen aus `gannot`.

### 6.19 Workflow-Schicht
- **Ersetzt:** Nextflow, conda, Docker, Java.
- **Ansatz:** die App orchestriert wie heute über `OperationModel`; die Werkzeuge streamen untereinander über Pipes, wo es geht. brca12-wgs bekommt eine Version 2.0 ohne Nextflow: `galign` + `gcall` + `gsv` + `gannot` auf die BRCA-Regionen, gleicher Bericht, gleiches Schema `genome.brca12.v1`, GIAB-Benchmark erneut laufen lassen.

---

## 7. Reihenfolge

Sortiert nach Nutzen und Risiko. Aufwand grob in Arbeitswochen konzentrierter Arbeit (Erfahrungswert HLA: Prototyp bis eingefrorene Rust-Version etwa eine Woche intensiver Arbeit, ohne Einbau). Die Zahlen sind Schätzungen, keine Zusagen.

| Welle | Inhalt | Nutzen | Risiko | Aufwand |
|---|---|---|---|---|
| **0+1 FASTQ → BAM** (erster Schritt) | genome-core-Workspace, gc-io, genome-bench, `genome fq2bam` = FASTQ-Prüfung und Trimmen + `galign` Stufe 1 + `gbam` (Sortieren, Duplikate, Index) + BAM-QC, Einbau in 3.x im Schattenbetrieb, Repo mit DOI | sehr hoch (Grundlage jeder Analyse, längster Rechenschritt) | mittel | 9 bis 12 |
| **2 Schnelle Gewinne** | `gannot` (VEP weg), restliches `gqc` (Kontamination, Identität), `ghap`, `glift`, `gmicro`, `gtel`, `gpgs`, `ghla` vollständig + Einbau | hoch (Stabilität, 20 GB weniger, weniger conda) | gering | 7 bis 9 |
| **3 Wissenschaftlicher Vorsprung** | `gpgx`, `gkir`, `gstr`, `gloci` (Lp(a) zuerst, dann SMN u. a.) | sehr hoch (Bereiche, in denen die Welt heute schwache Werkzeuge hat) | mittel | 10 bis 14 |
| **4 Kernpipeline** | `gcall` Stufe 1+2, `gcnv`, `gsv`, brca12-wgs 2.0 | hoch (danach kein conda mehr) | mittel bis hoch | 12 bis 16 |
| **5 Spitzenklasse** | `gcall` Stufe 3 (neuronales Modell, GPU über wgpu) | hoch | hoch | 8 bis 12 |
| **6 Forschung** | `galign` Stufe 2 (pangenom-bewusst), SV-Entdeckung mit lokaler Assemblierung | wissenschaftlich hoch | hoch | offen |

**Die App selbst** läuft als eigener Strang parallel (Abschnitt 12.6): Gerüst, Designsystem und Genom-Karte ab Welle 2, Gleichstand mit 3.x bis Ende Welle 4. Aufwand zusätzlich grob 16 bis 24 Wochen; das ist der größte Einzelposten des ganzen Plans.

Warum diese Reihenfolge:
- FASTQ → BAM zuerst (Festlegung pj): jede Analyse hängt an der BAM, es ist der längste Rechenschritt, und der Benchmark-Rahmen, gc-io und der Einbauweg entstehen dabei gleich mit.
- Annotation vor Variant Calling, weil VEP heute am meisten Ärger macht und die Wahrheit (ClinVar-HGVS) klar ist.
- Spezialgene vor Variant Calling, weil dort der Abstand zum Stand der Technik am größten ist, wie bei HLA. Sie laufen dann schon auf BAMs aus `fq2bam`.
- Variant Calling spät, weil es das größte Risiko trägt und vom Fundament (gc-io, bench) abhängt.
- Ab Ende Welle 4 braucht Genome kein conda mehr. Ältere Werkzeuge bleiben bis dahin als Rückfall im Schattenbetrieb.

---

## 8. Einbau in die App

In der Rust-App (Abschnitt 12) entfallen Swift-Dekodierer und Prozess-Adapter. Bis die Rust-App den Stand von 3.x erreicht, bekommen die neuen Werkzeuge aber auch einen Weg in die Swift-App, damit pj sie sofort nutzt. Für diese Übergangszeit gilt, nach dem Muster von KILDA- und brca12-Adapter:

1. **Ein Adaptertyp `NativeToolAdapter`** statt eines Adapters je Werkzeug: Programmname, Version, Commit-Hash, Schema-Version, Eingaben, erwartete Ausgaben. Die Commit-Prüfung aus `ExternalToolRuntime` (Ref-Pinning) bleibt, nur gegen gebündelte Programme statt Downloads.
2. **Schema-Dekodierer in Swift** je `genome.<werkzeug>.vN`, Schema-Tests mit Beispiel-JSON aus dem Rust-Repo (wie `brca12_results.sample.json`).
3. **Schattenbetrieb:** Einstellung je Bereich "Klassisch / Genome-eigen / beide". Bei "beide" laufen beide Werkzeuge, die App zeigt eine Abgleichsansicht und schreibt Abweichungen ins Protokoll. Erst wenn das Ablösetor (Abschnitt 2) bestanden ist, wird "Genome-eigen" Standard; das alte Werkzeug fällt mit dem nächsten Release weg.
4. **Fortschritt:** JSON-Zeilen auf stderr → `advanceStep`, damit keine Phase still wirkt.
5. **Manifest:** jeder Lauf schreibt Programmversion, Commit, Datenpaket-Prüfsummen und Schema ins bestehende Manifest.
6. **Berichte:** neue Abschnitte (KIR-Struktur, SMN, Ahnen-PCA) über den bestehenden `ReportDocumentRenderer`; jede Karte nennt Methode und Validierungsstand mit Link zum Forschungsbericht.
7. **Tests:** Vertragstests je Adapter im isolierten Test-Host, Architekturregeln bleiben gültig (Services rufen Programme, Views nie).
8. **Werkzeugverwaltung:** der Bereich "Tools" zeigt für eigene Programme Version und Validierungsstand statt Installationsknöpfen; ToolFinder sucht nur noch für Rückfall-Werkzeuge.

---

## 9. Veröffentlichung

Wie bei brca12-wgs und geplant bei hla-wgs:

- **Eigenes öffentliches Repo je wissenschaftlich eigenständigem Werkzeug:** hla-wgs, kir-wgs (oder in hla-wgs), pgx-wgs, loci-wgs, str-wgs, sv-cnv-wgs, annot (Name offen), gcall (Name offen), gqc. Gemeinsame Crates als öffentliches `genome-core`. Der Quellcode der App bleibt privat.
- **Je Repo:** README, LICENSE (Vorschlag MIT oder Apache-2.0, wie die Bibliotheken, auf denen wir aufbauen), CITATION.cff, CHANGELOG, CI, Release-Tag, Zenodo-DOI pro Release, Benchmark-Skripte zum Nachrechnen jeder Zahl.
- **Forschungsbericht je Werkzeug** auf pjlabs.dev im Forschungsbereich: Methode, Wahrheitssätze, Testdesign mit eingefrorenen Versionen, Ergebnisse mit Intervallen, Grenzen. Für HLA/KIR, PGx und den Variant Caller zusätzlich Preprint auf bioRxiv.
- **Lizenzen der Daten:** ClinVar, gnomAD, GIAB, HPRC, 1000G, PharmVar, IMGT/HLA, IPD-KIR, PGS Catalog sind frei nutzbar (mit Zitierpflicht). CADD, SpliceAI-Scores und AlphaMissense haben Einschränkungen (nicht kommerziell bzw. CC BY-NC-SA); ob Genome sie mitliefert, entscheiden wir je Datei.
- **Hinweis überall:** Forschungswerkzeug, kein Medizinprodukt (RUO).

---

## 10. Risiken

| Risiko | Gegenmaßnahme |
|---|---|
| Überanpassung an GIAB/HPRC (Werkzeug lernt die Testproben) | strikte Trennung, Einfrieren vor jeder Runde, Reserveproben, Platinum Pedigree als unabhängige Prüfung |
| Variant Caller erreicht DeepVariant nicht | Stufe 1+2 sind für sich schon ein Gewinn gegenüber bcftools; Entscheidung über DeepVariant-Modell-Abkürzung offen halten |
| Arbeitsumfang (ein Entwickler) | Wellen sind einzeln nutzbar; jede Welle endet mit einem Release, nichts hängt halbfertig |
| ExtremeSSD trennt sich unter Last | Arbeitsdaten auf interne SSD, Prüfsummen, atomare Veröffentlichung |
| Container-Netz (EBI, NCBI gesperrt) | diese Daten über den Mac holen und in genome-reflib spiegeln |
| Datenlizenzen | Liste je Datenpaket, Nachweis im Manifest |
| MGI-Daten verhalten sich anders als Illumina | MGI-Proben in jede Testrunde (1000G hat keine; pj's Daten und öffentliche DNBSEQ-GIAB-Läufe) |

---

## 11. Entscheidungen für dich

1. **Absicht freigeben:** Genome 4.0 nach diesem Plan, Start mit FASTQ → BAM (Abschnitt 0).
2. **Aligner:** reine Rust-Umsetzung des BWA-MEM-Verfahrens (Empfehlung, passt zu "alles in Rust", messbar gleich genau) oder gleich ein neues Verfahren (Minimizer wie minimap2 oder Pangenom, höheres Risiko). Ein neues Verfahren kann später als Stufe 2 dazukommen.
3. **Variant Calling Stufe 3:** eigenes neuronales Modell (eigenständig, länger, Paper-tauglich) oder veröffentlichte DeepVariant-Modelle über eigene Core-ML-Inferenz (schneller am Ziel, aber nicht algorithmisch eigen). Meine Empfehlung: zuerst Stufe 1+2 bauen, die Frage in Welle 4 mit Messwerten entscheiden.
4. **Repos:** ein Repo je Werkzeug (wie brca12-wgs und hla-wgs, je eigene DOI) oder ein großes öffentliches Repo mit mehreren DOIs. Meine Empfehlung: je Werkzeug, gemeinsame Crates in `genome-core`.
5. **Oberflächen-Bibliothek:** GPUI oder Iced, entschieden nach einem zweiwöchigen Vergleichsprototyp (Abschnitt 12.2). Meine Vorab-Empfehlung: GPUI.
6. **Rechenort für Benchmarks:** Container für Breitentests, Mac für Tempo und eigene Proben (Empfehlung) oder alles auf dem Mac.

**Erster Schritt nach Freigabe (Woche 1):** heutigen FASTQ→BAM-Weg der App auslesen (Werkzeuge, Versionen, Parameter, Referenz), Basismessung auf pj's Mac mit HG002 30x und pj's eigenen Lesungen (BWA-MEM2 + samtools sort/markdup/index: Wanduhr, RAM, Platte), genome-core-Workspace mit gc-io anlegen (Meilenstein M1 in Abschnitt 0).

---

## 12. Die App in Rust

### 12.1 Was das bedeutet, ehrlich

Genome 3.0 hat 1.134 Swift-Dateien und 655 Testdateien. Die App in Rust neu zu bauen ist der größte Einzelposten dieses Plans, größer als jedes einzelne Werkzeug. Der Gewinn:

- **Ein Code, eine Sprache:** Analysekern, Zustand, Berichte und Oberfläche teilen Typen. Keine Swift-Dekodierer für JSON-Schemata, keine doppelte Logik (heute rechnet z. B. der PGS-Abgleich in Swift anders als ein BAM-Pfad).
- **Tempo in der Oberfläche:** GPU-gerenderte Ansichten über Millionen Lesungen und Varianten, was SwiftUI-Listen nicht schaffen.
- **Prüfbare Schichtung:** die Regel "Views → Models → Services" wird vom Crate-Graphen erzwungen (eine UI-Crate kann die I/O-Crate gar nicht importieren). Das ersetzt die brüchigen String-Grep-Architekturtests.
- **Linux und Windows** aus demselben Code (Festlegung pj, Abschnitt 12.3).

Was verloren geht und ersetzt werden muss: SwiftUI-Systemsteuerelemente, WKWebView für HTML→PDF, Apple String Catalogs, PJDesignKit, Sparkle-artige Updates, Teile der Barrierefreiheit. Dafür steht unten je ein Ersatz.

### 12.2 Technik

| Bereich | heute (Swift) | Genome 4.0 (Rust) |
|---|---|---|
| Oberfläche | SwiftUI | **GPUI** (GPU-gerendert, aus dem Zed-Editor, läuft auf macOS, Linux und Windows, sehr schnell bei großen Datenmengen) oder **Iced** (wgpu, alle drei Systeme, stabile API, Elm-Architektur). Pflicht im Vergleich: gleiche Qualität auf allen drei Systemen. Entscheidung nach Prototyp: dieselbe Genom-Karte mit pj's BAM in beiden, gemessen an Bildrate, Speicher, Textqualität mit IBM Plex Mono und Barrierefreiheit |
| Genom-Ansicht | IGV extern / keine | eigene GPU-Zeichenfläche (wgpu: Metal, Vulkan, DirectX 12) für Karte, Abdeckung, Lesungen |
| Berichte, PDF | HTML + CSS → WKWebView → PDF | **Typst als Bibliothek** (Rust, Apache-2.0): Satz in Buchqualität, eingebettete IBM Plex Mono, Inhaltsverzeichnis, Quellenkarten, APA-Literatur; dieselbe Engine zeigt den Bericht in der App |
| Lokalisierung | String Catalogs | **Fluent** (Mozilla, Rust): Deutsch und Englisch, grammatikalisch korrekte Mehrzahl |
| Designsystem | PJDesignKit | Crate `pj-design`: Tokens (Terracotta `#D97757`, 14 px Radien, 36 px Aktionsflächen, Plex Mono 400/500, keine Fettschrift) als Rust-Konstanten, dieselben Werte exportiert für Typst und pjlabs.dev |
| Zustand | AppState + @Observable-Modelle | Crate `genome-engine`: Ereignis-gesteuerter Zustand, Läufe als dauerhafte Aufträge (neustartfest, wie im GPT-6-Auftrag gefordert) |
| Speicherung | Dateien + JSON | SQLite (rusqlite) für Projekte, Läufe, Befunde, Manifeste; Ergebnisdateien bleiben Dateien |
| MCP-Anbindung | eigene Swift-Implementierung | offizielles Rust-SDK `rmcp`, standardmäßig aus wie heute |
| Updates | eigene Update-Prüfung | eigener Updater mit Ed25519-Signaturen, Hosts in fester Allowlist |
| Tests | XCTest, SwiftSyntax-Regeln | `cargo test`, Snapshot-Tests der Oberfläche und der Typst-Berichte, Crate-Graph-Regeln mit `cargo-deny` |
| Wissensbasen | CSV/JSON im Bundle | unverändert übernommen (KnowledgeBaseCSV, HaplogroupContextCSV), Importskripte nach Rust |

### 12.3 macOS, Linux und Windows

Festlegung pj (18.11.): Genome 4.0 muss so gebaut sein, dass es künftig auch auf Linux und Windows läuft. Regeln ab dem ersten Commit:

- **CI auf allen drei Systemen** (macOS arm64, Linux x86-64 und arm64, Windows x86-64): Bauen, Tests, Lint bei jedem Push. Ein Commit, der auf einem System bricht, gilt als kaputt.
- **Nur reines Rust oder plattformneutrale Crates.** Keine Systemaufrufe direkt im Fachcode; alles Systemnahe (Pfade, Prozesse, Speicherabbild, Dateisperren, freier Speicher, RAM-Größe, Schlaf verhindern während langer Läufe) liegt in einer Crate `gc-sys` mit einer Umsetzung je System.
- **Pfade:** `PathBuf`, nie Zeichenketten mit `/`; lange Windows-Pfade (`\\?\`) und Nicht-ASCII-Namen (iCloud-Pfade, Umlaute) in Tests.
- **Prozesse:** Arbeitsprozesse über `std::process`, Abbruch plattformgerecht (Signale auf Unix, Job Objects auf Windows). Kein `fork`.
- **Dateiformate:** Zeilenende `\n` beim Schreiben, `\r\n` beim Lesen toleriert (FASTQ-Leser macht das schon).
- **SIMD:** gemeinsame Abstraktion mit Umsetzungen für NEON und AVX2/AVX-512, Auswahl zur Laufzeit; Ergebnis auf allen Systemen bitgleich (Test in CI).
- **GPU:** nur über wgpu (Metal, Vulkan, DirectX 12); keine Core-ML- oder Metal-only-Wege.
- **Benchmarks** laufen auf pj's Mac und in Linux (Container, CI); ein Windows-Lauf je Meilenstein.
- **Datenordner:** je System der übliche Ort (`~/Library/Application Support/Genome`, `$XDG_DATA_HOME/genome`, `%APPDATA%\Genome`), frei wählbar für externe Platten.

Das kostet anfangs etwa 10 % mehr Aufwand, verhindert aber, dass später ein Umbau nötig wird.

### 12.4 Crate-Schichten

```
pj-design          Designtokens
gc-*               Analysekern (Abschnitt 4)
genome-kb          Wissensbasen, Interpretation, ACMG-SF, CPIC-Regeln
genome-engine      Projekte, Läufe, Auftragswarteschlange, Manifest, SQLite
genome-report      Typst-Vorlagen, Befundkarten, Literatur
genome-ui          Oberfläche (darf nur genome-engine und genome-report importieren)
genome             ein Programm: Oberfläche, `genome worker …`, `genome mcp`
```

### 12.5 Was mitgenommen wird

Nicht neu erfunden, sondern übernommen: die Wissensbasen und ihre Belege, die Berichtsstruktur und Texte (über Fluent), die Fachregeln aus der Swift-App (CPIC, ACMG-SF, KIR-HLA-Liganden, PGS-Haploid-Regel vom 01.10., Probenbindung aus gvali R1 bis R12), die QMS-Struktur (Traceability, Risiken, SOUP-Liste wird kürzer), alle gelösten Fehlerfälle aus den Bug-Hunts als Rust-Tests.

### 12.6 Umstieg ohne Stillstand

1. Genome 3.x bleibt nutzbar und bekommt nur noch Fehlerbehebungen und die neuen Werkzeuge über Adapter (Abschnitt 8).
2. Die Rust-App wächst Bereich für Bereich: zuerst Projekte, Läufe, Werkzeugkern und Genom-Karte; dann Befunde und Berichte; dann Wissensbasen, PGS, SNP-Werkzeuge, Layout-Editor, MCP.
3. Für jeden Bereich gibt es einen Gleichstandstest: dieselbe Probe in 3.x und 4.0, Befunde und Bericht verglichen.
4. Projekte aus 3.x öffnet 4.0 direkt (Import der Ordner und Manifeste).
5. Wenn alle Bereiche den Gleichstand haben, wird 4.0 die Haupt-App und 3.x eingefroren.

---

## 13. UX-Konzept: das Genom als Ort, nicht als Formular

Leitidee: Heute ist Genome eine Werkzeugsammlung mit Knöpfen pro Analyse. Genome 4.0 zeigt **das eigene Genom** als Ort, an dem man sich bewegt, und jede Aussage lässt sich bis zur einzelnen Lesung zurückverfolgen. Die Gestaltung bleibt im PJ-Labs-Kanon: IBM Plex Mono, Terracotta nur für die eine wichtigste Handlung, ruhige Linien, viel Weißraum, keine Fettschrift.

### 13.1 Die Genom-Karte (Startfläche)

- Eine einzige zoombare Fläche: 24 Chromosomen als ruhige Linien → Region → Gen mit Exons → Abdeckung → einzelne Lesungen mit Basen. Stufenloser Zoom, GPU-gerendert, flüssig auch bei 30x.
- Befunde sind Punkte auf der Karte (HLA auf 6p, APOE auf 19, CYP2D6 auf 22, BRCA auf 13/17, Lp(a) auf 6q). Farbe ist nie Bewertung, sondern Kategorie; die Sicherheit zeigt ein Ring (voll = sicher, teilweise = Alternativen vorhanden).
- Die Karte ersetzt den externen IGV-Aufruf: bei maximalem Zoom ist sie ein vollwertiger Lesungsbetrachter, der auch die Haplotypen (Phase) farblich trennt.
- Mitochondrium und Y als eigene kleine Kreise mit Haplogruppen-Pfad im Baum.

### 13.2 Befehlszeile statt Menüs

`⌘K` öffnet eine Zeile, die alles versteht, was ein Genom-Mensch tippt: `rs429358`, `APOE`, `chr6:29.9M`, `CYP2D6`, `HLA-B`, `Lp(a)`, `Bericht`, `HG002 vergleichen`. Treffer springen auf der Karte hin und öffnen die Befundkarte. Alles lokal; ein externer KI-Assistent kann über MCP dieselben Befehle nutzen, wenn pj das einschaltet.

### 13.3 Befundkarten mit Beweiskette

Jede Befundkarte hat drei Ebenen, aufklappbar von oben nach unten:

1. **Aussage** in einem Satz, mit Sicherheitsring ("CYP2D6 *1/*4, intermediärer Metabolisierer").
2. **Warum:** die Merkmale, die zur Aussage führen (Allele, Kopienzahl, Alternativen mit Abstand), und die Leitlinie (CPIC, Quelle).
3. **Beweis:** die Lesungen an den entscheidenden Positionen direkt in der Karte, Werkzeug und Version, Datenstand, und der **Validierungsstand des Werkzeugs** ("auf 184 ungesehenen Proben 2937/2941 richtig, Bericht ↗").

Damit wird die Wissenschaftlichkeit, die die Werkzeuge haben, auch in der Oberfläche sichtbar. Ein Schalter je Karte wechselt zwischen verständlicher Sprache und Fachsprache.

### 13.4 Ergebnisse, die während der Analyse entstehen

- Der Lauf ist ein ruhiger Fluss-Graph (FASTQ → Alignment → … → Bericht), jeder Knoten zeigt Fortschritt, Tempo und geschätzte Restzeit aus echten Ereignissen.
- Weil die Werkzeuge streamen, erscheinen Befunde, sobald ihre Daten da sind: Geschlecht, Kontamination und mt-Haplogruppe nach wenigen Minuten; HLA, KIR und PGx aus einer frühen Zielregion-Vorauswahl, noch bevor das ganze Genom aligniert ist. Vorläufige Befunde sind als vorläufig markiert und werden still bestätigt oder korrigiert.
- Unterbrechung, Neustart, Abbruch: der Graph zeigt genau, welche Schritte fertig und gültig sind (neustartfest).

### 13.5 Vergleichen

- **Zwei Proben** nebeneinander auf derselben Karte (z. B. pj und HG002), Unterschiede hervorgehoben.
- **Familie/Trio:** Vererbung je Variante sichtbar, Mendel-Fehler markiert, Phasierung aus den Eltern.
- **Zwei Werkzeuge** (Schattenbetrieb aus Abschnitt 8) oder zwei Datenstände: Abweichungen als Liste mit Sprung zur Stelle.

### 13.6 "Was ist neu?"

Weil die Annotation in Minuten läuft, prüft Genome bei jedem neuen ClinVar-, PharmVar- oder PGS-Catalog-Stand das Genom neu und zeigt eine kurze Liste der geänderten Befunde mit Begründung ("ClinVar stuft Variante X seit 2027-01 als wahrscheinlich gutartig ein"). Kein stilles Überschreiben, jede Änderung bleibt in der Zeitachse.

### 13.7 Bericht als lebendes Dokument

Der Bericht ist eine Ansicht in der App, gesetzt mit derselben Typst-Engine wie das PDF: was man sieht, ist das PDF. Layout-Änderungen (heute der Layout-Editor) wirken sofort. Jede Karte im Bericht verlinkt zurück auf die Stelle in der Karte.

### 13.8 Ruhe, Privatsphäre, Zugang

- Eine dauerhafte, kleine Anzeige "lokal, offline", die zeigt, ob überhaupt eine Netzverbindung genutzt wird und wofür (nur Datenpakete, nur auf erlaubte Hosts).
- Tastatur zuerst, jede Handlung auch ohne Maus; Bildschirmleser-Unterstützung über AccessKit (Rust; VoiceOver, NVDA/Narrator, Orca) als Pflicht im Oberflächen-Vergleich.
- Hell und dunkel, Kontrast nach WCAG AA, Farbsehschwäche geprüft (Kategorien nie nur über Farbe).

### 13.9 Wie das UX geprüft wird

- Klickbarer Prototyp der Karte und der Befundkarte in Welle 0 mit pj's echter Probe, bevor die Oberfläche gebaut wird.
- Aufgabentests mit 3 bis 5 Menschen (z. B. "Finde deinen CYP2D6-Befund und zeige, warum er so lautet"), Zeit und Fehler gemessen, vor und nach jeder großen Änderung.
- Tempo-Ziele der Oberfläche: Start unter 1 s, Zoom flüssig mit 120 Bildern/s auf ProMotion, jede Suche unter 100 ms.
