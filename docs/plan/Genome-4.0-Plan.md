# Genome 4.0: eigene Rust-Werkzeuge für jeden Analyseschritt

Stand 18.11.2026, von pj gestartet (Ergänzungen: ganze App in Rust, neues UX; erster Schritt FASTQ → BAM, Abschnitt 0; läuft künftig auch auf Linux und Windows, Abschnitt 12.3; ein eigener Bericht je Werkzeug und daraus der Gesamtbericht, Abschnitt 14; Protokollierung mit voller Transparenz, Abschnitt 15).

**Stand der Umsetzung:** Das private Repo `genome4` steht mit Grundgerüst, FASTQ-Leser und BAM-Rundreise (Meilenstein M1). Die Basismessung der heutigen Kette läuft in der Nacht 18./19.11. auf pj's Mac. Der Plan wurde am 18.11. von vorne bis hinten durchgesehen und berichtigt (Liste am Ende).

## Kurzfassung

- Genome 3.0 ruft heute rund 30 fremde Werkzeuge über conda, Docker, Java und Nextflow auf. Genau dort entstehen die meisten Ausfälle (VEP-Segfaults, Manta nur über Docker, GATK-Versionsabfrage ohne Ende, CNVkit-SIGBUS von der SSD, 27,6 GB VEP-Cache).
- Genome 4.0 ersetzt sie durch **eine eigene Rust-Werkzeugfamilie `genome-core`**: signierte Programme, die direkt in der App liegen, ohne conda, ohne Docker, ohne Java, ohne Nextflow. Jedes Programm liefert JSON nach festem Schema, wie hla-wgs und brca12-wgs.
- Maßstab ist das HLA-Vorgehen: Wahrheitssätze aus Assemblies und Laborergebnissen, Version vor jeder Testrunde einfrieren, ungesehene Proben, Vergleich mit dem Original auf denselben Lesungen. Ein Werkzeug ersetzt das alte erst, wenn es dieses Tor besteht.
- **Ehrlich:** In acht Bereichen ist "besser als das Original" realistisch und messbar (Annotation, QC, PGx, KIR, Repeats, Lp(a), Haplogruppen, Spezialgene). Bei SV/CNV ist es mit Pangenom-Genotypisierung erreichbar, aber aufwendig. Beim Variant Calling ist GATK-Niveau gut erreichbar, DeepVariant-Niveau nur mit eigenem neuronalen Modell und viel Arbeit. Beim Aligner ist ein neues Verfahren, das BWA-MEM an Genauigkeit schlägt, nicht realistisch; wir setzen das bewährte BWA-MEM-Verfahren in reinem Rust neu um (gleich genau, schneller durch Verschmelzen der Schritte) und forschen erst später an einem Pangenom-Aligner.
- "Komplett unabhängig" heißt: keine fremden Programme mehr installieren. Referenzdaten (GRCh38, ClinVar, gnomAD, PharmVar, IMGT/HLA, PGS Catalog) bleiben fremde Daten. Die spiegeln wir mit Prüfsummen (wie genome-plus).
- **Alles in Rust, auch die App:** Genome 4.0 ist kein Swift-Programm mehr, sondern eine native Rust-App (GPU-gerenderte Oberfläche, Berichte mit Typst, Lokalisierung mit Fluent). Analysekern und Oberfläche teilen dieselben Crates. Abschnitt 12.
- **Neues UX:** eine zoombare Genom-Karte als Startfläche (vom Chromosom bis zur einzelnen Lesung), Befunde mit sichtbarer Beweiskette, Ergebnisse, die schon während der Analyse erscheinen, und ein Bericht, der in der App genau so aussieht wie im PDF. Abschnitt 13.
- **Modulare Berichte (Festlegung pj):** jedes Werkzeug (HLA, BRCA, KIR, Lp(a), PGx und alle anderen) erstellt einen eigenständigen Bericht; der Gesamtbericht wird nur aus diesen Bausteinen zusammengesetzt. Abschnitt 14.
- **Protokollierung (Festlegung pj):** jeder Lauf hält lückenlos fest, was er getan hat; Fehler erklären sich selbst (was, wo, warum, was tun). Abschnitt 15.
- **Erster Schritt (Festlegung pj): FASTQ → BAM komplett in Rust**, ein Programm `genome fq2bam` vom FASTQ bis zur indizierten, duplikatmarkierten BAM mit QC. Ausgearbeitet in Abschnitt 0.
- Danach: schnelle sichere Gewinne (Annotation, QC-Rest, Haplogruppen, HLA-Einbau), Spezialgene mit wissenschaftlichem Vorsprung (PGx, KIR, Repeats, Lp(a)), Variant Calling und SV/CNV, zuletzt Forschung (neuronaler Caller, Pangenom-Alignment).

---

## 0. Erster Schritt: FASTQ → BAM in Rust

Festlegung pj (18.11.): Genome 4.0 beginnt damit, wie aus einer FASTQ eine BAM wird. Ziel ist ein einziges Programm `genome fq2bam`, das ohne fremde Werkzeuge aus gepaarten FASTQ-Dateien eine sortierte, duplikatmarkierte, indizierte BAM (oder CRAM) mit vollständiger Qualitätsauswertung macht, gleich genau wie die heutige Kette (bwa 0.7.19 + samtools) und deutlich schneller.

### 0.1 Heute und morgen

| Schritt | heute (3.x) | `genome fq2bam` |
|---|---|---|
| FASTQ prüfen, trimmen | fastp als eigener Lauf, getrimmte Zwischen-FASTQ auf Platte (FastQC ist eine getrennte Aktion) | im Lesestrom, keine Zwischendatei |
| Alignment | bwa 0.7.19 (`bwa mem`), SAM als Text durch eine Pipe | BWA-MEM-Verfahren in Rust, Datensätze direkt im Speicher |
| Sortieren | samtools sort in derselben Pipe (3 Threads × 448 MB, Temp-Dateien) | im selben Prozess, Zwischenstände komprimiert und mit Prüfsumme im Temp-Ordner |
| Duplikate | samtools fixmate vor und markdup nach dem Sortieren, in derselben Pipe | beim Zusammenführen der sortierten Teile |
| Index, Kennzahlen | samtools index und flagstat (zwei weitere Lesedurchgänge); mosdepth ist eine getrennte Aktion | im selben Schreibdurchgang |
| Ergebnis | BAM + verstreute Logdateien | BAM/CRAM + Index + `fq2bam.json` (Schema `genome.fq2bam.v1`) mit Kennzahlen, Versionen, Prüfsummen + eigener Bericht (Abschnitt 14) |

Heute werden die Daten viermal komplett gelesen und zweimal komplett geschrieben (fastp liest und schreibt alles, die Pipe liest und schreibt alles, Index und flagstat lesen die BAM je noch einmal); dazu kommen die Temp-Dateien der Sortierung. `fq2bam` liest das FASTQ einmal und schreibt die BAM einmal.

### 0.2 Die Stufen im Programm

Alle Stufen laufen gleichzeitig als Datenstrom über alle Kerne:

1. **Einlesen:** FASTQ und FASTQ.gz (gepaart, verschachtelt, mehrere Lanes), auch unalignierte oder bereits alignierte BAM/CRAM als Eingang (zum Neu-Alignieren, z. B. pj's PJ.bam). Read Groups aus dem Kopf der Lesungen, für Illumina- und MGI/DNBSEQ-Schreibweise.
2. **Prüfen, bevor gerechnet wird:** Paare passen zusammen, Qualitätskodierung, Längen, abgeschnittene Dateien. Fehler brechen sofort mit klarer Meldung ab, nicht nach Stunden.
3. **Trimmen:** Adapter über die Überlappung der Paare (fastp-Prinzip) plus bekannte Illumina- und MGI-Adapter; Poly-G-Enden (zweifarbige Geräte). Kein hartes Qualitätstrimmen als Standard, weil das Alignment ohnehin weich clippt; die Wirkung beider Varianten wird gemessen und im Bericht begründet.
4. **Referenz und Index:** GRCh38 Analyse-Set mit Decoys (hs38d1, wie pj's Daten), mit oder ohne ALT. Suffix-Array (SA-IS) und FM-Index in Rust gebaut, als fertiges Datenpaket aus genome-reflib geladen und per Speicherabbild genutzt. Zwei Varianten: **schnell** (BWA-MEM2-artige Belegungstabellen, mehr RAM) und **kompakt** (unter 6 GB, für Macs mit 16 GB); Genome wählt nach Speicher.
5. **Alignment (BWA-MEM-Verfahren):** SMEM-Seeds mit Reseeding, Chaining und Chain-Filter, gebänderte Smith-Waterman-Erweiterung mit SIMD (NEON auf ARM, AVX2/AVX-512 auf x86, zur Laufzeit gewählt) über viele Lesungen gleichzeitig, Insert-Größen je Block geschätzt, Mate-Rescue, Paarbewertung, MAPQ, supplementäre Alignments, alternative Treffer (XA), ALT-bewusste Behandlung. **Ergebnis unabhängig von der Kernzahl** (feste Blockgröße, wie BWA `-K`), damit Läufe reproduzierbar sind.
6. **Sortieren:** Datensätze gehen ohne SAM-Text in Positionsfächer im Speicher; ist das RAM-Budget voll, werden sortierte, komprimierte Teile mit Prüfsumme in den Temp-Ordner geschrieben (intern, wenn Platz ist, sonst extern) und am Ende zusammengeführt; fällt die Platte aus, setzt der Lauf an den geprüften Teilen wieder an.
7. **Duplikate:** beim Zusammenführen markiert (Position, Orientierung, Mate, Bibliothek), optische Duplikate über Kachel-Koordinaten aus dem Lesungsnamen (Illumina und MGI).
8. **Schreiben:** BAM mit paralleler BGZF-Kompression (zlib-rs, reines Rust) oder CRAM 3.1; BAI/CSI/CRAI im selben Durchgang; MD/NM-Tags; `@PG` mit exakter Version und allen Parametern; Prüfsumme.
9. **Kennzahlen unterwegs:** Anteil aligniert und richtig gepaart, Duplikatrate, Insert-Größen, MAPQ-Verteilung, Abdeckung je Chromosom und Region (mosdepth-gleich), GC-Verzerrung, Geschlecht aus X/Y-Abdeckung, Fehlerrate je Zyklus.

Fortschritt geht als Ereignisse an die App (Stufe, Lesungen pro Sekunde, Restzeit), Abbruch ist jederzeit sauber, ein Neustart setzt an fertigen Teilen wieder an.

### 0.3 Validierung

Vergleichsmaßstab ist `bwa mem` 0.7.19, das Genome 3.x heute nutzt. BWA-MEM2 liefert nach Angabe seiner Autoren dieselben Alignments und dient als Tempo-Vergleich; sein Index braucht deutlich mehr RAM, deshalb läuft es voraussichtlich auf Linux und nicht auf pj's 16-GB-Mac.

Die Frage ist nicht "gleich wie bwa", sondern "landet jede Lesung dort, wo sie herkommt, und werden die Varianten danach mindestens so gut gefunden".

1. **Simulierte Lesungen mit bekannter Herkunft:** aus der diploiden HG002-Q100-Assembly und HPRC-Assemblies, mit Fehlerprofilen, die aus echten Illumina- und MGI-Läufen gelernt sind. Jede Lesung kennt ihre wahre Position, also wird die Fehlplatzierungsrate exakt messbar, getrennt nach GIAB-Schichten (Segmentduplikationen, Homopolymere, MHC, KIR, CYP2D6, LPA). Zusätzlich die **MAPQ-Kalibrierung**: stimmt "MAPQ 30" mit einer tatsächlichen Fehlerrate von 1:1000 überein. Dieselben Lesungen gehen durch bwa, BWA-MEM2 und minimap2.
2. **Echte Daten gegen bwa + samtools (die heutige Kette):** GIAB HG001 bis HG007 (Illumina 30x), HG002 auf MGI/DNBSEQ, 1000G-30x-Proben, pj's eigene Lesungen. Gemessen: Anteil gleicher Primärpositionen, jede Abweichungsgruppe nach Region erklärt, Duplikat-Übereinstimmung, flagstat-Zahlen.
3. **Wirkung auf die Varianten (der eigentliche Richter):** auf unserer BAM und auf der bwa-BAM laufen dieselben eingefrorenen Caller (DeepVariant und GATK HC als neutrale Richter). hap.py gegen GIAB v4.2.1, HG002 Q100 und CMRG, alle Schichten.
4. **Reproduzierbarkeit:** byteidentische BAM bei 1, 4 und 10 Kernen und bei Wiederholung auf demselben System; zwischen macOS, Linux und Windows identische Datensätze (die komprimierten Bytes können je nach Prozessor abweichen, der Inhalt nicht).
5. **Robustheit:** abgeschnittene Dateien, vertauschte Paare, gemischte Leselängen (2×100 bis 2×150), 10x bis 60x, BAM als Eingang.
6. **Testdesign wie bei HLA:** Version vor der Testrunde einfrieren; Entwicklung auf chr20 und HG001, Test auf den übrigen Proben.

### 0.4 Messziele

| Kennzahl | Ziel |
|---|---|
| Fehlplatzierung (simuliert, MAPQ ≥ 20) | höchstens so hoch wie bwa, in jeder Schicht |
| MAPQ-Kalibrierung | mindestens so genau wie bwa, besser angestrebt (gemessener Fehler je MAPQ-Stufe näher am Sollwert) |
| Gleiche Primärposition wie bwa (echt) | über 99,5 %, alle Abweichungen erklärt |
| SNV-/Indel-F1 danach (DeepVariant und HC als Richter) | nicht schlechter als mit bwa-BAM (95-%-Intervall), in schweren Regionen besser angestrebt |
| Duplikate | Übereinstimmung mit samtools markdup über 99,9 % |
| Wanduhr FASTQ bis indizierte BAM, 30x, pj's Mac | höchstens die Zeit, die `bwa mem` in der heutigen Kette allein braucht (Trimmen, Sortieren, Duplikate, Index, QC kosten nichts extra); Streckziel 1,5-mal schneller als die heutige Gesamtkette (11 h 40 min). Feste Zahlen folgen aus der Basismessung |
| Spitzen-RAM | unter 12 GB mit kompaktem Index (der Mac hat 16 GB, System und App brauchen den Rest) |
| Zwischenspeicher auf Platte | höchstens 1,2-mal die fertige BAM |

**Basis, ausgelesen am 18.11. (siehe `M1-Bestand-Genome3-FASTQ-BAM.md`):** Genome 3.x nutzt bwa 0.7.19 (BWA-MEM2 ist nicht installiert), fastp als eigenen Lauf, dann bwa mem | samtools fixmate | sort (3 Threads × 448 MB) | markdup, danach index und flagstat. Einzige belegte Laufzeit: HG001 30x in **11 h 40 min** auf dem M5 (10 Kerne, 16 GB). Zeiten je Schritt misst der Lauf in der Nacht 18./19.11. (ganzes HG001, Skript `bench/run_baseline.sh`, außerhalb der App, gleiche Befehle); die Zahlen kommen danach hierher. Die Ziele oben gelten gegen diese Kette; BWA-MEM2 ist zusätzlicher Tempo-Vergleich.

**Folge aus dem Bestand:** Die interne SSD hat nur 27 GiB frei. Zwischenstände der Sortierung müssen deshalb auf die externe Platte passen und einen Plattenausfall überstehen (Prüfsumme je Teil, Wiederaufnahme). Der kompakte Index ist für 16 GB RAM Pflicht.

### 0.5 Meilensteine

| M | Inhalt | Prüfung (geht weiter nur, wenn bestanden) | Dauer |
|---|---|---|---|
| M1 | Workspace, gc-io: FASTQ/BAM lesen und schreiben, BGZF parallel; Basismessung der heutigen Kette | Rundreise identisch zu samtools. **Stand 18.11.: erledigt** (844.606 Datensätze identisch, CI grün auf vier Systemen); Basismessung läuft | 1 Woche |
| M2 | SA-IS, FM-Index, beide Varianten; Protokoll-Grundlage `gc-log` | exakte Trefferzahlen und Positionen bei 1 Mio. Abfragen, geprüft gegen eine einfache Vollsuche und gegen bwa | 1 bis 2 Wochen |
| M3 | Seeds, Chaining, SIMD-Erweiterung, Einzel-Lesungen auf chr20 | über 99 % gleiche Positionen wie bwa, Tempo auf Kurs (sonst Halt und Neubewertung) | 2 bis 3 Wochen |
| M4 | Paarung, MAPQ, Supplementär, ALT, ganzes Genom | simulierte Fehlplatzierung höchstens wie bwa | 2 Wochen |
| M5 | Sortieren, Duplikate, Index, QC verschmolzen; Trimmen; erster eigenständiger Werkzeugbericht und Berichtsrahmen (Abschnitt 14); Einbau in Genome 3.x im Schattenbetrieb neben dem heutigen bwa-Weg | Duplikate und Kennzahlen gleich; App-Lauf auf pj's Probe; Bericht entsteht ohne App | 2 bis 3 Wochen |
| M6 | Eingefrorene Version, volle Testrunde (simuliert, GIAB, MGI, Varianten-Richter), Forschungsbericht, Repo mit DOI | alle Messziele aus 0.4 | 1 bis 2 Wochen |

Zusammen etwa 9 bis 13 Wochen. Größtes Risiko: die vielen kleinen Heuristiken von BWA-MEM (MAPQ, Rescue, Clipping-Strafen). Ziel ist darum nicht bitgleich mit bwa, sondern gemessen gleich gut gegen die Wahrheit.

Vergleichsläufe mit bwa: im Container gibt es bisher nur samtools (über pysam). bwa baue ich dort aus dem Quelltext, wenn das Netz es zulässt; sonst laufen die Vergleiche auf dem Mac, immer erst nach pj's Wort.

### 0.6 In der App

- In Genome 3.x kommt `fq2bam` als Alternative neben den heutigen Weg, zuerst im Schattenbetrieb (beide laufen, Abweichungen sichtbar). Schattenbetrieb heißt hier: je Testprobe einmal beide Wege, nicht bei jedem Lauf (der alte Weg dauert fast 12 Stunden). Nach M6 wird es Standard; bwa, samtools und fastp fallen für diesen Schritt weg.
- In der Rust-App ist es der erste Knoten im Fluss-Graphen (Abschnitt 13.4): Lesungen pro Sekunde, Restzeit, und Kennzahlen wie Duplikatrate und Geschlecht erscheinen schon während des Laufs.

### 0.7 Veröffentlichung

Eigenes Repo (Name offen, z. B. `fq2bam-rs`), MIT/Apache-2.0, Zenodo-DOI, Forschungsbericht auf pjlabs.dev mit dem Simulations- und GIAB-Vergleich. Das BWA-MEM-Verfahren wird in README und Bericht klar zitiert (Li 2013; Vasimuddin et al. 2019), weil es eine Neuumsetzung ist, keine neue Methode.

---

## 1. Ausgangslage: Werkzeuge in Genome 3.0

Quelle: README 3.0, gvali-Prüfung vom 16.11., brca12-Benchmark, HLA-Arbeit. Für FASTQ → BAM ist der Bestand am 18.11. direkt auf dem Mac ausgelesen (`M1-Bestand-Genome3-FASTQ-BAM.md`); die Versionen der übrigen Werkzeuge lese ich jeweils zu Beginn ihrer Welle aus.

| Bereich | Werkzeug heute | Laufzeitbasis | bekannte Probleme |
|---|---|---|---|
| FASTQ-QC, Trimmen | fastp, FastQC | conda, C++/Java | defekte fastp-Version als "Abgebrochen" (behoben) |
| Alignment | bwa 0.7.19 (BWA-MEM2 ist nicht installiert), minimap2 | conda | lange Laufzeit, Wiederholungslogik |
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
3. **Tempo und Speicher:** schneller als das Original auf pj's Mac und lauffähig auf einem Rechner mit 16 GB RAM (Spitzenbedarf unter 12 GB). Tempo wird zusätzlich auf Linux und Windows gemessen.
4. **Robustheit:** 30x, 20x, 15x, 10x (Downsampling); Illumina und MGI/DNBSEQ (pj's eigene Daten sind MGI); BAM und CRAM; GRCh38 mit und ohne ALT; männlich und weiblich.
5. **Ehrlichkeit im Ergebnis:** jede Aussage trägt eine Sicherheitsangabe und die Alternativen, wie beim HLA-Typisierer. Was das Werkzeug nicht kann, steht in der Ausgabe.
6. **Reproduzierbarkeit:** byteidentische Ausgabe bei gleichem Input, Version und Datenstand auf demselben System, inhaltsgleich zwischen den Systemen; Prüfsummen im Manifest.
7. **Schattenbetrieb bestanden:** in der App liefen altes und neues Werkzeug parallel auf echten Proben, die Abweichungen sind erklärt.

Kennzahlen je Bereich: Recall, Precision, F1 (Varianten, nach GIAB-Schichten), Genotyp-Konkordanz, Kopienzahl-Treffer, Allel-Treffer je Feld (HLA/KIR/PGx), Mendel-Fehlerrate in Trios, Laufzeit (Wanduhr und CPU), Spitzen-RAM.

---

## 3. Ehrliche Einordnung je Bereich

| Bereich | besser als Original realistisch? | warum |
|---|---|---|
| Annotation (VEP) | **ja, deutlich** | Regelwerk mit klarer Wahrheit (HGVS, Konsequenz); Gewinn bei Tempo, Größe, Stabilität |
| QC (fastp, FastQC, mosdepth, VerifyBamID2, Somalier) | **ja** | einfache Algorithmen; Gewinn durch einen einzigen Lesedurchgang statt mehrerer |
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
| Small-Variant-Calling | **GATK-Niveau ja, DeepVariant-Niveau schwer** | DeepVariant ist der Stand der Technik; erreichbar nur mit eigenem Modell, Training auf GIAB/Q100, eigener GPU-Inferenz über wgpu |
| Aligner | **gleich genau und schneller ja, genauer nein** | BWA-MEM/minimap2 sind ausgereift; ein neuer Algorithmus wäre Jahre Arbeit für kaum messbaren Gewinn. Sinnvoll: das BWA-MEM-Verfahren in Rust neu umsetzen (Abschnitt 0), Gewinn durch Verschmelzen der Schritte; später pangenom-bewusst |

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
| `gc-log` | Ereignisstrom und Lauf-Protokoll für alle Werkzeuge (Abschnitt 15) |
| `gc-sys` | alles Systemnahe mit einer Umsetzung je System (Abschnitt 12.3) |
| `gc-report` | Berichtsbaustein: gemeinsames Format und Satz, mit dem jedes Werkzeug seinen eigenen Bericht erzeugt (Abschnitt 14) |
| `gc-schema` | JSON-Schemata `genome.<werkzeug>.vN` und Manifest |
| Programme | `fq2bam` (bündelt Trimmen, `galign`, `gbam`, BAM-QC), `gqc`, `galign`, `gbam`, `gcall`, `gsv`, `gcnv`, `gannot`, `ghla` (= pjlabshla), `gkir`, `gpgx`, `gloci` (Lp(a), SMN, …), `gstr`, `ghap`, `gtel`, `gmicro`, `gpgs`, `glift` |

Grundsätze:
- **Ein Programm je Aufgabe, gemeinsame Bausteine.** Jedes Programm ist allein lauffähig (für die Veröffentlichung) und nutzt dieselben Crates (für die App).
- **Aufruf aus der App:** Die Rust-App linkt die Crates direkt (gleiche Typen, kein JSON-Umweg für Zustände), startet schwere Läufe aber als **Arbeitsprozess desselben Programms** (`genome worker <werkzeug>`). Begründung: ein Absturz oder Speicherüberlauf reißt die Oberfläche nicht mit, Abbruch ist sauber, und dieselben Werkzeuge bleiben als eigenständige Programme veröffentlichbar. Fortschritt kommt als typisierte Ereignisse über eine Pipe (behebt die "stillen Phasen", die Hermes bei Aldy fand).
- **macOS, Linux und Windows gleichwertig** (Abschnitt 12.3): SIMD für ARM (NEON) und x86 (AVX2/AVX-512) mit Auswahl zur Laufzeit, alle Kerne, GPU über wgpu (Metal, Vulkan, DirectX 12) für neuronale Teile. Kein Code darf nur auf einem System laufen.
- **Ausgabe:** JSON nach Schema plus Standardformate (VCF, BAM/CRAM, TSV), damit Ergebnisse auch außerhalb von Genome nutzbar sind.
- **Ein Bericht je Werkzeug:** jedes Programm erzeugt aus seinem Ergebnis einen eigenständigen Bericht (PDF und Baustein-Datei), ohne die App. Der Gesamtbericht entsteht nur aus diesen Bausteinen (Abschnitt 14).
- **Determinismus:** feste Zufallssaaten, stabile Sortierung (beim HLA-Port nötig gewesen, siehe typer_v16b).

### 4.2 Verteilung

- Je System ein signiertes Paket (macOS: notarisiertes `Genome.app`/DMG; Windows: signiertes MSI; Linux: AppImage und .deb/.rpm) mit einem Rust-Programm, das Oberfläche und alle Werkzeuge enthält; die veröffentlichten Einzelprogramme sind dieselben Crates mit eigener Kommandozeile. Keine Installation, kein conda, kein Docker. Das löst das Kernproblem aus dem GPT-6-Auftrag (zuverlässige Werkzeugverwaltung) an der Wurzel.
- Große Daten (Referenz, Annotation, Pangenom-Panels, Modelle) lädt Genome wie heute, aber als **versionierte, geprüfte Pakete** aus genome-plus/genome-reflib (Prüfsumme, atomarer Austausch, Wiederaufnahme).
- Zielgröße Annotationsdaten: unter 6 GB statt 27,6 GB VEP-Cache.

### 4.3 Arbeitsdaten

Zwischenergebnisse auf die interne SSD, wenn dort Platz ist (auf pj's Mac sind nur 27 GiB frei, also meist extern); jeder Zwischenstand trägt eine Prüfsumme, damit ein Ausfall der ExtremeSSD erkannt wird und der Lauf wieder aufsetzen kann. Endergebnisse mit Prüfsumme auf den gewählten Speicherort.

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

Für jedes: was ersetzt wird, Ansatz, Validierung, Messziele, Risiko, Veröffentlichung. Basiswerte (Laufzeit, Speicher, Genauigkeit der Originale) messe ich jeweils zu Beginn der Welle auf pj's Mac (für FASTQ → BAM läuft das mit HG001 30x); die Ziele sind deshalb relativ zum Original angegeben.

### 6.1 `gc-io` und `gbam`: BAM/CRAM-Verarbeitung
- **Ersetzt:** samtools (view, sort, index, fixmate, markdup, depth, faidx), htslib, bgzip, tabix, Teile von bcftools (norm, view, index), Picard-Funktionen.
- **Ansatz:** noodles als Grundlage. Kern ist ein **verschmolzener Durchgang**: Aligner-Ausgabe direkt sortieren, Duplikate markieren (samblaster-Logik, mate-aware), indizieren, Abdeckung und QC-Zahlen mitschreiben. Heute ist das eine Pipe aus vier Prozessen plus zwei weitere Lesedurchgänge für Index und flagstat.
- **Validierung:** byteweiser Vergleich der Datensätze mit samtools-Ausgabe (gleiche Flags, gleiche Duplikate) auf 20 Proben; CRAM-Rundreise gegen samtools; Fuzzing mit kaputten Dateien.
- **Ziele:** Datensätze identisch zu samtools markdup; Sortieren, Duplikate, Index und Kennzahlen verlängern den Lauf gegenüber dem reinen Alignment kaum (Zielwert nach der Basismessung); RAM unter 8 GB.
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
- **Ersetzt:** bwa und minimap2 als fremde Programme.
- **Ehrliche Einschätzung:** Ein neuer Alignment-Algorithmus, der BWA-MEM2 an Genauigkeit deutlich schlägt, ist nicht realistisch; das Feld hat 15 Jahre Optimierung hinter sich. Realistisch ist eine **reine Rust-Umsetzung des BWA-MEM-Verfahrens**, die gleich genau ist, auf ARM und x86 schnell läuft (die App nutzt heute das ältere bwa; BWA-MEM2 ist vor allem für x86 optimiert) und mit Trimmen, Sortieren, Duplikaten und QC in einem Prozess verschmilzt.
- **Ansatz in zwei Stufen:**
  1. **Stufe 1 (erster Schritt, Abschnitt 0):** BWA-MEM-Verfahren in Rust (FM-Index, SMEM-Seeds, Chaining, SIMD-Smith-Waterman, Paarung, MAPQ). minimap2 und BWA-MEM2 dienen nur noch als Vergleich im Benchmark.
  2. **Stufe 2 (Welle 6, Forschung):** pangenom-bewusstes Alignment gegen GRCh38 plus HPRC-Haplotypen in schwierigen Regionen (MHC, KIR, CYP2D6, LPA, SMN, Segmentduplikationen), mit Rückprojektion auf GRCh38. Vorbild: vg Giraffe + DeepVariant zeigte weniger Fehler in schweren Regionen. Hier liegt der echte wissenschaftliche Gewinn, nicht in der Grundgeschwindigkeit.
- **Validierung:** wie in Abschnitt 0.3: Konkordanz der Primärplatzierung mit bwa (Ziel über 99,5 %, alle Abweichungen erklärt), vor allem aber Wirkung auf das Variant Calling (F1 nach GIAB-Schichten). Stufe 2: F1 in CMRG-Genen und Q100-schweren Regionen.
- **Ziele:** Stufe 1 wie in Abschnitt 0.4. Stufe 2 messbar weniger Fehler in schweren Regionen bei höchstens 1,5-facher Laufzeit.
- **Risiko:** Stufe 1 mittel (viele Heuristiken, Abschnitt 0.5), Stufe 2 hoch.
- **Veröffentlichung:** Stufe 2 wäre ein eigenes Paper.

### 6.4 `gcall`: kleine Varianten (SNV/Indel)
- **Ersetzt:** bcftools mpileup/call (App), GATK HaplotypeCaller (brca12-wgs).
- **Stand der Technik:** DeepVariant (CNN auf Pileup-Bildern) ist bei Kurzlesungen führend, GATK HC liegt deutlich dahinter bei Indels, bcftools noch weiter. DeepVariant läuft auf dem Mac nur auf der CPU und ist dort langsam.
- **Ansatz in drei Stufen:**
  1. **Haplotyp-Caller:** aktive Regionen, lokale De-Bruijn-Assemblierung, Pair-HMM (SIMD), Genotyp-Likelihoods, gVCF mit Referenz-Konfidenz. Ziel HC-Niveau. Löst sofort Java aus brca12-wgs.
  2. **Lernender Filter:** Gradient-Boosting auf Merkmalen je Kandidat (Allelbalance, Strang, Basenqualität, Kontext, Mappability), trainiert auf GIAB HG001/HG005-7, getestet auf HG002-4 (nie auf Trainingschromosomen).
  3. **Neuronales Modell:** eigenes Netz auf Pileup-Tensoren, trainiert auf GIAB + Q100 + Platinum Pedigree, Inferenz in reinem Rust (Burn, das ein wgpu-Backend hat) auf der GPU jedes Systems, CPU als Rückfall. Hier ist DeepVariant-Niveau das Ziel; dazu Merkmale aus Stufe 2 von `galign` (Pangenom-Kontext).
  - **Abkürzung, falls gewünscht:** Die veröffentlichten DeepVariant-Modelle (BSD-3) über eine eigene Rust-Inferenz betreiben (Modell umgewandelt, GPU über wgpu). Das wäre schnell und genau, aber algorithmisch nicht eigen. Entscheidung für pj (Abschnitt 11).
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
- **Validierung:** alle ClinVar-Varianten mit hinterlegtem HGVS-Ausdruck (mehrere Millionen; genaue Zahl mit dem eingefrorenen ClinVar-Stand): HGVS-Übereinstimmung mit ClinVar und mit VEP; Konsequenz-Konkordanz mit VEP; jeder Unterschied wird einzeln aufgelöst (wer hat recht, belegt mit Transkript). Zusätzlich Mutalyzer-Prüfung einer Stichprobe.
- **Ziele:** HGVS-Übereinstimmung mit ClinVar mindestens so hoch wie VEP; 5 Mio. Varianten (ein WGS) in unter 2 Minuten; Daten unter 6 GB; kein Absturz auf 1000G-Gesamt-VCF.
- **Risiko:** gering bis mittel (viele Sonderfälle, aber klare Wahrheit).
- **Gewinn für pj sofort:** beseitigt Segfaults, 27,6 GB Cache und Perl.
- **Veröffentlichung:** eigenes Repo mit DOI, Forschungsbericht mit dem ClinVar-Abgleich.

### 6.8 `ghla`: HLA (fertig, Einbau offen)
- **Ersetzt:** HLA*LA und T1K.
- **Stand:** pjlabshla v16 in Rust, 48 s, 1,3 GB; 184 ungesehene HPRC-Proben 2-Feld 2937/2941 (T1K 2921/2941), DRB3/4/5-Kopien fast fehlerfrei; Laborwahrheit 1000G läuft.
- **Offen:** Panel-Mapping und Stufe A in Rust (ein Programm, JSON), alle Gene inkl. DRA, E, F, G, H, K, V, GRCh37/no-ALT, Trio-Test, 4-Feld-Kalibrierung, Einbau in die App, Repo "hla-wgs" mit DOI.
- **Rolle im Plan:** erstes fertiges eigenes Werkzeug und Vorbild für das Testdesign. Pilot für Einbau, Signierung, Schema, Bericht und Schattenbetrieb ist jetzt `fq2bam` (Abschnitt 0.6); `ghla` folgt in Welle 2 auf demselben Weg.

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
- **Ansatz:** in der Übergangszeit orchestriert die Swift-App wie heute über `OperationModel`, in der Rust-App die Auftragswarteschlange in `genome-engine`; die Werkzeuge streamen untereinander über Pipes, wo es geht. brca12-wgs bekommt eine Version 2.0 ohne Nextflow: `galign` + `gcall` + `gsv` + `gannot` auf die BRCA-Regionen, eigener BRCA-Bericht wie heute (Abschnitt 14), gleiches Schema `genome.brca12.v1`, GIAB-Benchmark erneut laufen lassen.

---

## 7. Reihenfolge

Sortiert nach Nutzen und Risiko. Aufwand grob in Arbeitswochen konzentrierter Arbeit (Erfahrungswert HLA: Prototyp bis eingefrorene Rust-Version etwa eine Woche intensiver Arbeit, ohne Einbau). Die Zahlen sind Schätzungen, keine Zusagen.

| Welle | Inhalt | Nutzen | Risiko | Aufwand |
|---|---|---|---|---|
| **0+1 FASTQ → BAM** (erster Schritt) | genome-core-Workspace, gc-io, genome-bench, `genome fq2bam` = FASTQ-Prüfung und Trimmen + `galign` Stufe 1 + `gbam` (Sortieren, Duplikate, Index) + BAM-QC, Berichtsrahmen `gc-report` mit dem ersten Werkzeugbericht, Einbau in 3.x im Schattenbetrieb, Repo mit DOI | sehr hoch (Grundlage jeder Analyse, längster Rechenschritt) | mittel | 9 bis 13 |
| **2 Schnelle Gewinne** | `gannot` (VEP weg), restliches `gqc` (Kontamination, Identität), `ghap`, `glift`, `gmicro`, `gtel`, `gpgs`, `ghla` vollständig + Einbau | hoch (Stabilität, 20 GB weniger, weniger conda) | gering | 7 bis 9 |
| **3 Wissenschaftlicher Vorsprung** | `gpgx`, `gkir`, `gstr`, `gloci` (Lp(a) zuerst, dann SMN u. a.) | sehr hoch (Bereiche, in denen die Welt heute schwache Werkzeuge hat) | mittel | 10 bis 14 |
| **4 Kernpipeline** | `gcall` Stufe 1+2, `gcnv`, `gsv`, brca12-wgs 2.0 | hoch (danach kein conda mehr) | mittel bis hoch | 12 bis 16 |
| **5 Spitzenklasse** | `gcall` Stufe 3 (neuronales Modell, GPU über wgpu) | hoch | hoch | 8 bis 12 |
| **6 Forschung** | `galign` Stufe 2 (pangenom-bewusst), SV-Entdeckung mit lokaler Assemblierung | wissenschaftlich hoch | hoch | offen |

**Die App selbst** läuft als eigener Strang parallel (Abschnitt 12.6): Gerüst, Designsystem und Genom-Karte ab Welle 2, Gleichstand mit 3.x bis Ende Welle 4. Aufwand zusätzlich grob 16 bis 24 Wochen; das ist der größte Einzelposten des ganzen Plans.

**Summe:** Wellen 0 bis 5 zusammen 46 bis 64 Arbeitswochen, die App 16 bis 24, zusammen also rund 62 bis 88 Wochen konzentrierter Arbeit (grob anderthalb Jahre); Welle 6 ist offen. Am unsichersten ist die Schätzung für die App (1.134 Swift-Dateien als Vorlage).

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
6. **Berichte:** jedes neue Werkzeug liefert seinen eigenen Bericht und Berichtsbaustein mit (Abschnitt 14); die Swift-App zeigt in der Übergangszeit dessen Inhalt (KIR-Struktur, SMN, Ahnen-PCA) über den bestehenden `ReportDocumentRenderer`; jede Karte nennt Methode und Validierungsstand mit Link zum Forschungsbericht.
7. **Tests:** Vertragstests je Adapter im isolierten Test-Host, Architekturregeln bleiben gültig (Services rufen Programme, Views nie).
8. **Werkzeugverwaltung:** der Bereich "Tools" zeigt für eigene Programme Version und Validierungsstand statt Installationsknöpfen; ToolFinder sucht nur noch für Rückfall-Werkzeuge.

---

## 9. Veröffentlichung

Wie bei brca12-wgs und geplant bei hla-wgs:

- **Eigenes öffentliches Repo je wissenschaftlich eigenständigem Werkzeug:** hla-wgs, kir-wgs (oder in hla-wgs), pgx-wgs, loci-wgs, str-wgs, sv-cnv-wgs, fq2bam (Name offen), annot (Name offen), gcall (Name offen), gqc. Gemeinsame Crates als öffentliches `genome-core`. Entwickelt wird im privaten Repo `genome4`; die öffentlichen Repos entstehen je Werkzeug zum Release als Auszug daraus. Der Quellcode der App bleibt privat.
- **Je Repo:** README, LICENSE (Vorschlag MIT oder Apache-2.0, wie die Bibliotheken, auf denen wir aufbauen), CITATION.cff, CHANGELOG, CI, Release-Tag, Zenodo-DOI pro Release, Benchmark-Skripte zum Nachrechnen jeder Zahl.
- **Forschungsbericht je Werkzeug** auf pjlabs.dev im Forschungsbereich: Methode, Wahrheitssätze, Testdesign mit eingefrorenen Versionen, Ergebnisse mit Intervallen, Grenzen. Für HLA/KIR, PGx und den Variant Caller zusätzlich Preprint auf bioRxiv.
- **Lizenzen der Daten:** ClinVar, gnomAD, GIAB, HPRC, 1000G, PharmVar, IMGT/HLA, IPD-KIR, PGS Catalog sind frei nutzbar (mit Zitierpflicht; IMGT/HLA und IPD-KIR nur unverändert weitergeben, einzelne PGS-Scores tragen eigene Lizenzen). CADD und die vorberechneten SpliceAI-Scores sind nur nicht kommerziell frei; AlphaMissense und REVEL prüfen wir am eingefrorenen Stand. Ob Genome solche Dateien mitliefert, entscheiden wir je Datei.
- **Hinweis überall:** Forschungswerkzeug, kein Medizinprodukt (RUO).

---

## 10. Risiken

| Risiko | Gegenmaßnahme |
|---|---|
| Überanpassung an GIAB/HPRC (Werkzeug lernt die Testproben) | strikte Trennung, Einfrieren vor jeder Runde, Reserveproben, Platinum Pedigree als unabhängige Prüfung |
| Variant Caller erreicht DeepVariant nicht | Stufe 1+2 sind für sich schon ein Gewinn gegenüber bcftools; Entscheidung über DeepVariant-Modell-Abkürzung offen halten |
| Arbeitsumfang (ein Entwickler) | Wellen sind einzeln nutzbar; jede Welle endet mit einem Release, nichts hängt halbfertig |
| ExtremeSSD trennt sich unter Last, interne SSD hat nur 27 GiB frei | Zwischenstände klein halten (höchstens 1,2-mal die BAM), Prüfsumme je Teil, Wiederaufnahme nach Ausfall, atomare Veröffentlichung; auf Dauer hilft nur eine zuverlässige Platte |
| 16 GB RAM auf pj's Mac | kompakter Index, festes RAM-Budget je Lauf, Spitzenbedarf in jedem Benchmark gemessen |
| BWA-MEM-Nachbau erreicht Tempo oder Genauigkeit nicht | Haltepunkt nach M3 (Abschnitt 0.5); Rückfall: der heutige bwa-Weg bleibt, die verschmolzene Kette dahinter ist trotzdem ein Gewinn |
| App-Neubau in Rust ist größer als geschätzt | 3.x bleibt nutzbar, Umstieg Bereich für Bereich mit Gleichstandstest (Abschnitt 12.6) |
| Oberflächen-Bibliothek unreif (Barrierefreiheit, Windows/Linux) | Vergleichsprototyp vor der Festlegung, Oberfläche hinter eigener Schicht, damit ein Wechsel möglich bleibt |
| Linux und Windows nur in CI, nie im echten Gebrauch | ein echter Lauf je System und Meilenstein, nicht nur Tests |
| Neue Werkzeuge verlangen neue QMS-Nachweise | Validierungsbericht je Werkzeug ist zugleich der Nachweis; Traceability je Welle nachziehen |
| CI-Minuten des privaten Repos (macOS zählt zehnfach) | CI nur bei Pull Requests und auf main, lange Benchmarks nie in CI |
| Container-Netz (EBI, NCBI gesperrt) | diese Daten über den Mac holen und in genome-reflib spiegeln |
| Datenlizenzen | Liste je Datenpaket, Nachweis im Manifest |
| MGI-Daten verhalten sich anders als Illumina | MGI-Proben in jede Testrunde (1000G hat keine; pj's Daten und öffentliche DNBSEQ-GIAB-Läufe) |

---

## 11. Entscheidungen für dich

1. **Absicht freigeben:** Genome 4.0 nach diesem Plan, Start mit FASTQ → BAM (Abschnitt 0). **Entschieden am 18.11.** ("starte").
2. **Aligner:** reine Rust-Umsetzung des BWA-MEM-Verfahrens. **Läuft so**, solange du nichts anderes sagst; ein neues Verfahren (Pangenom) kommt später als Stufe 2.
3. **Variant Calling Stufe 3:** eigenes neuronales Modell (eigenständig, länger, Paper-tauglich) oder veröffentlichte DeepVariant-Modelle über eigene Rust-Inferenz (schneller am Ziel, aber nicht algorithmisch eigen). Meine Empfehlung: zuerst Stufe 1+2 bauen, die Frage am Ende von Welle 4 mit Messwerten entscheiden.
4. **Repos:** ein Repo je Werkzeug (wie brca12-wgs und hla-wgs, je eigene DOI) oder ein großes öffentliches Repo mit mehreren DOIs. Meine Empfehlung: je Werkzeug, gemeinsame Crates in `genome-core`.
5. **Oberflächen-Bibliothek:** GPUI oder Iced, entschieden nach einem zweiwöchigen Vergleichsprototyp (Abschnitt 12.2). Meine Vorab-Empfehlung: GPUI.
6. **Rechenort für Benchmarks:** Container für Breitentests, Mac für Tempo und eigene Proben (Empfehlung) oder alles auf dem Mac.

**Stand und nächster Schritt:** Bestand der heutigen Kette ausgelesen, Workspace mit gc-io angelegt, BAM-Rundreise gegen samtools identisch (M1). Offen aus M1: Ergebnis der Basismessung (ganzes HG001, Nacht 18./19.11.). Danach M2: SA-IS und FM-Index.

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
| Oberfläche | SwiftUI | **GPUI** (GPU-gerendert, aus dem Zed-Editor, läuft auf macOS, Linux und Windows, sehr schnell bei großen Datenmengen; Barrierefreiheit und API-Stabilität sind die offenen Fragen) oder **Iced** (wgpu, alle drei Systeme, stabile API, Elm-Architektur). Pflicht im Vergleich: gleiche Qualität auf allen drei Systemen. Entscheidung nach Prototyp: dieselbe Genom-Karte mit pj's BAM in beiden, gemessen an Bildrate, Speicher, Textqualität mit IBM Plex Mono und Barrierefreiheit |
| Genom-Ansicht | IGV extern / keine | eigene GPU-Zeichenfläche (wgpu: Metal, Vulkan, DirectX 12) für Karte, Abdeckung, Lesungen |
| Berichte, PDF | HTML + CSS → WKWebView → PDF | **Typst als Bibliothek** (Rust, Apache-2.0): Satz in Buchqualität, eingebettete IBM Plex Mono, Inhaltsverzeichnis, Quellenkarten, APA-Literatur; dieselbe Engine zeigt den Bericht in der App |
| Lokalisierung | String Catalogs | **Fluent** (Mozilla, Rust): Deutsch und Englisch, grammatikalisch korrekte Mehrzahl |
| Designsystem | PJDesignKit | Crate `pj-design`: Tokens (Terracotta `#D97757`, 14 px Radien, 36 px Aktionsflächen, Plex Mono 400/500, keine Fettschrift) als Rust-Konstanten, dieselben Werte exportiert für Typst und pjlabs.dev |
| Zustand | AppState + @Observable-Modelle | Crate `genome-engine`: Ereignis-gesteuerter Zustand, Läufe als dauerhafte Aufträge (neustartfest, wie im GPT-6-Auftrag gefordert) |
| Speicherung | Dateien + JSON | SQLite (rusqlite) für Projekte, Läufe, Befunde, Manifeste; Ergebnisdateien bleiben Dateien |
| MCP-Anbindung | eigene Swift-Implementierung | offizielles Rust-SDK `rmcp`, standardmäßig aus wie heute |
| Updates | eigene Update-Prüfung | eigener Updater mit Ed25519-Signaturen, Hosts in fester Allowlist |
| Tests | XCTest, SwiftSyntax-Regeln | `cargo test`, Snapshot-Tests der Oberfläche und der Typst-Berichte; die Schichtung erzwingen die Cargo-Abhängigkeiten selbst, `cargo-deny` prüft Lizenzen und unerwünschte Crates |
| Wissensbasen | CSV/JSON im Bundle | unverändert übernommen (KnowledgeBaseCSV, HaplogroupContextCSV), Importskripte nach Rust |

### 12.3 macOS, Linux und Windows

Festlegung pj (18.11.): Genome 4.0 muss so gebaut sein, dass es künftig auch auf Linux und Windows läuft. Regeln ab dem ersten Commit:

- **CI auf allen drei Systemen** (macOS arm64, Linux x86-64 und arm64, Windows x86-64): Bauen, Tests, Lint bei jedem Push. Ein Commit, der auf einem System bricht, gilt als kaputt.
- **Nur reines Rust oder plattformneutrale Crates.** Keine Systemaufrufe direkt im Fachcode; alles Systemnahe (Pfade, Prozesse, Speicherabbild, Dateisperren, freier Speicher, RAM-Größe, Schlaf verhindern während langer Läufe) liegt in einer Crate `gc-sys` mit einer Umsetzung je System.
- **Pfade:** `PathBuf`, nie Zeichenketten mit `/`; lange Windows-Pfade (`\\?\`) und Nicht-ASCII-Namen (iCloud-Pfade, Umlaute) in Tests.
- **Prozesse:** Arbeitsprozesse über `std::process`, Abbruch plattformgerecht (Signale auf Unix, Job Objects auf Windows). Kein `fork`.
- **Dateiformate:** Zeilenende `\n` beim Schreiben, `\r\n` beim Lesen toleriert (FASTQ-Leser macht das schon).
- **SIMD:** gemeinsame Abstraktion mit Umsetzungen für NEON und AVX2/AVX-512, Auswahl zur Laufzeit; Rechenergebnis (Alignments, Kennzahlen) auf allen Systemen gleich (Test in CI).
- **GPU:** nur über wgpu (Metal, Vulkan, DirectX 12); keine Core-ML- oder Metal-only-Wege.
- **Benchmarks** laufen auf pj's Mac und in Linux (Container, CI); ein Windows-Lauf je Meilenstein.
- **Datenordner:** je System der übliche Ort (`~/Library/Application Support/Genome`, `$XDG_DATA_HOME/genome`, `%APPDATA%\Genome`), frei wählbar für externe Platten.

Das kostet anfangs etwa 10 % mehr Aufwand, verhindert aber, dass später ein Umbau nötig wird.

### 12.4 Crate-Schichten

```
pj-design          Designtokens
gc-*               Analysekern (Abschnitt 4), darin gc-report: Berichtsbaustein je Werkzeug
genome-kb          Wissensbasen, Interpretation, ACMG-SF, CPIC-Regeln
genome-engine      Projekte, Läufe, Auftragswarteschlange, Manifest, SQLite
genome-report      setzt den Gesamtbericht aus den Bausteinen der Werkzeuge zusammen (Abschnitt 14)
genome-ui          Oberfläche (darf nur pj-design, genome-engine und genome-report importieren)
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

`⌘K` (Windows und Linux: `Strg+K`) öffnet eine Zeile, die alles versteht, was ein Genom-Mensch tippt: `rs429358`, `APOE`, `chr6:29.9M`, `CYP2D6`, `HLA-B`, `Lp(a)`, `Bericht`, `HG002 vergleichen`. Treffer springen auf der Karte hin und öffnen die Befundkarte. Alles lokal; ein externer KI-Assistent kann über MCP dieselben Befehle nutzen, wenn pj das einschaltet.

### 13.3 Befundkarten mit Beweiskette

Jede Befundkarte hat drei Ebenen, aufklappbar von oben nach unten:

1. **Aussage** in einem Satz, mit Sicherheitsring ("CYP2D6 *1/*4, intermediärer Metabolisierer").
2. **Warum:** die Merkmale, die zur Aussage führen (Allele, Kopienzahl, Alternativen mit Abstand), und die Leitlinie (CPIC, Quelle).
3. **Beweis:** die Lesungen an den entscheidenden Positionen direkt in der Karte, Werkzeug und Version, Datenstand, und der **Validierungsstand des Werkzeugs** ("auf 184 ungesehenen Proben 2937/2941 richtig, Bericht ↗").

Damit wird die Wissenschaftlichkeit, die die Werkzeuge haben, auch in der Oberfläche sichtbar. Ein Schalter je Karte wechselt zwischen verständlicher Sprache und Fachsprache.

### 13.4 Ergebnisse, die während der Analyse entstehen

- Der Lauf ist ein ruhiger Fluss-Graph (FASTQ → Alignment → … → Bericht), jeder Knoten zeigt Fortschritt, Tempo und geschätzte Restzeit aus echten Ereignissen.
- Weil die Werkzeuge streamen, erscheinen Befunde, sobald ihre Daten da sind: Geschlecht, Kontamination und mt-Haplogruppe nach wenigen Minuten; HLA, KIR und PGx aus einer frühen Zielregion-Vorauswahl, noch bevor das ganze Genom aligniert ist. Vorläufige Befunde sind als vorläufig markiert und werden still bestätigt oder korrigiert. Das ist ein Ziel, kein gesicherter Stand: es setzt voraus, dass eine Vorauswahl der Lesungen vor dem Alignment dasselbe Ergebnis liefert wie danach. Das wird je Werkzeug gemessen; wo es nicht stimmt, erscheint der Befund erst nach dem Alignment.
- Unterbrechung, Neustart, Abbruch: der Graph zeigt genau, welche Schritte fertig und gültig sind (neustartfest).

### 13.5 Vergleichen

- **Zwei Proben** nebeneinander auf derselben Karte (z. B. pj und HG002), Unterschiede hervorgehoben.
- **Familie/Trio:** Vererbung je Variante sichtbar, Mendel-Fehler markiert, Phasierung aus den Eltern.
- **Zwei Werkzeuge** (Schattenbetrieb aus Abschnitt 8) oder zwei Datenstände: Abweichungen als Liste mit Sprung zur Stelle.

### 13.6 "Was ist neu?"

Weil die Annotation in Minuten läuft, prüft Genome bei jedem neuen ClinVar-, PharmVar- oder PGS-Catalog-Stand das Genom neu und zeigt eine kurze Liste der geänderten Befunde mit Begründung ("ClinVar stuft Variante X seit 2027-01 als wahrscheinlich gutartig ein"). Kein stilles Überschreiben, jede Änderung bleibt in der Zeitachse.

### 13.7 Bericht als lebendes Dokument

Der Bericht ist eine Ansicht in der App, gesetzt mit derselben Typst-Engine wie das PDF: was man sieht, ist das PDF. Layout-Änderungen (heute der Layout-Editor) wirken sofort. Jede Karte im Bericht verlinkt zurück auf die Stelle in der Karte. Der Gesamtbericht besteht aus den Bausteinen der einzelnen Werkzeuge; jeder Werkzeugbericht lässt sich auch einzeln öffnen und ausgeben (Abschnitt 14).

### 13.8 Ruhe, Privatsphäre, Zugang

- Eine dauerhafte, kleine Anzeige "lokal, offline", die zeigt, ob überhaupt eine Netzverbindung genutzt wird und wofür (nur Datenpakete, nur auf erlaubte Hosts).
- Tastatur zuerst, jede Handlung auch ohne Maus; Bildschirmleser-Unterstützung (VoiceOver, NVDA/Narrator, Orca) über AccessKit. Wie gut die jeweilige Oberflächen-Bibliothek das anbindet, wird im Prototyp geprüft und nicht vorausgesetzt; es ist Pflichtkriterium im Vergleich.
- Hell und dunkel, Kontrast nach WCAG AA, Farbsehschwäche geprüft (Kategorien nie nur über Farbe).

### 13.9 Wie das UX geprüft wird

- Klickbarer Prototyp der Karte und der Befundkarte zu Beginn des App-Strangs (Welle 2) mit pj's echter Probe, bevor die Oberfläche gebaut wird.
- Aufgabentests mit 3 bis 5 Menschen (z. B. "Finde deinen CYP2D6-Befund und zeige, warum er so lautet"), Zeit und Fehler gemessen, vor und nach jeder großen Änderung.
- Tempo-Ziele der Oberfläche: Start unter 1 s, Zoom flüssig in der Bildwiederholrate des Bildschirms (60 bis 120 Bilder/s), jede Suche unter 100 ms.

---

## 14. Berichte: ein Bericht je Werkzeug, daraus der Gesamtbericht

Festlegung pj (18.11.): Genome 4.0 ist so modular, dass jedes Werkzeug (HLA, BRCA, KIR, Lp(a), PGx und alle anderen) einen eigenständigen Bericht erstellt. Der Gesamtbericht entsteht aus diesen Berichten.

### 14.1 Grundsatz

- **Das Werkzeug besitzt seinen Bericht.** Texte, Tabellen, Abbildungen, Methode, Grenzen und Literatur zu HLA liegen beim HLA-Werkzeug, nicht in der App. Wer nur `hla-wgs` aus dem öffentlichen Repo nutzt, bekommt denselben HLA-Bericht wie in Genome.
- **Der Gesamtbericht rechnet nichts.** Er setzt nur zusammen. Jede Zahl im Gesamtbericht stammt aus genau einem Werkzeugbericht und ist dorthin rückverfolgbar.
- **Drei Ebenen je Werkzeug:**
  1. Ergebnis `<werkzeug>.json` (Schema `genome.<werkzeug>.vN`): die Fakten, maschinenlesbar.
  2. Berichtsbaustein `<werkzeug>.report.json` (gemeinsames Schema `genome.report.v1`): die Darstellung, abgeleitet nur aus dem Ergebnis und der Wissensbasis.
  3. Bericht als PDF (und als Ansicht in der App), gesetzt aus dem Baustein.
- Aufruf: `genome <werkzeug> … --report` erzeugt alle drei; `genome <werkzeug> report ergebnis.json` setzt den Bericht neu, ohne neu zu rechnen (z. B. nach neuem Layout oder neuer Sprache).

### 14.2 Der Berichtsbaustein (`genome.report.v1`)

Jeder Baustein hat dieselben Teile, damit alle Berichte gleich aufgebaut sind:

| Teil | Inhalt |
|---|---|
| Kopf | Werkzeug, Version, Commit, Probe, Prüfsumme der Eingabe, Datenstand (z. B. IMGT/HLA-, PharmVar-, ClinVar-Version), Datum, Laufzeit |
| Kernaussagen | die Befundkarten: Aussage in einem Satz, Sicherheit, Alternativen, Bedeutung mit Leitlinie und Quelle |
| Ergebnisse | Tabellen und Abbildungen. Abbildungen als Daten mit Diagrammangabe, nicht als Bild, damit Einzel- und Gesamtbericht sie im selben Stil zeichnen |
| Qualität des Laufs | was für dieses Werkzeug zählt (z. B. Abdeckung im Gen, Anteil eindeutiger Lesungen) und ob die Voraussetzungen erfüllt sind |
| Methode | Verfahren in kurzen Worten, Verweis auf den Forschungsbericht |
| Validierungsstand | Zahlen der eingefrorenen Testrunde ("auf 184 ungesehenen Proben 2937/2941 richtig") |
| Grenzen | was das Werkzeug nicht kann oder in dieser Probe nicht entscheiden konnte |
| Literatur | Quellen in APA, mit DOI |

Texte stehen als Fluent-Schlüssel mit Werten im Baustein, damit Deutsch und Englisch ohne neuen Lauf möglich sind. Jeder Baustein trägt eine Prüfsumme.

### 14.3 Welche Berichte es gibt

| Bericht | Werkzeug |
|---|---|
| Alignment und Datenqualität | `fq2bam` (erster Bericht, Meilenstein M5) |
| Probenqualität: Kontamination, Identität, Geschlecht, Verwandtschaft | `gqc` |
| Kleine Varianten, klinisch relevante Varianten, ACMG-Zusatzbefunde | `gcall` + `gannot` |
| Strukturelle Varianten und Kopienzahl | `gsv`, `gcnv` |
| HLA | `ghla` |
| KIR | `gkir` |
| Pharmakogenetik | `gpgx` |
| Lp(a), SMN1/2 und jedes weitere Spezialgen je ein eigener Bericht | `gloci` |
| Repeat-Expansionen | `gstr` |
| Haplogruppen mt/Y | `ghap` |
| Telomere | `gtel` |
| Viren (EBV u. a.) | `gmicro` |
| Polygene Scores | `gpgs` |
| BRCA1/2 | brca12-wgs 2.0 |

### 14.4 Der Gesamtbericht

`genome report` (Crate `genome-report`) liest die Bausteine und das Manifest einer Probe und setzt daraus:

1. Deckblatt und Probenangaben.
2. **Zusammenfassung:** die wichtigsten Kernaussagen aller Werkzeuge, nach festen Regeln geordnet.
3. **Ein Kapitel je Werkzeug**, wählbar in Kurz- oder Vollform (der heutige Layout-Editor steuert Auswahl und Reihenfolge).
4. **Querschnittskapitel**, wo eine Aussage zwei Werkzeuge braucht: KIR-HLA-Liganden (`gkir` + `ghla`), Arzneimittel-Übersicht (`gpgx` + `ghla`, z. B. HLA-B*57:01), erbliche Tumorprädisposition (BRCA + ACMG-Zusatzbefunde). Diese Kapitel sind eigene kleine Module mit eigenem Baustein; sie lesen nur die Ergebnisse der beteiligten Werkzeuge und nennen sie als Quelle.
5. Methoden-Anhang, Versionstabelle aller Werkzeuge und Datenstände, gemeinsame Literaturliste ohne Doppelte.

Regeln:
- Fehlt ein Werkzeuglauf oder ist er fehlgeschlagen, steht das im Kapitel ("nicht gelaufen", mit Grund). Nichts verschwindet still.
- Passt ein Baustein nicht zur Probe (andere Eingabe-Prüfsumme) oder ist sein Datenstand älter als der der anderen, wird das angezeigt.
- Schlechte Probenqualität (z. B. Kontamination) erscheint als Warnung in jedem betroffenen Einzelbericht, nicht nur im Gesamtbericht: die Werkzeuge nehmen das QC-Ergebnis als Eingabe.
- Läuft ein Werkzeug neu, ändert sich nur sein Baustein; der Gesamtbericht wird in Sekunden neu gesetzt. "Was ist neu?" (Abschnitt 13.6) vergleicht Bausteine.

### 14.5 Technik und Prüfung

- Crate `gc-report`: Typen des Bausteins, Diagramme, Satz mit Typst und `pj-design`. Jedes Werkzeug nutzt sie, deshalb sehen Einzelbericht und Gesamtbericht gleich aus.
- Tests je Werkzeug: Beispiel-Ergebnis → Baustein → PDF als Snapshot, in CI auf allen drei Systemen, ohne App.
- Test des Gesamtberichts: jede Zahl muss in einem Baustein vorkommen; ein absichtlich fehlender und ein veralteter Baustein müssen sichtbar gemeldet werden.
- Jedes öffentliche Werkzeug-Repo liefert einen Beispielbericht (HG002) mit.

### 14.6 Reihenfolge und Aufwand

- Der Rahmen (`gc-report`, Schema, Satz) entsteht mit `fq2bam` in M5; der Alignment-Bericht ist der Pilot. Aufwand etwa eine Woche, in M5 eingerechnet.
- Jedes weitere Werkzeug bringt seinen Bericht in seiner Welle mit. Inhalte und Texte der heutigen Berichte (HLA, BRCA, PGx, Lp(a) usw.) werden übernommen, nicht neu erfunden.
- Die Berichtsarbeit war bisher im App-Strang eingerechnet. Sie wandert von dort zu den Werkzeugen; die Gesamtsumme in Abschnitt 7 bleibt etwa gleich.
- In der Übergangszeit zeigt die Swift-App die Inhalte der Bausteine über ihren heutigen Renderer (Abschnitt 8).

---

## 15. Protokollierung: volle Transparenz, einfache Fehlersuche

Festlegung pj (18.11.): Das Logging muss komplette Transparenz liefern und die Fehlerbehebung einfach machen. Das gilt für jedes Werkzeug und für die App, ab dem ersten Werkzeug.

### 15.1 Grundsätze

- **Eine Quelle:** Jeder Lauf schreibt einen einzigen Ereignisstrom (Crate `gc-log`, Schema `genome.log.v1`). Fortschritt in der App, Protokolldatei und Fehlermeldung sind dieselben Ereignisse, nicht drei getrennte Wege.
- **Zwei Formen derselben Ereignisse:** `run.log.jsonl` (eine JSON-Zeile je Ereignis, maschinenlesbar) und eine gut lesbare Textform. Beide liegen beim Ergebnis.
- **Sofort auf Platte:** jede Zeile wird direkt geschrieben. Nach Absturz, Stromausfall oder Plattenausfall ist bis zur letzten Zeile sichtbar, was lief.
- **Nichts Stilles:** kein Schritt ohne Anfang, Ende und Dauer; bei langen Schritten ein Lebenszeichen alle paar Sekunden mit Fortschritt, Tempo, RAM und freiem Plattenplatz. Daran sieht man sofort, wo ein Lauf hängt.

### 15.2 Was jeder Lauf festhält

| Teil | Inhalt |
|---|---|
| Laufkopf | Lauf-ID, vollständige Befehlszeile, alle Parameter einschließlich der Standardwerte, Programmversion, Commit, System, Prozessor und gewählte SIMD-Stufe, Kerne, RAM, freier Plattenplatz |
| Eingaben | jede Datei mit Pfad, Größe, Änderungszeit, Prüfsumme; Referenz und Datenpakete mit Version und Prüfsumme |
| Entscheidungen | jede automatische Wahl mit Grund ("kompakter Index, weil 16 GB RAM", "Temp-Ordner extern, weil intern nur 27 GiB frei") |
| Stufen | Anfang, Ende, Dauer, Datensätze hinein und heraus, gelesene und geschriebene Bytes, Spitzen-RAM, Temp-Dateien |
| Bilanz | die Zahlen müssen aufgehen (Lesungen hinein = Lesungen heraus + begründet verworfene). Geht die Bilanz nicht auf, ist das ein Fehler, kein Hinweis |
| Warnungen | je Art gezählt, die ersten Beispiele im Wortlaut, danach nur noch die Zahl (keine Flut) |
| Abschluss | Zusammenfassung, Endstatus als eigene Datei (fertig oder fehlgeschlagen mit Fehlercode), der genaue Befehl zum Wiederholen |

### 15.3 Fehler, die sich selbst erklären

Jeder Fehler hat:
- einen **festen Fehlercode** (z. B. `FQ2BAM.FASTQ.PAAR_UNGLEICH`), der in der Dokumentation erklärt ist,
- **was** passiert ist und **wo** (Stufe, Datei, Zeile oder Datensatz, Name der Lesung),
- die **Ursachenkette** bis zum Systemfehler (bei Dateifehlern immer Pfad und Meldung des Systems),
- die **wahrscheinliche Ursache** und **was zu tun ist** ("Datei ist abgeschnitten: neu kopieren und Prüfsumme vergleichen").

Kein nacktes "fehlgeschlagen", kein Exit-Code ohne Text. Programmabstürze werden abgefangen und mit Aufrufstapel protokolliert. Fehlt bei einem Lauf sowohl "fertig" als auch "fehlgeschlagen", erkennt Genome das als Abbruch von außen (Platte weg, Rechner aus) und sagt das so.

### 15.4 Tiefe nach Bedarf

- Stufen: Fehler, Warnung, Info (Standard), Debug, Trace; je Stufe des Programms einzeln einstellbar.
- **Eine Lesung verfolgen:** `--trace-read <Name>` schreibt für diese eine Lesung den ganzen Entscheidungsweg (Seeds, Ketten, Bewertungen, MAPQ, Duplikat-Entscheidung). Damit lässt sich jede Abweichung von bwa oder samtools an der einzelnen Lesung erklären. Dasselbe Prinzip gilt später je Variante, Allel oder Locus.
- Debug- und Trace-Zeilen kosten ohne Einschalten keine Laufzeit.

### 15.5 In der App

- Protokollansicht je Lauf: Zeitachse der Stufen, Filter nach Stufe und Schwere, Sprung vom Fehler zur Ursache.
- Jeder Befund und jeder Berichtsbaustein (Abschnitt 14) trägt die Lauf-ID; von dort führt ein Klick zu den Protokollzeilen, die ihn erzeugt haben.
- **Diagnosepaket:** ein Knopf packt Protokolle, Manifest und Systemangaben in eine Datei für die Fehlersuche, ohne Genomdaten.

### 15.6 Schutz der Daten

Protokolle enthalten standardmäßig keine Sequenzen und keine Genotypen, nur Zahlen, Pfade und Namen von Lesungen. Sequenzen erscheinen nur bei ausdrücklich eingeschaltetem Trace. Protokolle verlassen den Rechner nie von selbst.

### 15.7 Prüfung

- Für jeden Fehlercode gibt es einen Test, der ihn auslöst und Text, Ort und Handlungsempfehlung prüft.
- Jede Protokollzeile wird in CI gegen das Schema geprüft.
- Test "Stecker gezogen": ein mitten im Lauf beendeter Prozess muss ein lesbares Protokoll bis zur letzten Stufe hinterlassen und als Abbruch erkannt werden.
- Grundlage `gc-log` entsteht in M2 und wird von da an von jedem Werkzeug genutzt.

---

## Durchsicht vom 18.11.2026: was berichtigt wurde

- **Aligner widersprüchlich:** Kurzfassung und Abschnitt 3 sagten noch "als Bibliothek einbinden", Abschnitt 0 und 6.3 "in Rust neu umsetzen". Jetzt überall die Rust-Umsetzung; Risiko, Zielwerte und Welle in 6.3 an Abschnitt 0 und 7 angeglichen.
- **Falsche Vergleichsbasis:** der Plan verglich gegen BWA-MEM2, die App nutzt aber bwa 0.7.19. Maßstab ist jetzt bwa; BWA-MEM2 nur noch Tempo-Vergleich.
- **Heutige Kette falsch beschrieben:** FastQC und mosdepth gehören nicht zur Kette, die Schritte laufen als eine Pipe und nicht als getrennte Läufe; "fünf- bis sechsmal gelesen und geschrieben" stimmte nicht.
- **Rechner:** Kernzahlen (10 statt 16), RAM-Ziel (unter 12 statt unter 16 GB), interne SSD zu voll für Zwischenstände.
- **Nur-Apple-Reste entfernt:** Core-ML-Inferenz, "Apple-Chip", "ProMotion", `⌘K` ohne Entsprechung auf Windows und Linux.
- **Sachfehler:** Candle hat kein wgpu-Backend (Burn hat eines); `cargo-deny` prüft keine Schichtung; byteidentische BAM gilt je System, zwischen Systemen nur inhaltsgleich; ClinVar-Zahl; Lizenzangaben zu IMGT/HLA, PGS-Scores und AlphaMissense.
- **Zu sicher formuliert:** frühe HLA-, KIR- und PGx-Befunde vor Ende des Alignments und die Bildschirmleser-Unterstützung sind jetzt als zu prüfende Ziele markiert.
- **Veralteter Stand:** "noch nichts gebaut", "Welle 0", HG002 als Messprobe, bereits getroffene Entscheidungen in Abschnitt 11.
- **Lücken geschlossen:** `gc-sys` und `gc-report` in der Crate-Liste, Gesamtaufwand (rund 62 bis 88 Wochen), sieben fehlende Risiken, `fq2bam` als Einbau-Pilot statt `ghla`.
- **Neu:** Abschnitt 14 (modulare Berichte).
- **Neu:** Abschnitt 15 (Protokollierung).
