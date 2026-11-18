# M1: Bestand der heutigen FASTQ→BAM-Kette in Genome 3.x

Ausgelesen am 18.11.2026 auf pj's Mac, nur gelesen, nichts geändert. Grundlage für die Basismessung und die Messziele von `genome fq2bam` (Plan Abschnitt 0).

## Hardware
- Apple M5, 10 Kerne (4 Performance, 6 Effizienz), 16 GB RAM, macOS 27.2.1.
- Interne SSD 460 GiB, nur 27 GiB frei. ExtremeSSD 1,8 TiB, 786 GiB frei (am 17.11. viermal kurz ausgefallen).

## Kette heute (Einstellungen pj: Aligner bwa, Trimmen an, Temp auf ExtremeSSD, Profil maximum)
1. `bwa index <ref.fa>` falls Index fehlt (Full Analysis Set: 5929 s laut bwa_index.log).
2. `fastp -i R1 -o trimmed_R1 -I R2 -O trimmed_R2 --detect_adapter_for_pe --thread 8 --verbose -h … -j …` als eigener Lauf mit getrimmten Dateien im Temp-Ordner (Streaming abgeschaltet).
3. Eine Shell-Pipe:
   `bwa mem -t 10 -R '<RG>' ref.fa R1 R2 | samtools fixmate -m -u -@ 10 - - | samtools sort -@ 3 -m 448M -T <tmp> - | samtools markdup [-d D --read-coords … --coords-order txy] -f <stats> - <bam>`
   - `-d` nach Sequenzer: NovaSeq/NextSeq 2500, HiSeq/MiSeq 100, MGI 300.
   - Read Group: `ID/SM/LB/PU = <Probe>`, `PL:ILLUMINA` fest, auch bei MGI.
   - Sort ist bei 16 GB RAM auf 3 Threads × 448 MB gedrosselt.
4. `samtools index -@ 10` und `samtools flagstat -@ 10` parallel, dann atomare Veröffentlichung.
- mosdepth und FastQC gehören nicht zur Kette (eigene Aktionen).
- Quellen: AlignmentModel+AlignFASTQOnce.swift, AlignmentPipelineCommandBuilder.swift:108-185, AlignmentCommandBuilder.swift, ExecutionResourceBudget.swift, SequencerDetector.swift:43-51.

## Versionen (arm64, /Volumes/ExtremeSSD/Genome/Tools/conda/bin)
bwa 0.7.19-r1273, minimap2 2.31, samtools 1.23.1, fastp 1.4.0, mosdepth 0.3.14, picard 3.5.0, FastQC 0.13.0. **bwa-mem2 ist nicht installiert.**

## Referenzen (/Volumes/ExtremeSSD/Genome/Ref)
hs38.fa (195 Contigs), **hs38d1.fa** (2580 Contigs, Decoys, ohne ALT; damit liefen HG1 und PJ), hg38.fa (455 Contigs, 261 ALT), hs37d5.fa, hg19.fa. In der Bibliothek liegt **kein Aligner-Index**; ein bwa-Index existiert nur für das Full Analysis Set unter brca12-dev/vollbericht/ref (ohne .alt-Datei).

## Testdaten
- PJ (MGI DNB T7): PJ_1.fq.gz 44,2 GiB, PJ_2.fq.gz 37,6 GiB; PJ.bam 66,1 GiB (hs38d1, ohne Duplikatmarkierung). Teilmengen unter ~/Downloads/pjlabshla.
- GIAB HG001 bis HG007 NovaSeq 30x FASTQ (je rund 25 + 26 GB) unter /Volumes/ExtremeSSD/Entwicklergenome. Keine fertigen BAMs.

## Einzige belegte Laufzeit
HG001 NovaSeq 30x (53,1 GB FASTQ.gz, 713 Mio. Lesungen), Genome 2.0.0, 12.09.2026, bwa mem 0.7.19, fastp, samtools sort/markdup, 10 Threads, hs38d1, Platte damals ExtremePro: **700 min 18 s (11 h 40 min)** gesamt, 10,0 % Duplikate. Zeiten je Schritt sind nicht erhalten.

## Folgen für Genome 4.0
- Die Basis ist **bwa 0.7.19**, nicht BWA-MEM2. Messziele gegen diese Kette (11 h 40 min) und zusätzlich gegen BWA-MEM2 als Stand der Technik.
- 16 GB RAM: BWA-MEM2 braucht für das Humangenom deutlich mehr Arbeitsspeicher als bwa; der kompakte Index (unter 6 GB) ist für diesen Mac Pflicht.
- Interne SSD hat nur 27 GiB frei: Zwischenstände der Sortierung (30x: rund 60 bis 100 GB) passen dort nicht hin. `fq2bam` muss mit externem Temp umgehen (Prüfsummen je Teil, Wiederaufnahme nach Plattenausfall) und so wenig wie möglich zwischenspeichern.
- Auffälligkeiten in 3.x (nicht geändert): `PL:ILLUMINA` auch bei MGI; Hilfetext behauptet bwa-mem2 werde bevorzugt; PJ.bam ohne Duplikatmarkierung; minimap2-Index wird bei jedem Lauf neu gebaut.
