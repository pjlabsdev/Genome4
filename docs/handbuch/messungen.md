# Verzeichnis der Messungen

Jede Messung liegt als Dokument unter `docs/bench/`, ihre Skripte und Rohdaten unter `bench/`.

| Meilenstein | Frage | Ergebnis in einem Satz | Dokument | Skripte und Rohdaten |
|---|---|---|---|---|
| M1 | Wie macht Genome 3.x heute aus FASTQ eine BAM? | fastp, dann `bwa mem` 0.7.19 in einer Pipe mit samtools 1.23.1; einzige belegte Laufzeit 11 h 40 min für HG001 30x | [Bestand](../bench/M1-Bestand-Genome3-FASTQ-BAM.md) | `bench/m1/run_baseline.sh` |
| M1 | Schreibt `gc-io` eine BAM identisch zu samtools? | ja, 844.606 Datensätze und Kopf identisch (GIAB HG002, chr20:10-12 Mb) | [BAM-Rundreise](../bench/M1-bam-rundreise.md) | `bench/m1/compare_bam.py` |
| M1 | Wie lange braucht jeder Schritt der heutigen Kette? | Messung angesetzt (ganzes HG001, Nacht 18./19.11.2026) | folgt | `bench/m1/run_baseline.sh` |
| M2 | Sind Suffix-Array und FM-Index richtig? | chr20: Suffix-Array vollständig geprüft, 1 Mio. Abfragen ohne Abweichung, gegen bwa 65 von 1 Mio. rohen Trefferzahlen verschieden und alle erklärt | [FM-Index](../bench/M2-fm-index.md) | `bench/m2/` |
| M2 | Lässt sich der Index für das ganze Genom mit festem RAM bauen und ist er richtig? | GRCh38: Bau 18,5 min auf 4 Kernen, Spitze 8,68 GB RAM, Datei 4,84 GB; Stichprobe (1 Mio. Abfragen) und Vollprüfung (6,2 Mrd. Zeichen) ohne Abweichung | [Index ganzes Genom](../bench/M2-index-ganzes-genom.md) | `bench/m2/ergebnisse/` |
| M3 | Findet die Seed-Suche genau die längsten exakten Treffer? | chr20, 200.000 erzeugte Reads: 471.822 Treffer einzeln nachgeprüft, 0 Abweichungen; 14.700 Reads je Sekunde auf einem Kern (nicht optimiert) | [Seed-Suche](../bench/M3-seed-suche.md) | `bench/m3/` |
| M3 | Führen Seeds und Ketten zum wahren Ort? | chr20, 200.000 erzeugte Reads: 99,95 % mit einer Kette am wahren Ort, bei 96,1 % ist sie eindeutig die schwerste; die 74 Reads ohne Seed am wahren Ort sind einzeln erklärt | [Ketten](../bench/M3-ketten.md) | `bench/m3/` |
| M3 | Liegt der beste Alignment-Bereich am wahren Ort? | chr20, 200.000 erzeugte Reads: 99,67 % (96,24 % eindeutig, 3,43 % gleich gut wie ein anderer Ort); 32 Reads (0,016 %) verpasst | [Bereiche](../bench/M3-bereiche.md) | `bench/m3/` |
