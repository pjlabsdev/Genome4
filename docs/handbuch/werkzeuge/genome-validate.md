# `genome-validate`

Stand: gebaut

## Zweck
Ein eigenes Prüfprogramm neben `genome`. Es erzeugt Daten mit bekannter Antwort und misst die
Werkzeuge daran. Zahlen zur Genauigkeit in Messberichten und PRs kommen aus diesem Programm,
nicht aus Skripten. Ein Prüfer, der vom Geprüften abhängt, belegt nichts; deshalb nutzt es nur
`gc-io` zum Lesen und Schreiben, keinen Baustein des Ausrichtens.

## Aufruf

```
genome-validate reads-erzeugen --referenz <fasta> --ausgabe <prefix> [Optionen]
genome-validate ausrichtung-pruefen --contigs <prefix>.contigs.tsv --bam name=pfad [--bam ...]
                                    [--bereich name=bed ...] [--toleranz 20] [--json datei]
```

`reads-erzeugen`:

| Option | Vorgabe | Bedeutung |
|---|---|---|
| `--paare` | 1.000.000 | Zahl der Read-Paare |
| `--laenge` | 150 | Read-Länge |
| `--insert`, `--streuung` | 400, 60 | Länge der Stücke, normalverteilt |
| `--saat` | 1 | Saat; gleiche Saat und Referenz geben auf jedem System dieselben Dateien |
| `--profil` | `illumina` | Fehlerprofil `illumina` oder `mgi` |
| `--snv-rate` | 0,001 | SNVs je Base der Probe |
| `--indel-rate` | 0,000125 | kleine Indels je Base der Probe |
| `--max-indel` | 10 | längstes Indel der Probe |

`ausrichtung-pruefen`: die erste BAM wird mit jeder weiteren Read für Read verglichen.
`--bereich` zählt zusätzlich nur Reads, deren wahres linkes Ende in der BED-Datei liegt (zum
Beispiel die schweren Bereiche aus den GIAB-Schichten).

Beispiel:

```
genome-validate reads-erzeugen --referenz chr20.fa --ausgabe w/c20 --paare 200000 --saat 3
genome fq2bam chr20.kompakt.gcfm --fastq w/c20_1.fq.gz --fastq2 w/c20_2.fq.gz --aus w/c20.g4.bam
genome-validate ausrichtung-pruefen --contigs w/c20.contigs.tsv --bam genome4=w/c20.g4.bam \
    --bam bwa=w/c20.bwa.bam --json w/wahrheit.json
```

## Eingaben
- Referenz als FASTA, auch gzip. Für die Prüfung muss die BAM gegen dieselbe Referenz
  ausgerichtet sein (Contigs werden über den Namen zugeordnet).
- BAMs in beliebiger Sortierung. Read 1 und Read 2 werden über die Flags 64 und 128 erkannt.

## Ausgaben

`reads-erzeugen` schreibt:

| Datei | Inhalt |
|---|---|
| `<prefix>_1.fq.gz`, `<prefix>_2.fq.gz` | Read-Paare, Name `w.<nr>.<contig1>.<pos1>.<contig2>.<pos2>` |
| `<prefix>.contigs.tsv` | Nummer, Name und Länge jedes Contigs |
| `<prefix>.varianten.tsv` | alle Varianten der Probe wie VCF geschrieben (Indels mit Ankerbase), Genotyp phasiert |
| `<prefix>.json` | Einstellungen und Zahlen, Schema `genome.validate.reads.v1` |

`ausrichtung-pruefen` schreibt eine Tabelle auf die Konsole und mit `--json` eine Datei mit Schema
`genome.validate.mapping.v1`:

| Feld | Bedeutung |
|---|---|
| `toleranz` | größter erlaubter Abstand zum wahren linken Ende in Basen |
| `bams[].gesamt` | `reads`, `unmapped`, `wrong`, `mapq20` (Reads mit Ort und MAPQ ab 20), `wrong_mapq20` |
| `bams[].stufen[]` | je MAPQ-Stufe `reads`, `wrong` und `expected_wrong`, die Zahl falscher Reads, die die MAPQ verspricht (Summe von 10^(-MAPQ/10)) |
| `bams[].je_mapq[]` | dieselben Zahlen je einzelner MAPQ |
| `bams[].bereiche[]` | Name, Basen und die Zahlen von `gesamt` für jede BED-Datei |
| `bams[].uebersprungen_sekundaer_supplementaer` | nicht gezählte Datensätze |
| `bams[].namen_ohne_wahrheit` | Datensätze ohne lesbaren Namen (sollte 0 sein) |
| `vergleiche[]` | je Paar aus erster und weiterer BAM: beide richtig, nur eine richtig, beide nicht; sicher (MAPQ ab 20) falsch in der einen und richtig in der anderen |

## Verfahren
Siehe [Wahrheit erzeugen und Ausrichtung prüfen](../verfahren/wahrheit-simulieren.md).

## Validierung
Tests im Quelltext (`cargo test -p genome-validate`):
- Stück eines Haplotyps ohne Varianten gleich der Referenz; SNV, Lücke und Einfügung wirken nur
  auf ihren Haplotyp, die Referenzorte jeder Base stimmen; ein Stück, das in einer Lücke beginnt,
  überspringt sie; Stücke mit N oder über das Contig-Ende werden abgelehnt.
- Die gemeldete Qualität hält, was sie verspricht: gemessene Fehler je Stufe gleich 10^(-Q/10).
- Gleiche Saat gibt dieselbe Zufallsfolge (fester Wert im Test).
- Weich und hart abgeschnittene Basen werden zum linken Ende zurückgerechnet; BED-Bereiche werden
  verschmolzen; der Vergleich zählt jeden Read einmal.

Messung: [M4 Wahrheit](../../bench/M4-wahrheit.md).

## Grenzen
- Das Fehlerprofil ist ein Modell (vier Qualitätsstufen, zum Read-Ende hin öfter niedrig, seltene
  Einfügungen und Lücken von einer Base), nicht aus echten Läufen gelernt. Echte Fehler hängen
  zusätzlich von der Sequenz ab (Homopolymere, GGC-Motive) und kommen gehäuft vor.
- Die Probe hat nur SNVs und kleine Indels. Strukturvarianten, Kopienzahl und Sequenz, die in der
  Referenz fehlt, gibt es nicht; dadurch ist die Ausrichtung hier leichter als bei echten Menschen.
- Richtig heißt: Contig gleich und linkes Ende höchstens `--toleranz` Basen daneben. Ein Read mit
  richtigem Ort und falschem CIGAR zählt als richtig.
- Erzeugt wird aus der Referenz, gegen die auch ausgerichtet wird. Reads aus Contigs wie chrEBV
  oder den Decoys entstehen nach ihrem Anteil an der Länge.

## Fehler und was zu tun ist
Das Programm bricht mit einer Meldung ab, wenn die Referenz leer ist, die Stücke kürzer als ein
Read wären, eine BAM keinen einzigen Read mit Wahrheit im Namen hat oder die Contig-Liste nicht
fortlaufend nummeriert ist.

## Literatur
Vorgehen wie bei wgsim und Mason (simulierte Reads mit bekanntem Ort), siehe
[Literatur](../literatur.md).
