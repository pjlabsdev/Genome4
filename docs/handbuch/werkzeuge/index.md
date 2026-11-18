# `genome index`

Stand: gebaut und geprüft. Messungen: [chr20 gegen bwa](../../bench/M2-fm-index.md),
[ganzes Genom](../../bench/M2-index-ganzes-genom.md).

## Zweck
Baut den FM-Index über das Referenzgenom, mit dem der Aligner exakte Treffer (Seeds) findet.
Ersetzt `bwa index`.

## Aufruf

```bash
genome index build <referenz.fa[.gz]> <index.gcfm> [Optionen]
genome index locate <index.gcfm> <SEQUENZ> [--limit N]
genome index count <index.gcfm> <sequenzen.txt>
genome index smem <index.gcfm> --sequenz <SEQUENZ> | --fastq <reads.fq[.gz]> [Optionen]
genome index ketten <index.gcfm> --fastq <reads.fq[.gz]> [--wahrheit-im-namen] [--max-orte N] [--fehlende DATEI] [--tabelle DATEI]
```

Optionen von `build`:

| Option | Standard | Bedeutung |
|---|---|---|
| `--variante kompakt\|schnell` | `kompakt` | Größe gegen Tempo, siehe [FM-Index](../verfahren/fm-index.md) |
| `--contig NAME` | alle | nur diese Contigs aufnehmen (mehrfach möglich) |
| `--verfahren auto\|direkt\|blockweise` | `auto` | `direkt` hält das volle Suffix-Array im RAM (bis 300 Mio. Basen), `blockweise` arbeitet mit festem RAM für das ganze Genom |
| `--block-mio N` | 96 | blockweise: Millionen Suffixe je Block (16 MB RAM je Million) |
| `--kerne N` | 0 = alle | Kerne für den blockweisen Bau |
| `--pruefen N` | 100000 | Zahl der Prüfabfragen (0 = aus), siehe Validierung |
| `--vollpruefung` | aus | die fertige BWT vollständig gegen die Referenz prüfen |
| `--log-stufe` | `info` | Ausführlichkeit des Protokolls |

Optionen von `smem` (Seeds eines Reads, siehe [Seed-Suche](../verfahren/seed-suche.md)):

| Option | Standard | Bedeutung |
|---|---|---|
| `--sequenz SEQ` | | ein Read; Ausgabe sind seine Seeds mit Orten |
| `--fastq DATEI` | | alle Reads der Datei; Ausgabe ist nur die Zusammenfassung |
| `--min-laenge N` | 19 | kürzester gemeldeter Treffer |
| `--orte N` | 20 | mit `--sequenz`: höchstens so viele Orte je Treffer |
| `--pruefen` | aus | mit `--fastq`: jeden Treffer mit der einfachen Suche nachprüfen |

## Eingaben
FASTA, unkomprimiert oder gzip. Groß- und Kleinbuchstaben gelten gleich. Alles außer A, C, G, T
(N und Mehrdeutigkeits-Codes) wird als N-Strecke vermerkt.

## Ausgaben
- `<index.gcfm>`: [Indexdatei](../formate/gcfm.md), mit Prüfsumme.
- `<index.gcfm>.log.jsonl`: [Lauf-Protokoll](../protokollierung.md).
- `build`: JSON, Schema `genome.index-build.v1`: `variante`, `verfahren`, `referenz_basen`,
  `contigs`, `n_strecken`, `zeilen`, `basen_je_zaehlblock`, `sa_abstand`, `datei_bytes`,
  `bits_je_base`, `selbstpruefung`, `vollpruefung`, `seconds`, `lauf` (Lauf-ID).
- `locate`: JSON, Schema `genome.index-locate.v1`: `fundstellen_roh` und `treffer` mit `contig`,
  `pos` (1-basiert, linkes Ende auf dem Vorwärtsstrang), `strang`.
- `smem`: JSON, Schema `genome.index-smem.v1`. Mit `--sequenz`: `seeds` mit `read_von`, `read_bis`
  (1-basiert, einschließlich), `laenge`, `fundstellen_roh`, `orte`. Mit `--fastq`: `reads`, `basen`,
  `seeds`, `seed_basen`, `fundstellen_roh`, `reads_ohne_seed`, `sekunden`, `reads_je_sekunde`,
  `kerne`, `nachgeprueft`, `abweichungen`.
- `ketten`: JSON, Schema `genome.index-ketten.v1`, Zusammenfassung über alle Reads: `reads`,
  `seeds`, `ketten`, `ketten_nach_filter`, `bereiche`, `reads_ohne_kette`, `sekunden`, `reads_je_sekunde`, `max_orte`; mit
  `--wahrheit-im-namen` zusätzlich `wahrheit` (siehe [Seeds und Ketten](../verfahren/seeds-und-ketten.md)).
  `--tabelle` schreibt je Read den primären Bereich: Name, Contig, Position, Strang, MAPQ, Score,
  gemeldeter Score, Read-Anfang, Read-Ende, `sub`, `sub_n`, Anteil in sehr häufigen Seeds,
  CIGAR, `NM`, `MD`. Die Position ist die des fertigen Alignments.
  Ein Prüfbefehl für die Entwicklung; er entfällt, wenn `fq2bam` steht.
- `count`: je Zeile `Sequenz<TAB>rohe Trefferzahl` (beide Stränge).

Rohe Fundstellen und gefilterte Treffer unterscheiden sich: siehe [Glossar](../glossar.md).

## Verfahren
1. FASTA lesen, Referenz als 2-Bit-Folge mit Contig-Tabelle und N-Strecken.
2. Text = Vorwärtsstrang + reverses Komplement + Endzeichen.
3. Sortierte Suffixe: [SA-IS](../verfahren/suffix-array-sais.md) (direkt) oder
   [blockweise](../verfahren/blockweiser-bau.md).
4. BWT, Zähltabelle und Stichproben: [FM-Index](../verfahren/fm-index.md).
5. Prüfen, dann atomar schreiben.

## Validierung

| Prüfung | wann | was sie belegt |
|---|---|---|
| Suffix-Array vollständig geprüft | direkter Bau mit `--pruefen > 0` | jede Position kommt genau einmal vor, Reihenfolge stimmt |
| Abfragen gegen binäre Suche im vollen Suffix-Array | direkter Bau | Suche und Ortung des FM-Index stimmen mit einem unabhängigen Weg überein |
| Stichproben gegen die Referenz | blockweiser Bau | Stücke der Referenz werden an ihrer Herkunft gefunden; an jeder georteten Stelle steht das Muster |
| Vollprüfung (`--vollpruefung`) | beide | die BWT ergibt rückwärts Base für Base die Referenz und bestätigt jede Suffix-Array-Stichprobe |
| blockweise gleich direkt | Tests und chr20 | beide Bauwege liefern dieselbe Datei, Byte für Byte |
| Vergleich mit bwa | Messung | gleiche Trefferzahlen und Positionen, Abweichungen erklärt |

Zahlen stehen in den beiden Messdokumenten oben.

## Grenzen
- Der Index wird beim Laden vollständig eingelesen; Speicherabbild folgt mit `gc-sys`.
- Die Suche ist noch nicht auf Tempo ausgelegt (Seeds für den Aligner kommen in M3).
- ALT-Contigs werden wie jedes andere Contig behandelt; die ALT-bewusste Behandlung gehört zum Aligner.
- Die N-Füllung unterscheidet sich von der von bwa. Rohe Trefferzahlen können deshalb bei
  Wiederholungen neben Lücken um einzelne Fundstellen abweichen; gefilterte Treffer nicht.

## Fehler und was zu tun ist
`INDEX.CONTIG_FEHLT`, `INDEX.BAU_FEHLGESCHLAGEN`, `INDEX.PRUEFUNG_FEHLGESCHLAGEN`: siehe
[Ereignisse und Fehlercodes](../ereignisse-und-fehlercodes.md). Eine beschädigte oder abgeschnittene
Indexdatei wird beim Laden an der Prüfsumme erkannt.

## Literatur
Burrows & Wheeler 1994; Ferragina & Manzini 2000; Nong, Zhang & Chan 2009; Kärkkäinen 2007;
Li & Durbin 2009. Siehe [Literatur](../literatur.md).
