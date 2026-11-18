# M3: CIGAR, NM und MD gegen bwa, mit Indels

Stand: 18.11.2026. Verfahren:
[Globales Alignment und CIGAR](../handbuch/verfahren/globales-alignment.md).

## Aufbau
- Referenz: chr20 aus GRCh38, einzeln ausgerichtete Reads zu 150 Basen, beide Stränge.
- Satz A (wie bisher): 200.000 Reads, 1 % Austauschfehler, keine Indels
  (`simulate_reads.py chr20.fa 200000 150 0.01 1`).
- Satz B (neu): 200.000 Reads, 1 % Austauschfehler, dazu Insertionen und Deletionen von 1 bis 10
  Basen mit Rate 0,002 je Base (`simulate_reads.py chr20.fa 200000 150 0.01 2 0.002`). 48.648
  Reads tragen mindestens ein Indel, 6.262 davon zwei oder mehr.
- Genome 4.0: `genome index ketten ... --tabelle`; bwa: bwa mem aus pybwa 2.3.2 mit Vorgaben,
  Vergleich mit `bench/m3/compare_bwa_align.py`. Verglichen wird das primäre Alignment je Read.

## Ergebnis

| | Satz A (ohne Indels) | Satz B (mit Indels) |
|---|---:|---:|
| Reads | 200.000 | 200.000 |
| Ort und Strang gleich wie bwa | 199.993 | 199.990 |
| davon Score gleich | 199.993 | 199.990 |
| davon CIGAR gleich | 199.993 | 199.990 |
| davon `NM` gleich | 199.993 | 199.990 |
| davon `MD` gleich | 199.993 | 199.990 |
| davon MAPQ gleich | 199.990 | 199.987 |
| Reads, bei denen bwa eine Lücke meldet | 26 | 45.544 |
| davon CIGAR gleich | 26 | 45.544 |
| anderer Ort, bei beiden MAPQ 0 | 7 | 10 |
| anderer Ort, sonst | 0 | 0 |

Wo der Ort übereinstimmt, stimmen Score, CIGAR, `NM` und `MD` in beiden Sätzen ausnahmslos.

## Die Abweichungen in Satz B
Alle 13 Reads liegen in chr20:26,4 bis 29,3 Mb (Satelliten-DNA am Zentromer), wie in Satz A.

- 10 Reads an anderem Ort, bei beiden Programmen MAPQ 0. Bei 8 ist der Score gleich; bei 2 meldet
  Genome 4.0 einen höheren Score als bwa (141 gegen 136, 131 gegen 130). Dort unterscheiden sich
  also die Kandidaten, nicht nur die Wahl zwischen gleich guten. Ursache nicht geklärt.
- 3 Reads mit gleicher Lage und anderer MAPQ (59 gegen 60, 60 gegen 21, 31 gegen 26). Nicht geklärt.

## Gegen die Wahrheit (Satz B)
Der beste Bereich liegt bei 199.283 von 200.000 Reads (99,64 %) am wahren Ort (Toleranz 100
Basen), davon 192.284 eindeutig. Bei 717 Reads liegt er woanders. Die Aufteilung dieser 717 nach
Score am wahren Ort gilt für Satz B nicht: Sie rechnet ohne Lücken und unterschätzt Reads mit
Indels.

## Was der Vergleich nicht zeigt
- Verglichen mit dem in pybwa eingebauten bwa, nicht mit dem Programm bwa 0.7.19.
- Nur das primäre Alignment. Sekundäre und ergänzende Alignments (geteilte Reads), `SA`, `XA`
  und Paare sind nicht verglichen.
- Erzeugte Reads mit gleichmäßig verteilten Fehlern, keine echten Reads, keine Qualitätswerte.
- Indels bis 10 Basen. Längere Lücken und geteilte Reads sind nicht im Satz.
- Tempo: Satz B 92 s bei Genome 4.0, 73 s bei bwa, je ein Kern, nichts optimiert.

Rohdaten: `bench/m3/ergebnisse/vergleich_bwa.json` (Satz A), `vergleich_bwa_indel.json` und
`ketten_indel.json` (Satz B).
