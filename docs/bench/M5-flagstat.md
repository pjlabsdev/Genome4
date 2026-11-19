# M5: Zählung nach Flags

Stand: 19.11.2026. Werkzeug: [`genome flagstat`](../handbuch/werkzeuge/flagstat.md).

## Frage
Gibt `genome flagstat` für jede Zeile dieselbe Zahl wie `samtools flagstat`?

## Aufbau
- Prüfer: `samtools flagstat -O json` 1.24 aus pysam 0.24.1.
- Dateien: die markierte BAM aus [M5 Duplikate](M5-duplikate.md) (260.536 Datensätze), die
  sortierte BAM mit 7 Contigs aus [M5 Sortieren](M5-sortieren.md) (300.000) und eine BAM, die
  das Skript erzeugt: 400.000 Datensätze mit zufälligen Flags (alle zwölf Bits), zufälligem
  Contig und Partner-Contig (auch keines) und MAPQ aus 0, 3, 4, 5, 6, 30, 60.
- `bench/m5/compare_flagstat.py`, Ergebnis `bench/m5/flagstat_vergleich.json`. Verglichen werden
  je Datei 16 Zahlen in zwei Gruppen (Qualitätsprüfung bestanden oder nicht).

## Ergebnis

| Datei | Datensätze (bestanden / nicht) | verglichene Zahlen | Abweichungen |
|---|---|---|---|
| chr20, markiert | 260.536 / 0 | 32 | 0 |
| 7 Contigs | 300.000 / 0 | 32 | 0 |
| zufällige Flags | 200.222 / 199.778 | 32 | 0 |

## Grenzen
Die Zählung ist einfach; die Messung zeigt, dass jede Regel von samtools getroffen ist, auch bei
unsinnigen Kombinationen von Flags. Laufzeit nicht gemessen.

## Wiederholen
```bash
python bench/m5/compare_flagstat.py target/release/genome arbeit bench/m5/flagstat_vergleich.json \
  idx_chr20.bam idx_7contigs.bam
```
