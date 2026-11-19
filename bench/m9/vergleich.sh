#!/usr/bin/env bash
# Vergleich genome trim gegen fastp 1.4.0 an den ersten N Paaren einer Probe.
# Aufruf: vergleich.sh FASTP GENOME R1 R2 N ZIELORDNER
# fastp läuft genau wie in Genome 3 (--detect_adapter_for_pe --thread 8), dazu nur -j/-h.
set -euo pipefail
fastp=$1; genome=$2; r1=$3; r2=$4; n=$5; ziel=$6
mkdir -p "$ziel"
"$fastp" --version 2>&1 | tee "$ziel/fastp-version.txt"
lines=$((n * 4))
# Ausschnitt; gzip -dc liest auch unkomprimierte Dateien nicht, daher je nach Endung.
cut_head() { case "$1" in *.gz) gzip -dc "$1" | head -n "$lines" ;; *) head -n "$lines" "$1" ;; esac; }
cut_head "$r1" > "$ziel/in_1.fq" || true
cut_head "$r2" > "$ziel/in_2.fq" || true
/usr/bin/time -l "$fastp" -i "$ziel/in_1.fq" -I "$ziel/in_2.fq" -o "$ziel/fastp_1.fq" -O "$ziel/fastp_2.fq" \
  --detect_adapter_for_pe --thread 8 -j "$ziel/fastp.json" -h "$ziel/fastp.html" 2> "$ziel/fastp.err"
/usr/bin/time -l "$genome" trim --fastq "$ziel/in_1.fq" --fastq2 "$ziel/in_2.fq" \
  --aus "$ziel/genome_1.fq" --aus2 "$ziel/genome_2.fq" > "$ziel/genome.json" 2> "$ziel/genome.err"
python3 "$(dirname "$0")/trim-vergleich.py" "$ziel/fastp_1.fq" "$ziel/fastp_2.fq" \
  "$ziel/genome_1.fq" "$ziel/genome_2.fq" "$ziel/fastp.json" "$ziel/genome_1.fq.log.jsonl" > "$ziel/vergleich.json"
rm -f "$ziel/in_1.fq" "$ziel/in_2.fq" "$ziel/fastp_1.fq" "$ziel/fastp_2.fq" "$ziel/genome_1.fq" "$ziel/genome_2.fq"
echo "fertig: $ziel/vergleich.json"
