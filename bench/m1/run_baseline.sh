#!/bin/bash
# run_baseline.sh: Basismessung der FASTQ->BAM-Kette von Genome 3.x ausserhalb der App.
#
# Aufruf:
#   run_baseline.sh REF_FASTA R1 R2 PROBE WORKDIR [THREADS=10] [OPT_DIST=2500]
#
# Bildet die App-Kette nach (Genome 3.0, Stand 2026-11-18):
#   AlignmentCommandBuilder.swift:56-95           fastp-Trim (Argumentreihenfolge, --thread min(8, THREADS))
#   AlignmentPipelineCommandBuilder.swift:108-185 bwa mem | fixmate | sort | markdup
#   AlignmentPipelineCommandBuilder+Markdup.swift:8-29  markdup -d / --read-coords / --coords-order
#   AlignmentModel+ReferenceOutput.swift:120      Read Group
#   SequencerDetector.swift:96-97,106-108         Regex illuminaModernCoords, NovaSeq 6000 -> -d 2500, txy
#   HardwareDetector.swift:95-129                 sort -@ 3 -m 448M (9 bis 24 GB verfuegbarer RAM), kein -u bei externem Temp
#   AlignmentCommandBuilder.swift:7-33            samtools index / flagstat mit -@ THREADS
#
# Schreibt nur unter WORKDIR. Werkzeuge werden nur gelesen.

set -Eeuo pipefail

usage() {
    echo "Aufruf: $0 REF_FASTA R1 R2 PROBE WORKDIR [THREADS=10] [OPT_DIST=2500]" >&2
}
if [ $# -lt 5 ] || [ $# -gt 7 ]; then usage; exit 2; fi

REF_FASTA=$1
R1=$2
R2=$3
PROBE=$4
WORKDIR=$5
THREADS=${6:-10}
OPT_DIST=${7:-2500}

VOLUME=/Volumes/ExtremeSSD
TOOLS=$VOLUME/Genome/Tools/conda/bin
BWA=$TOOLS/bwa
SAMTOOLS=$TOOLS/samtools
FASTP=$TOOLS/fastp
TIME=/usr/bin/time

# Wortlich aus SequencerDetector.swift (illuminaModernCoords, NovaSeq 6000).
READ_COORDS=':([0-9]+:[a-zA-Z0-9]{9}:[0-9]:[0-9]+):([0-9]+):([0-9]+)'
COORDS_ORDER=txy
# ComputeProfile auf 16-GB-Mac (9 bis 24 GB verfuegbar): sort -@ min(3, THREADS), -m 1344/3 = 448M.
SORT_THREADS=3
SORT_MEM=448M
if [ "$THREADS" -lt "$SORT_THREADS" ]; then SORT_THREADS=$THREADS; fi
# AlignmentCommandBuilder.fastpThreadCount: max(1, min(8, threads)).
FASTP_THREADS=$THREADS
if [ "$FASTP_THREADS" -gt 8 ]; then FASTP_THREADS=8; fi
# Read Group wie AlignmentModel+ReferenceOutput.swift:120 (\t bleibt Backslash-t, bwa wandelt um).
RG="@RG\tID:${PROBE}\tSM:${PROBE}\tPL:ILLUMINA\tLB:${PROBE}\tPU:${PROBE}"

REF=$WORKDIR/ref/ref.fa
TMP=$WORKDIR/tmp
OUT=$WORKDIR/out
LOGS=$WORKDIR/logs
T1=$TMP/trimmed_R1.fq.gz
T2=$TMP/trimmed_R2.fq.gz
SORT_PREFIX=$TMP/genome_sort_$PROBE
BAM=$OUT/$PROBE.markdup.bam
TIMING=$LOGS/timing.tsv
DISK=$LOGS/disk.tsv
INFO=$LOGS/run_info.txt

CURRENT_STEP=init
STEP_ERR_FILES=""

now()  { perl -MTime::HiRes=time -e 'printf "%.3f", time'; }
iso()  { date +%Y-%m-%dT%H:%M:%S%z; }
secs() { perl -e 'printf "%.3f", $ARGV[1] - $ARGV[0]' "$1" "$2"; }

check_volume() {
    if [ ! -d "$VOLUME/Genome" ]; then
        echo "$(iso) ExtremeSSD fehlt, Abbruch" >&2
        exit 3
    fi
}

on_err() {
    local rc=$?
    trap - ERR
    local pids
    pids=$(jobs -p || true)
    if [ -n "$pids" ]; then kill $pids 2>/dev/null || true; fi
    if [ ! -d "$VOLUME/Genome" ]; then exit "$rc"; fi
    {
        echo "step=$CURRENT_STEP"
        echo "exit_code=$rc"
        echo "time=$(iso)"
        for f in $STEP_ERR_FILES; do
            echo "--- letzte 50 Zeilen: $f"
            tail -n 50 "$f" 2>/dev/null || true
        done
    } > "$WORKDIR/FAILED"
    exit "$rc"
}
trap on_err ERR

cmdline() {
    local s="" a
    for a in "$@"; do s="$s $(printf '%q' "$a")"; done
    echo "${s# }"
}

disk_snapshot() {
    local label=$1 line used avail tmpkb
    line=$(df -k "$VOLUME" | tail -n 1)
    used=$(echo "$line" | awk '{print $3}')
    avail=$(echo "$line" | awk '{print $4}')
    tmpkb=$(du -sk "$TMP" 2>/dev/null | cut -f1)
    printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$(iso)" "$used" "$avail" "$tmpkb" >> "$DISK"
}

STEP_T0=""
STEP_ISO0=""
step_begin() {
    check_volume
    CURRENT_STEP=$1
    STEP_ERR_FILES=$2
    STEP_ISO0=$(iso)
    STEP_T0=$(now)
}
step_end() {
    local note=${1:-}
    local t1 iso1
    t1=$(now)
    iso1=$(iso)
    printf '%s\t%s\t%s\t%s\t%s\n' "$CURRENT_STEP" "$STEP_ISO0" "$iso1" "$(secs "$STEP_T0" "$t1")" "$note" >> "$TIMING"
    disk_snapshot "after_$CURRENT_STEP"
}

# ---------- Vorpruefung ----------
check_volume
for d in "$WORKDIR/ref" "$TMP" "$OUT" "$LOGS"; do
    if [ ! -d "$d" ]; then echo "Ordner fehlt: $d" >&2; exit 2; fi
done
for f in "$REF_FASTA" "$R1" "$R2"; do
    if [ ! -r "$f" ]; then echo "Eingabe nicht lesbar: $f" >&2; exit 2; fi
done
for t in "$BWA" "$SAMTOOLS" "$FASTP"; do
    if [ ! -x "$t" ]; then echo "Werkzeug fehlt: $t" >&2; exit 2; fi
done
for m in DONE FAILED; do
    if [ -e "$WORKDIR/$m" ]; then echo "$WORKDIR/$m existiert schon, Abbruch ohne Ueberschreiben" >&2; exit 2; fi
done

RUN_T0=$(now)
printf 'step\tstart\tend\tseconds\tnote\n' > "$TIMING"
printf 'label\ttime\tdf_used_kb\tdf_avail_kb\ttmp_kb\n' > "$DISK"
disk_snapshot before_start

CMD_FASTP=("$FASTP" -i "$R1" -o "$T1" -h "$OUT/${PROBE}_fastp.html" -j "$OUT/${PROBE}_fastp.json" --thread "$FASTP_THREADS" --verbose -I "$R2" -O "$T2" --detect_adapter_for_pe)
CMD_BWA=("$BWA" mem -t "$THREADS" -R "$RG" "$REF" "$T1" "$T2")
CMD_FIXMATE=("$SAMTOOLS" fixmate -m -u -@ "$THREADS" - -)
CMD_SORT=("$SAMTOOLS" sort -@ "$SORT_THREADS" -m "$SORT_MEM" -T "$SORT_PREFIX" -)
CMD_MARKDUP=("$SAMTOOLS" markdup -d "$OPT_DIST" --read-coords "$READ_COORDS" --coords-order "$COORDS_ORDER" -f "$OUT/${PROBE}_markdup.txt" - "$BAM")
CMD_INDEX=("$SAMTOOLS" index -@ "$THREADS" "$BAM")
CMD_FLAGSTAT=("$SAMTOOLS" flagstat -@ "$THREADS" "$BAM")

{
    echo "== Lauf"
    echo "probe=$PROBE"
    echo "workdir=$WORKDIR"
    echo "threads=$THREADS fastp_threads=$FASTP_THREADS sort_threads=$SORT_THREADS sort_mem=$SORT_MEM opt_dist=$OPT_DIST"
    echo "start=$(iso)"
    echo "script_sha256=$(shasum -a 256 "$0" | awk '{print $1}')"
    echo
    echo "== Hardware"
    echo "cpu=$(sysctl -n machdep.cpu.brand_string)"
    echo "hw.ncpu=$(sysctl -n hw.ncpu)"
    echo "hw.perflevel0.physicalcpu=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo n/a)"
    echo "hw.perflevel1.physicalcpu=$(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null || echo n/a)"
    echo "hw.memsize=$(sysctl -n hw.memsize)"
    echo "os=$(sw_vers -productVersion 2>/dev/null || true) $(uname -r)"
    echo "vm_stat (Start):"
    vm_stat | head -n 6
    echo
    echo "== Werkzeuge"
    echo "bwa: $("$BWA" 2>&1 | grep -i '^Version' || true)"
    echo "samtools: $("$SAMTOOLS" --version | head -n 2 | tr '\n' ' ')"
    echo "fastp: $("$FASTP" --version 2>&1 | head -n 1)"
    echo
    echo "== Eingaben"
    echo "ref_fasta=$REF_FASTA bytes=$(stat -L -f %z "$REF_FASTA")"
    echo "r1=$R1 bytes=$(stat -L -f %z "$R1")"
    echo "r2=$R2 bytes=$(stat -L -f %z "$R2")"
    case "$R1" in
        *.gz) echo "r1_first_header=$(gzip -dc "$R1" 2>/dev/null | head -n 1 || true)" ;;
        *)    echo "r1_first_header=$(head -n 1 "$R1")" ;;
    esac
    echo
    echo "== Befehle"
    echo "a) $(cmdline "$BWA" index "$REF")   (nur wenn Indexdateien fehlen)"
    echo "b) $(cmdline "${CMD_FASTP[@]}")"
    echo "c) ulimit -n 8192; $(cmdline "${CMD_BWA[@]}") | $(cmdline "${CMD_FIXMATE[@]}") | $(cmdline "${CMD_SORT[@]}") | $(cmdline "${CMD_MARKDUP[@]}")"
    echo "d) $(cmdline "${CMD_INDEX[@]}") & $(cmdline "${CMD_FLAGSTAT[@]}") > $(printf '%q' "$OUT/${PROBE}_flagstat.txt")"
    echo "Jede Stufe umhuellt mit: $TIME -l -o $LOGS/<stufe>.time; Werkzeug-stderr in $LOGS/<stufe>.stderr"
} > "$INFO" 2>&1

# ---------- a) Referenz und bwa-Index ----------
step_begin a_ref_index "$LOGS/bwa_index.stderr"
ln -sfn "$REF_FASTA" "$REF"
if [ -e "$REF_FASTA.fai" ]; then ln -sfn "$REF_FASTA.fai" "$REF.fai"; fi
missing=0
for ext in amb ann bwt pac sa; do
    if [ ! -s "$REF.$ext" ]; then missing=1; fi
done
if [ "$missing" -eq 1 ]; then
    "$TIME" -l -o "$LOGS/bwa_index.time" "$BWA" index "$REF" 2> "$LOGS/bwa_index.stderr"
    step_end built
else
    step_end skipped_index_present
fi

# ---------- b) fastp ----------
step_begin b_fastp "$LOGS/fastp.stderr"
"$TIME" -l -o "$LOGS/fastp.time" "${CMD_FASTP[@]}" 2> "$LOGS/fastp.stderr"
step_end

# ---------- c) bwa mem | fixmate | sort | markdup ----------
step_begin c_pipe "$LOGS/bwa_mem.stderr $LOGS/fixmate.stderr $LOGS/sort.stderr $LOGS/markdup.stderr"
ulimit -n 8192
"$TIME" -l -o "$LOGS/bwa_mem.time" "${CMD_BWA[@]}" 2> "$LOGS/bwa_mem.stderr" \
    | "$TIME" -l -o "$LOGS/fixmate.time" "${CMD_FIXMATE[@]}" 2> "$LOGS/fixmate.stderr" \
    | "$TIME" -l -o "$LOGS/sort.time" "${CMD_SORT[@]}" 2> "$LOGS/sort.stderr" \
    | "$TIME" -l -o "$LOGS/markdup.time" "${CMD_MARKDUP[@]}" 2> "$LOGS/markdup.stderr"
step_end

# ---------- d) index und flagstat parallel ----------
step_begin d_index_flagstat "$LOGS/index.stderr $LOGS/flagstat.stderr"
"$TIME" -l -o "$LOGS/index.time" "${CMD_INDEX[@]}" 2> "$LOGS/index.stderr" &
P_INDEX=$!
"$TIME" -l -o "$LOGS/flagstat.time" "${CMD_FLAGSTAT[@]}" > "$OUT/${PROBE}_flagstat.txt" 2> "$LOGS/flagstat.stderr" &
P_FLAGSTAT=$!
wait "$P_INDEX"
wait "$P_FLAGSTAT"
step_end

RUN_T1=$(now)
TOTAL=$(secs "$RUN_T0" "$RUN_T1")

# ---------- Abschluss (nicht Teil der Messung) ----------
CURRENT_STEP=finish
STEP_ERR_FILES=""
check_volume
"$SAMTOOLS" quickcheck -v "$BAM"
rm -f "$T1" "$T2"
rm -f "$SORT_PREFIX".*
disk_snapshot after_cleanup
{
    echo "total_seconds=$TOTAL"
    echo "finished=$(iso)"
} > "$WORKDIR/DONE"
echo "end=$(iso)" >> "$INFO"
