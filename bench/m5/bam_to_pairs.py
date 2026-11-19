#!/usr/bin/env python3
"""Echte Read-Paare aus einer BAM als zwei FASTQ-Dateien, in der Reihenfolge der BAM.

Nur Hauptzeilen. Ein Paar kommt nur heraus, wenn beide Reads in der Datei stehen (bei einem
Ausschnitt fehlen Partner außerhalb). Reads auf dem Gegenstrang werden zurückgedreht, so dass
Basen und Qualitäten wie im Sequenziergerät stehen.
Aufruf: bam_to_pairs.py ein.bam r1.fq r2.fq
"""
import sys
import pysam

COMP = str.maketrans("ACGTN", "TGCAN")
first, second, order = {}, {}, []
for r in pysam.AlignmentFile(sys.argv[1], check_sq=False):
    if r.is_secondary or r.is_supplementary or not r.is_paired or r.query_sequence is None:
        continue
    seq, qual = r.query_sequence, pysam.qualities_to_qualitystring(r.query_qualities)
    if r.is_reverse:
        seq, qual = seq.translate(COMP)[::-1], qual[::-1]
    side = first if r.is_read1 else second
    if r.query_name in side:
        continue
    side[r.query_name] = (seq, qual)
    other = second if r.is_read1 else first
    if r.query_name in other:
        order.append(r.query_name)
with open(sys.argv[2], "w") as f1, open(sys.argv[3], "w") as f2:
    for name in order:
        for f, side in ((f1, first), (f2, second)):
            seq, qual = side[name]
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")
print(len(order), "Paare")
