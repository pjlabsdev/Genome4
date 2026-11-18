#!/usr/bin/env python3
"""Vergleicht den primären Bereich je Read mit dem primären Alignment von bwa mem.

bwa läuft über pybwa (eingebautes bwa) mit Vorgaben, ein Read nach dem anderen (single-end).
Unsere Tabelle stammt aus `genome index ketten --tabelle`.
Aufruf: compare_bwa_align.py bwa-prefix reads.fq unser.tsv ergebnis.json
"""
import json
import sys
import time

import pybwa
import pysam

prefix, fastq, ours_path, out_path = sys.argv[1:5]
ours = {}
for line in open(ours_path):
    name, contig, pos, strand, mapq, score, truesc, qb, qe = line.rstrip("\n").split("\t")[:9]
    ours[name] = (contig, int(pos), strand, int(mapq), int(score), int(truesc), int(qb), int(qe))

queries = list(pysam.FastxFile(fastq))
t = time.time()
with pybwa.BwaMem(prefix=prefix) as aligner:
    results = aligner.align(queries)
bwa_seconds = time.time() - t

stats = dict(reads=0, gleich_ort_strang=0, gleich_auch_mapq=0, gleich_auch_score=0,
             anderer_ort_beide_mapq0=0, anderer_ort_sonst=0, nur_wir_ohne=0, nur_bwa_ohne=0,
             mapq_abweichung_bei_gleichem_ort=0, score_abweichung_bei_gleichem_ort=0)
mapq_diffs, examples = {}, []
for q, recs in zip(queries, results):
    stats["reads"] += 1
    prim = [r for r in recs if not r.is_secondary and not r.is_supplementary]
    b = prim[0] if prim else None
    o = ours.get(q.name)
    b_un = b is None or b.is_unmapped
    o_un = o is None or o[0] == "*"
    if b_un or o_un:
        stats["nur_wir_ohne"] += int(o_un and not b_un)
        stats["nur_bwa_ohne"] += int(b_un and not o_un)
        continue
    b_strand = "-" if b.is_reverse else "+"
    same = b.reference_name == o[0] and b.reference_start + 1 == o[1] and b_strand == o[2]
    if same:
        stats["gleich_ort_strang"] += 1
        if b.mapping_quality == o[3]:
            stats["gleich_auch_mapq"] += 1
        else:
            stats["mapq_abweichung_bei_gleichem_ort"] += 1
            d = o[3] - b.mapping_quality
            mapq_diffs[d] = mapq_diffs.get(d, 0) + 1
        if b.get_tag("AS") == o[4]:
            stats["gleich_auch_score"] += 1
        else:
            stats["score_abweichung_bei_gleichem_ort"] += 1
            if len(examples) < 20:
                examples.append(dict(read=q.name, art="score", bwa_as=b.get_tag("AS"), wir=o[4], cigar=b.cigarstring))
    elif b.mapping_quality == 0 and o[3] == 0:
        stats["anderer_ort_beide_mapq0"] += 1
    else:
        stats["anderer_ort_sonst"] += 1
        if len(examples) < 40:
            examples.append(dict(read=q.name, art="ort", bwa=[b.reference_name, b.reference_start + 1, b_strand,
                                 b.mapping_quality, b.get_tag("AS"), b.cigarstring], wir=list(o)))
stats["bwa_sekunden"] = round(bwa_seconds, 1)
stats["mapq_differenzen_wir_minus_bwa"] = dict(sorted(mapq_diffs.items()))
json.dump(dict(ergebnis=stats, beispiele=examples), open(out_path, "w"), indent=1, ensure_ascii=False)
print(json.dumps(stats, ensure_ascii=False))
