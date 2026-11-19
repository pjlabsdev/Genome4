import collections, json, sys, pysam
bam = pysam.AlignmentFile(sys.argv[1])
c = collections.Counter(); prev = (-1, -1); order = True
for r in bam:
    key = (r.reference_id if r.reference_id >= 0 else 1 << 30, r.reference_start)
    if key < prev: order = False
    prev = key
    if r.is_secondary or r.is_supplementary: c["neben"] += 1; continue
    contig, p1, p2 = r.query_name.split(":")[:3]
    want = int(p1 if r.is_read1 else p2)
    q = r.mapping_quality
    band = "mapq0" if q == 0 else "mapq1_29" if q < 30 else "mapq30plus"
    c[band] += 1
    # weiche Enden zurückrechnen
    start = r.reference_start + 1 - (r.cigartuples[0][1] if r.cigartuples and r.cigartuples[0][0] == 4 else 0)
    ok = (not r.is_unmapped) and r.reference_name == contig and abs(start - want) <= 5
    c[band + "_richtig"] += ok
    c["gesamt"] += 1; c["richtig"] += ok
json.dump({"sortiert": order, **c}, sys.stdout, indent=1, sort_keys=True)
