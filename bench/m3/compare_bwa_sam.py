#!/usr/bin/env python3
"""Vergleicht unsere SAM-Ausgabe Zeile für Zeile mit bwa mem (pybwa, Vorgaben, einzelne Reads).

Je Read werden alle Zeilen verglichen (primär und ergänzend), Feld für Feld: Flag, Contig,
Position, MAPQ, CIGAR, Basen, Qualitäten und die Angaben NM, MD, AS, XS, SA, XA.
Zuordnung nach Reihenfolge der Reads, nicht nach Namen.
Aufruf: compare_bwa_sam.py bwa-prefix reads.fq unser.sam ergebnis.json
"""
import json
import sys

import pybwa
import pysam

prefix, fastq, sam_path, out_path = sys.argv[1:5]
TAGS = ("NM", "MD", "AS", "XS", "SA", "XA")


def ours_by_read(path):
    group, last = [], None
    for line in open(path):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        tags = {t[:2]: t[5:] for t in f[11:]}
        rec = (int(f[1]), f[2], int(f[3]), int(f[4]), f[5], f[9], f[10]) + tuple(tags.get(t) for t in TAGS)
        # Zeilen eines Reads stehen hintereinander; die erste ist nie ergänzend.
        if last is not None and not (rec[0] & 0x800):
            yield last, group
            group = []
        last = f[0]
        group.append(rec)
    if group:
        yield last, group


def bwa_rec(r):
    tags = dict(r.get_tags())
    cigar = r.cigarstring or "*"
    qual = pysam.qualities_to_qualitystring(r.query_qualities) if r.query_qualities is not None else "*"
    return (r.flag, r.reference_name or "*", r.reference_start + 1, r.mapping_quality, cigar,
            r.query_sequence or "*", qual) + tuple(None if tags.get(t) is None else str(tags[t]) for t in TAGS)


queries = list(pysam.FastxFile(fastq))
with pybwa.BwaMem(prefix=prefix) as aligner:
    results = aligner.align(queries)

FIELDS = ("flag", "contig", "pos", "mapq", "cigar", "seq", "qual") + TAGS
stats = dict(reads=0, reads_alle_zeilen_gleich=0, reads_andere_zeilenzahl=0, zeilen_bwa=0, zeilen_wir=0,
             reads_mit_ergaenzenden_zeilen_bwa=0, reads_mit_xa_bwa=0, reads_ohne_ort_bwa=0)
diff_fields, examples = {}, []
ours = ours_by_read(sam_path)
for q, recs in zip(queries, results):
    name, mine = next(ours)
    assert name == q.name, f"Reihenfolge weicht ab: {name} / {q.name}"
    theirs = [bwa_rec(r) for r in recs]
    stats["reads"] += 1
    stats["zeilen_bwa"] += len(theirs)
    stats["zeilen_wir"] += len(mine)
    stats["reads_mit_ergaenzenden_zeilen_bwa"] += int(len(theirs) > 1)
    stats["reads_mit_xa_bwa"] += int(any(t[-1] is not None for t in theirs))
    stats["reads_ohne_ort_bwa"] += int(theirs[0][0] & 4 != 0)
    if theirs == mine:
        stats["reads_alle_zeilen_gleich"] += 1
        continue
    if len(theirs) != len(mine):
        stats["reads_andere_zeilenzahl"] += 1
        fields = ["zeilenzahl"]
    else:
        fields = sorted({FIELDS[i] for a, b in zip(theirs, mine) for i in range(len(FIELDS)) if a[i] != b[i]})
    key = "+".join(fields)
    diff_fields[key] = diff_fields.get(key, 0) + 1
    if len(examples) < 60:
        examples.append(dict(read=q.name, felder=fields, bwa=theirs, wir=mine))
stats["abweichende_felder"] = dict(sorted(diff_fields.items(), key=lambda x: -x[1]))
json.dump(dict(ergebnis=stats, beispiele=examples), open(out_path, "w"), indent=1, ensure_ascii=False)
print(json.dumps(stats, ensure_ascii=False))
