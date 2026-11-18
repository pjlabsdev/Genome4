#!/usr/bin/env python3
"""Vergleicht unsere SAM-Ausgabe für Read-Paare Zeile für Zeile mit bwa mem (pybwa, Vorgaben).

Je Read werden alle Zeilen verglichen, Feld für Feld: Flag, Contig, Position, MAPQ, CIGAR,
Contig und Position des Partners, TLEN, Basen, Qualitäten und die Angaben NM, MD, MC, MQ, AS,
XS, SA, XA. Zuordnung nach Reihenfolge. pybwa liefert die Flag-Bits 2 (richtiges Paar), 64 und
128 (Read 1, Read 2) nicht; sie bleiben deshalb aus dem Vergleich heraus. Beide Seiten schätzen die Insertgröße je Block von
`block_basen` Basen (bwa -K).
Aufruf: compare_bwa_pairs.py bwa-prefix r1.fq r2.fq unser.sam ergebnis.json [block_basen]
"""
import json
import sys

import pybwa
import pysam

prefix, fq1, fq2, sam_path, out_path = sys.argv[1:6]
block_bases = int(sys.argv[6]) if len(sys.argv) > 6 else 10_000_000
TAGS = ("NM", "MD", "MC", "MQ", "AS", "XS", "SA", "XA")
FIELDS = ("flag", "contig", "pos", "mapq", "cigar", "partner_contig", "partner_pos", "tlen", "seq", "qual") + TAGS


def ours_by_read(path):
    group, key = [], None
    for line in open(path):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        tags = {t[:2]: t[5:] for t in f[11:]}
        flag = int(f[1])
        rec = (flag & ~0xC2, f[2], int(f[3]), int(f[4]), f[5], f[6], int(f[7]), int(f[8]), f[9], f[10]) + tuple(tags.get(t) for t in TAGS)
        k = (f[0], flag & 0xC0)
        if key is not None and k != key:
            yield key[0], group
            group = []
        key = k
        group.append(rec)
    if group:
        yield key[0], group


def bwa_rec(r):
    tags = dict(r.get_tags())
    qual = pysam.qualities_to_qualitystring(r.query_qualities) if r.query_qualities is not None else "*"
    mate = "*" if r.next_reference_name is None else ("=" if r.next_reference_name == r.reference_name else r.next_reference_name)
    return (r.flag & ~0xC2, r.reference_name or "*", r.reference_start + 1, r.mapping_quality, r.cigarstring or "*",
            mate, r.next_reference_start + 1, r.template_length,
            r.query_sequence or "*", qual) + tuple(None if tags.get(t) is None else str(tags[t]) for t in TAGS)


queries = []
for a, b in zip(pysam.FastxFile(fq1), pysam.FastxFile(fq2)):
    queries += [a, b]
# pybwa teilt selbst in Blöcke: ein Read, dann weitere, bis die Basenzahl erreicht und die Zahl
# der Reads gerade ist. Bei gleich langen Reads sind das dieselben Blöcke wie bei uns (Paare,
# bis die Basenzahl erreicht ist). Die laufende Nummer der Reads zählt über die Blöcke weiter;
# sie entscheidet zwischen gleich guten Orten.
opt = pybwa.BwaMemOptions(interleaved_paired_end=True, chunk_size=block_bases)
with pybwa.BwaMem(prefix=prefix) as aligner:
    results = aligner.align(queries, opt)

stats = dict(reads=0, paare=0, reads_alle_zeilen_gleich=0, paare_ganz_gleich=0, reads_andere_zeilenzahl=0,
             zeilen_bwa=0, zeilen_wir=0, reads_ohne_ort_bwa=0,
             reads_mit_ergaenzenden_zeilen_bwa=0)
diff_fields, examples = {}, []
ours = ours_by_read(sam_path)
pair_ok = True
for n, (q, recs) in enumerate(zip(queries, results)):
    name, mine = next(ours)
    assert name == q.name, f"Reihenfolge weicht ab: {name} / {q.name}"
    theirs = [bwa_rec(r) for r in recs]
    stats["reads"] += 1
    stats["zeilen_bwa"] += len(theirs)
    stats["zeilen_wir"] += len(mine)
    stats["reads_ohne_ort_bwa"] += int(theirs[0][0] & 4 != 0)
    stats["reads_mit_ergaenzenden_zeilen_bwa"] += int(len(theirs) > 1)
    same = theirs == mine
    if n % 2 == 0:
        pair_ok = same
    else:
        stats["paare"] += 1
        stats["paare_ganz_gleich"] += int(pair_ok and same)
    if same:
        stats["reads_alle_zeilen_gleich"] += 1
        continue
    if len(theirs) != len(mine):
        stats["reads_andere_zeilenzahl"] += 1
        fields = ["zeilenzahl"]
    else:
        fields = sorted({FIELDS[i] for a, b in zip(theirs, mine) for i in range(len(FIELDS)) if a[i] != b[i]})
    key = "+".join(fields)
    diff_fields[key] = diff_fields.get(key, 0) + 1
    if len(examples) < 80:
        examples.append(dict(read=q.name, nr=n, felder=fields, bwa=theirs, wir=mine))
stats["abweichende_felder"] = dict(sorted(diff_fields.items(), key=lambda x: -x[1]))
json.dump(dict(block_basen=block_bases, ergebnis=stats, beispiele=examples), open(out_path, "w"), indent=1, ensure_ascii=False)
print(json.dumps(stats, ensure_ascii=False))
