#!/usr/bin/env python3
"""Wertet SAM-Zeilen simulierter Paare gegen den wahren Ort aus (Namen aus simulate_pairs.py).

Gezählt werden die primären Zeilen der Reads mit bekanntem Ort (alle außer dem fremden Read 2
der Paare mit Art x). Richtig heißt: Contig stimmt und das linke Ende (ohne abgeschnittene
Basen gerechnet) liegt höchstens 20 Basen neben dem wahren. Einzeln ausgerichtete Reads
(zwei SAM-Dateien, Read 1 und Read 2) lassen sich genauso auswerten.
Aufruf: evaluate_pairs.py ergebnis.json paare.sam
        evaluate_pairs.py ergebnis.json read1.sam read2.sam
"""
import json
import re
import sys

out_path, paths = sys.argv[1], sys.argv[2:]
stats = dict(reads_mit_wahrem_ort=0, richtig=0, falsch=0, ohne_ort=0, falsch_mapq_ab_30=0, mapq_ab_30=0,
             richtig_mapq_0=0, fremde_reads=0, fremde_reads_mit_ort=0, richtiges_paar_flag=0)
by_kind = {}


def handle(f, read_nr):
    flag = int(f[1])
    if flag & 0x900:
        return
    contig, p1, p2, strand, kind, _ = f[0].rsplit(":", 5)
    if flag & 0x80:
        read_nr = 2
    elif flag & 0x40:
        read_nr = 1
    truth = int(p1 if read_nr == 1 else p2)
    if kind == "x" and read_nr == 2:
        stats["fremde_reads"] += 1
        stats["fremde_reads_mit_ort"] += int(not flag & 4)
        return
    k = by_kind.setdefault(kind, dict(reads=0, richtig=0, falsch=0, ohne_ort=0))
    stats["reads_mit_wahrem_ort"] += 1
    k["reads"] += 1
    stats["richtiges_paar_flag"] += int(flag & 2 != 0)
    if flag & 4:
        stats["ohne_ort"] += 1
        k["ohne_ort"] += 1
        return
    mapq = int(f[4])
    lead = re.match(r"(\d+)S", f[5])
    start = int(f[3]) - (int(lead.group(1)) if lead else 0)
    ok = f[2] == contig and abs(start - truth) <= 20
    stats["richtig" if ok else "falsch"] += 1
    k["richtig" if ok else "falsch"] += 1
    stats["mapq_ab_30"] += int(mapq >= 30)
    stats["falsch_mapq_ab_30"] += int(mapq >= 30 and not ok)
    stats["richtig_mapq_0"] += int(mapq == 0 and ok)


for nr, path in enumerate(paths):
    for line in open(path):
        if not line.startswith("@"):
            handle(line.rstrip("\n").split("\t"), nr + 1)
stats["nach_art"] = by_kind
json.dump(stats, open(out_path, "w"), indent=1, ensure_ascii=False)
print(json.dumps(stats, ensure_ascii=False))
