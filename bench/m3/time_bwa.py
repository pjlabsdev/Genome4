#!/usr/bin/env python3
"""Laufzeit von bwa mem (pybwa) auf einem Kern, zum Vergleich mit `genome index sam --kerne 1`.

Aufruf:
    time_bwa.py <bwa-index-prefix> <reads_1.fq> [<reads_2.fq>] --out ergebnis.json

Gemessen wird nur der Aufruf des Aligners, nicht das Lesen der FASTQ. In der Zeit steckt, dass
pybwa je Zeile ein pysam-Objekt anlegt; das Programm bwa schreibt stattdessen Text. Der Vergleich
ist deshalb eine Größenordnung, keine genaue Messung.
"""
import argparse
import importlib.metadata
import json
import time

import pybwa
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("prefix")
ap.add_argument("fastq")
ap.add_argument("fastq2", nargs="?")
ap.add_argument("--block-basen", type=int, default=10_000_000)
ap.add_argument("--out", required=True)
args = ap.parse_args()

if args.fastq2:
    queries = []
    for a, b in zip(pysam.FastxFile(args.fastq), pysam.FastxFile(args.fastq2)):
        queries += [a, b]
else:
    queries = list(pysam.FastxFile(args.fastq))

opt = pybwa.BwaMemOptions(
    interleaved_paired_end=bool(args.fastq2), chunk_size=args.block_basen, threads=1
)
with pybwa.BwaMem(prefix=args.prefix) as aligner:
    start = time.time()
    results = aligner.align(queries, opt)
    seconds = time.time() - start

result = {
    "programm": "bwa mem aus pybwa " + importlib.metadata.version("pybwa"),
    "paare": bool(args.fastq2),
    "reads": len(queries),
    "kerne": 1,
    "sekunden": round(seconds, 1),
    "reads_je_sekunde": round(len(queries) / seconds),
}
print(result)
with open(args.out, "w") as f:
    json.dump(result, f, indent=2, ensure_ascii=False)
    f.write("\n")
