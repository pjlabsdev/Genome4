#!/usr/bin/env python3
"""Mehrere Kerne und BAM: Laufzeit je Kernzahl, Gleichheit der Ausgabe, BAM gegen SAM.

Aufruf:
    measure_cores.py <genome> <index.gcfm> <reads_1.fq> [<reads_2.fq>] --out ergebnis.json

Für 1, 2 und 4 Kerne wird `genome index sam` gestartet. Verglichen wird die Ausgabe ohne die
`@PG`-Zeile (dort steht der Aufruf). Danach ein Lauf mit `--bam`; die BAM wird mit samtools
(pysam) nach SAM gewandelt und mit der SAM-Ausgabe verglichen.
"""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
import time

import pysam


def records(text):
    return [line for line in text.splitlines() if not line.startswith("@PG")]


def digest(lines):
    h = hashlib.sha256()
    for line in lines:
        h.update(line.encode())
        h.update(b"\n")
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genome")
    ap.add_argument("index")
    ap.add_argument("fastq")
    ap.add_argument("fastq2", nargs="?")
    ap.add_argument("--kerne", default="1,2,4")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    base = [args.genome, "index", "sam", args.index, "--fastq", args.fastq]
    if args.fastq2:
        base += ["--fastq2", args.fastq2]

    runs = []
    reference = None
    for cores in [int(x) for x in args.kerne.split(",")]:
        start = time.time()
        out = subprocess.run(
            base + ["--kerne", str(cores)], check=True, capture_output=True, text=True
        ).stdout
        seconds = time.time() - start
        lines = records(out)
        if reference is None:
            reference = lines
        runs.append(
            {
                "kerne": cores,
                "sekunden": round(seconds, 1),
                "zeilen": sum(1 for l in lines if not l.startswith("@")),
                "sha256_ohne_pg": digest(lines),
                "gleich_wie_erster_lauf": lines == reference,
            }
        )
        print(runs[-1], flush=True)

    with tempfile.TemporaryDirectory() as tmp:
        bam = os.path.join(tmp, "aus.bam")
        cores = runs[-1]["kerne"]
        start = time.time()
        subprocess.run(
            base + ["--kerne", str(cores), "--bam", bam], check=True, capture_output=True
        )
        seconds = time.time() - start
        size = os.path.getsize(bam)
        pysam.quickcheck(bam)
        back = records(pysam.view("-h", bam))
    bam_result = {
        "kerne": cores,
        "sekunden": round(seconds, 1),
        "bytes": size,
        "samtools_quickcheck": "bestanden",
        "zeilen_gleich_wie_sam": back == reference,
        "abweichende_zeilen": sum(1 for a, b in zip(back, reference) if a != b)
        + abs(len(back) - len(reference)),
    }
    print(bam_result, flush=True)

    result = {
        "eingabe": [os.path.basename(args.fastq)]
        + ([os.path.basename(args.fastq2)] if args.fastq2 else []),
        "paare": bool(args.fastq2),
        "kerne_im_rechner": os.cpu_count(),
        "laeufe": runs,
        "bam": bam_result,
    }
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
        f.write("\n")


if __name__ == "__main__":
    main()
