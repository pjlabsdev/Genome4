#!/usr/bin/env python3
"""Vorher und nachher: zwei gebaute Programme auf denselben Paaren, je mehrere Läufe.

Je Lauf die Sekunden des Ausrichtens aus dem Protokoll (ALIGN.FERTIG, ohne Laden des Index) und
die Prüfsumme der SAM-Ausgabe ohne die Kopfzeile @PG. Gleiche Prüfsumme heißt: jede Zeile gleich.

Aufruf:
  measure_tempo.py <vorher> <nachher> <index.gcfm> <reads_1.fq> <reads_2.fq> <kerne> [läufe]
"""
import hashlib
import json
import statistics
import subprocess
import sys
import tempfile


def run(binary, index, fq1, fq2, cores):
    with tempfile.NamedTemporaryFile(suffix=".jsonl") as log:
        cmd = [binary, "index", "sam", index, "--fastq", fq1, "--fastq2", fq2,
               "--kerne", str(cores), "--log", log.name]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        md5 = hashlib.md5()
        for line in proc.stdout:
            if not line.startswith(b"@PG"):
                md5.update(line)
        if proc.wait() != 0:
            sys.exit(f"Lauf gescheitert: {' '.join(cmd)}")
        seconds = None
        for line in open(log.name, encoding="utf-8"):
            event = json.loads(line)
            if event["code"] == "ALIGN.FERTIG":
                seconds = event["fields"]["sekunden"]
        return seconds, md5.hexdigest()


def main():
    before, after, index, fq1, fq2, cores = sys.argv[1:7]
    runs = int(sys.argv[7]) if len(sys.argv) > 7 else 3
    out = {"index": index, "reads": [fq1, fq2], "kerne": int(cores), "laeufe": runs}
    sums = set()
    for name, binary in (("vorher", before), ("nachher", after)):
        results = [run(binary, index, fq1, fq2, cores) for _ in range(runs)]
        seconds = [round(s, 3) for s, _ in results]
        sums.update(m for _, m in results)
        out[name] = {"sekunden": seconds, "median": statistics.median(seconds)}
    out["ausgabe_gleich"] = len(sums) == 1
    out["md5_ohne_pg"] = sorted(sums)
    out["zeit_nachher_zu_vorher"] = round(out["nachher"]["median"] / out["vorher"]["median"], 3)
    json.dump(out, sys.stdout, ensure_ascii=False, indent=2)
    print()


if __name__ == "__main__":
    main()
