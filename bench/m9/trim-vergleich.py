#!/usr/bin/env python3
"""Vergleicht die Ausgabe von fastp 1.4.0 mit der von `genome trim`, Read für Read.

Aufruf: trim-vergleich.py fastp_R1 fastp_R2 genome_R1 genome_R2 [fastp.json genome.log.jsonl]
Dateien roh oder gzip. Gibt JSON aus: Paare je Seite, nur bei fastp, nur bei genome, gleiche
Folge und Qualität, abweichende Länge, und bis zu 20 Beispiele je Art.
"""
import gzip
import json
import sys


def records(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as f:
        while True:
            head = f.readline()
            if not head:
                return
            seq = f.readline().rstrip(b"\r\n")
            f.readline()
            qual = f.readline().rstrip(b"\r\n")
            name = head[1:].split()[0]
            if name.endswith(b"/1") or name.endswith(b"/2"):
                name = name[:-2]
            yield name.decode(), seq, qual


def compare(fastp_path, genome_path):
    a = {n: (s, q) for n, s, q in records(fastp_path)}
    out = {"fastp": len(a), "genome": 0, "nur_genome": 0, "gleich": 0,
           "andere_laenge": 0, "gleiche_laenge_andere_folge": 0, "beispiele": []}
    seen = set()
    for n, s, q in records(genome_path):
        out["genome"] += 1
        seen.add(n)
        if n not in a:
            out["nur_genome"] += 1
            if len(out["beispiele"]) < 20:
                out["beispiele"].append({"art": "nur_genome", "name": n, "laenge": len(s)})
            continue
        fs, fq = a[n]
        if fs == s and fq == q:
            out["gleich"] += 1
        elif len(fs) != len(s):
            out["andere_laenge"] += 1
            if len(out["beispiele"]) < 20:
                out["beispiele"].append({"art": "andere_laenge", "name": n,
                                         "fastp": len(fs), "genome": len(s)})
        else:
            out["gleiche_laenge_andere_folge"] += 1
    missing = [n for n in a if n not in seen]
    out["nur_fastp"] = len(missing)
    for n in missing[:20]:
        out["beispiele"].append({"art": "nur_fastp", "name": n, "laenge": len(a[n][0])})
    return out


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    result = {"read1": compare(sys.argv[1], sys.argv[3]),
              "read2": compare(sys.argv[2], sys.argv[4])}
    if len(sys.argv) >= 7:
        fj = json.load(open(sys.argv[5]))
        result["fastp_adapter"] = fj.get("adapter_cutting", {})
        result["fastp_filter"] = fj.get("filtering_result", {})
        result["genome_entscheidungen"] = [
            e for e in map(json.loads, open(sys.argv[6]))
            if str(e.get("code", "")).startswith("TRIM.")
        ]
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
