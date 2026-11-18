#!/usr/bin/env python3
"""Erklärt Reads, bei denen kein Seed am wahren Ort liegt.

Für jeden solchen Read: die fehlerfreien Abschnitte ab 19 Basen am wahren Ort bestimmen und
prüfen, ob sich der Abschnitt im Read (also mit der fehlerhaften Nachbarbase) verlängern lässt
und anderswo im Genom vorkommt. Dann ist der Abschnitt kein längster Treffer des Reads, und die
Seed-Suche meldet zu Recht den längeren Treffer, der nicht am wahren Ort liegt.
Aufruf: explain_missing.py genome index.gcfm ref.fa fehlende.tsv
"""
import subprocess
import sys
import tempfile

genome, index, fa, missing = sys.argv[1:5]
seq = "".join(l.strip().upper() for l in open(fa) if not l.startswith(">"))
comp = str.maketrans("ACGT", "TGCA")
queries, plan = [], []
for line in open(missing):
    name, read = line.split()
    _, pos, strand, _ = name.split(":")
    ref = seq[int(pos) - 1:int(pos) - 1 + len(read)]
    if strand == "-":
        ref = ref.translate(comp)[::-1]
    runs, start = [], None
    for i in range(len(read) + 1):
        same = i < len(read) and read[i] == ref[i]
        if same and start is None:
            start = i
        if not same and start is not None:
            if i - start >= 19:
                runs.append((start, i))
            start = None
    entry = []
    for a, b in runs:
        ext = []
        if a > 0:
            ext.append(read[a - 1:b])
        if b < len(read):
            ext.append(read[a:b + 1])
        entry.append((a, b, len(queries), len(ext)))
        queries.append(read[a:b])
        queries.extend(ext)
    plan.append((name, entry))
with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
    f.write("\n".join(queries) + "\n")
out = subprocess.run([genome, "index", "count", index, f.name], capture_output=True, text=True, check=True).stdout
counts = [int(l.split("\t")[1]) for l in out.splitlines()]
explained = thinned = other = 0
for name, entry in plan:
    # erklärt, wenn jeder fehlerfreie Abschnitt im Read verlängerbar ist (kein längster Treffer)
    all_ext = all(any(counts[q + 1 + k] > 0 for k in range(n)) for _, _, q, n in entry)
    many = any(counts[q] > 500 for _, _, q, _ in entry)
    if entry and all_ext:
        explained += 1
    elif many:
        thinned += 1
    else:
        other += 1
        print("unerklärt", name, [(a, b, counts[q]) for a, b, q, _ in entry])
print(f"reads={len(plan)} verlaengerbar={explained} ausgeduennt={thinned} unerklaert={other}")
