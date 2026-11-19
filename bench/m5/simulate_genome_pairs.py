#!/usr/bin/env python3
"""Lesepaare über alle Contigs einer (auch gepackten) FASTA ziehen, gewichtet nach Länge.

Stücke mit normalverteilter Länge, von beiden Enden gelesen (Read 1 und Read 2 auf
entgegengesetzten Strängen). Nur Substitutionen mit `fehlerrate` je Base, Stücke mit N
werden verworfen. Braucht die ganze Referenz im RAM (ein Byte je Base).

Name: contig:pos1:pos2:strang:p:nummer mit den linken Enden beider Reads auf dem
Vorwärtsstrang (1-basiert) und dem Strang von Read 1.
Aufruf: simulate_genome_pairs.py ref.fa[.gz] anzahl laenge mittel streuung fehlerrate saat \
        r1.fq r2.fq
"""
import bisect
import gzip
import random
import sys

fa, n, length = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
mean, sd, rate = map(float, sys.argv[4:7])
seed, out1, out2 = int(sys.argv[7]), sys.argv[8], sys.argv[9]
rng = random.Random(seed)

names, seqs, parts = [], [], []
opener = gzip.open if fa.endswith(".gz") else open
with opener(fa, "rb") as src:
    for line in src:
        if line.startswith(b">"):
            if parts:
                seqs.append(b"".join(parts).upper())
                parts = []
            names.append(line[1:].split()[0].decode())
        else:
            parts.append(line.rstrip())
seqs.append(b"".join(parts).upper())
del parts
ends, total = [], 0
for s in seqs:
    total += len(s)
    ends.append(total)

COMP = bytes.maketrans(b"ACGT", b"TGCA")
OTHER = {65: b"CGT", 67: b"AGT", 71: b"ACT", 84: b"ACG"}


def with_errors(read):
    if rate <= 0:
        return read
    out, i = bytearray(read), int(rng.expovariate(rate))
    while i < len(out):
        out[i] = rng.choice(OTHER[out[i]])
        i += 1 + int(rng.expovariate(rate))
    return bytes(out)


qual = b"F" * length
done = 0
with open(out1, "wb") as f1, open(out2, "wb") as f2:
    while done < n:
        size = max(length, int(rng.gauss(mean, sd) + 0.5))
        g = rng.randrange(total)
        c = bisect.bisect_right(ends, g)
        seq = seqs[c]
        pos = g - (ends[c] - len(seq))
        if pos + size > len(seq):
            continue
        piece = seq[pos:pos + size]
        if piece.strip(b"ACGT"):
            if b"N" in piece or piece.translate(None, b"ACGT"):
                continue
        left, right = piece[:length], piece[-length:].translate(COMP)[::-1]
        forward = rng.random() < 0.5
        r1, r2 = (left, right) if forward else (right, left)
        p1, p2 = (pos + 1, pos + size - length + 1) if forward else (pos + size - length + 1, pos + 1)
        head = f"@{names[c]}:{p1}:{p2}:{'+' if forward else '-'}:p:{done}".encode()
        f1.write(head + b"/1\n" + with_errors(r1) + b"\n+\n" + qual + b"\n")
        f2.write(head + b"/2\n" + with_errors(r2) + b"\n+\n" + qual + b"\n")
        done += 1
