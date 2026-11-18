#!/usr/bin/env python3
"""Einfache Reads aus einer FASTA ziehen: feste Länge, beide Stränge, Substitutionsfehler.

Nur für die Prüfung der Seed-Suche. Kein realistisches Fehlermodell (keine Indels,
keine Qualitätswerte); der wahre Ort steht im Namen: contig:pos(1-basiert):strang:fehler.
Aufruf: simulate_reads.py ref.fa anzahl laenge fehlerrate saat > reads.fq
"""
import random
import sys

fa, n, length, rate, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
rng = random.Random(seed)
name, parts = None, []
for line in open(fa):
    if line.startswith(">"):
        if name is not None:
            break
        name = line[1:].split()[0]
    else:
        parts.append(line.strip().upper())
seq = "".join(parts)
comp = str.maketrans("ACGT", "TGCA")
out = sys.stdout
made = 0
while made < n:
    pos = rng.randrange(len(seq) - length)
    piece = seq[pos:pos + length]
    if "N" in piece:
        continue
    strand = "+"
    if rng.random() < 0.5:
        piece = piece.translate(comp)[::-1]
        strand = "-"
    bases = list(piece)
    errors = 0
    for i in range(length):
        if rng.random() < rate:
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
            errors += 1
    out.write(f"@{name}:{pos + 1}:{strand}:{errors}\n{''.join(bases)}\n+\n{'I' * length}\n")
    made += 1
