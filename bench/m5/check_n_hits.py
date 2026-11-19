#!/usr/bin/env python3
"""Prüft, ob ein Read in der Referenz von bwa eine Stelle trifft, die in der FASTA aus N besteht.

bwa füllt N-Strecken beim Indexbau mit Zufallsbasen (`.pac`). Ein Stück eines Reads kann dort
zufällig passen; bwa meldet das dann als zweitbesten Treffer (`XS`). Gesucht werden alle
Stücke von `k` Basen des Reads auf beiden Strängen im Text von bwa, gemeldet werden die, die
in der FASTA mindestens ein N berühren.
Aufruf: check_n_hits.py ref.fa reads.txt [k]     (ref.fa.pac daneben; reads.txt: Name Basen)
"""
import sys

fa, reads = sys.argv[1], sys.argv[2]
k = int(sys.argv[3]) if len(sys.argv) > 3 else 19
ref = "".join(line.strip() for line in open(fa) if not line.startswith(">")).upper()
table = {i: "".join("ACGT"[(i >> s) & 3] for s in (6, 4, 2, 0)) for i in range(256)}
pac = "".join(table[x] for x in open(fa + ".pac", "rb").read())[:len(ref)]
comp = str.maketrans("ACGT", "TGCA")
explained = total = 0
for line in open(reads):
    name, seq = line.split()
    hits = []
    for x in range(len(seq) - k + 1):
        piece = seq[x:x + k]
        if "N" in piece:
            continue
        for strand, q in (("+", piece), ("-", piece.translate(comp)[::-1])):
            p = pac.find(q)
            while p >= 0:
                if "N" in ref[p:p + k]:
                    hits.append((x, strand, p + 1))
                p = pac.find(q, p + 1)
    total += 1
    explained += bool(hits)
    print(name, "Treffer in aufgefüllten N-Strecken:", len(hits), hits[:3])
print(f"{explained} von {total} Reads treffen eine aufgefüllte N-Strecke")
