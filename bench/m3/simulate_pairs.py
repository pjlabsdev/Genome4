#!/usr/bin/env python3
"""Einfache Read-Paare aus einer FASTA ziehen: DNA-Stücke mit normalverteilter Länge, von
beiden Enden gelesen (Read 1 vorwärts, Read 2 vom anderen Ende auf dem Gegenstrang).

Kein realistisches Fehlermodell. Substitutionen je Base mit `fehlerrate`, dazu im Stück
Insertionen und Deletionen von 1 bis 10 Basen mit `indelrate`. Damit auch die Sonderfälle
vorkommen, ist bei `anteil_fremd` der Paare Read 2 eine Zufallsfolge (findet keinen Ort) und
bei `anteil_fern` stammt Read 2 aus einem anderen, unabhängig gezogenen Stück.

Name: contig:pos1:pos2:strang:art:nummer. pos1 und pos2 sind die linken Enden von Read 1 und
Read 2 auf dem Vorwärtsstrang (1-basiert, ohne Verschiebung durch Indels im Read), strang ist
der Strang von Read 1, art ist p (Paar), x (Read 2 fremd) oder f (Read 2 fern).
Aufruf: simulate_pairs.py ref.fa anzahl laenge mittel streuung fehlerrate indelrate \
        anteil_fremd anteil_fern saat r1.fq r2.fq
"""
import random
import sys

fa, n, length = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
mean, sd, rate, indel_rate, part_alien, part_far = map(float, sys.argv[4:10])
seed, out1, out2 = int(sys.argv[10]), sys.argv[11], sys.argv[12]
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


def fragment():
    """Ein Stück ohne N: (linkes Ende, Folge mit Indels, Länge auf der Referenz)."""
    while True:
        size = max(length, int(rng.gauss(mean, sd) + 0.5))
        pos = rng.randrange(len(seq) - size - 64)
        piece = seq[pos:pos + size + 64]
        if "N" in piece:
            continue
        bases, i = [], 0
        while i < size:
            if 5 <= i and rng.random() < indel_rate:
                k = rng.randint(1, 10)
                if rng.random() < 0.5:
                    i += k
                else:
                    bases.extend(rng.choice("ACGT") for _ in range(k))
                    continue
            bases.append(piece[i])
            i += 1
        if len(bases) < length:
            continue
        return pos, "".join(bases), i


def with_errors(read):
    bases = list(read)
    for i in range(len(bases)):
        if rng.random() < rate:
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


with open(out1, "w") as f1, open(out2, "w") as f2:
    for made in range(n):
        pos, frag, ref_len = fragment()
        left, right = frag[:length], frag[-length:].translate(comp)[::-1]
        pos_left, pos_right = pos + 1, pos + ref_len - length + 1
        kind = "p"
        x = rng.random()
        if x < part_alien:
            kind = "x"
            right = "".join(rng.choice("ACGT") for _ in range(length))
        elif x < part_alien + part_far:
            kind = "f"
            pos2, frag2, ref_len2 = fragment()
            right = frag2[-length:].translate(comp)[::-1]
            pos_right = pos2 + ref_len2 - length + 1
        if rng.random() < 0.5:
            r1, r2, strand, p1, p2 = left, right, "+", pos_left, pos_right
        else:
            r1, r2, strand, p1, p2 = right, left, "-", pos_right, pos_left
        if kind == "x" and strand == "-":
            # Der fremde Read soll immer Read 2 sein.
            r1, r2, p1, p2 = r2, r1, p2, p1
        title = f"{name}:{p1}:{p2}:{strand}:{kind}:{made}"
        q = "I" * length
        f1.write(f"@{title}\n{with_errors(r1)}\n+\n{q}\n")
        f2.write(f"@{title}\n{with_errors(r2)}\n+\n{q}\n")
