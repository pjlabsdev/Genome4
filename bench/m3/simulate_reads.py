#!/usr/bin/env python3
"""Einfache Reads aus einer FASTA ziehen: feste Länge, beide Stränge, Substitutionsfehler.

Kein realistisches Fehlermodell (keine Qualitätswerte); der wahre Ort steht im Namen:
contig:pos(1-basiert):strang:fehler. Mit einer Indelrate (6. Angabe, je Base) kommen
Insertionen und Deletionen von 1 bis 10 Basen dazu; der Name bekommt dann zwei weitere Felder,
die Zahl der Indels und eine laufende Nummer, damit kein Name doppelt vorkommt. Ohne Indelrate
bleibt die Ausgabe genau wie bisher (auch mit den wenigen doppelten Namen).
Aufruf: simulate_reads.py ref.fa anzahl laenge fehlerrate saat [indelrate] > reads.fq
"""
import random
import sys

fa, n, length, rate, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
indel_rate = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
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
    span = length + 64 if indel_rate > 0 else length
    pos = rng.randrange(len(seq) - span)
    piece = seq[pos:pos + span]
    if "N" in piece:
        continue
    indels = 0
    if indel_rate > 0:
        # Indels auf dem Vorwärtsstrang setzen, nicht in den ersten und letzten 5 Basen.
        made_bases, i = [], 0
        while len(made_bases) < length:
            if 5 <= len(made_bases) < length - 5 and rng.random() < indel_rate:
                size = rng.randint(1, 10)
                indels += 1
                if rng.random() < 0.5:
                    i += size  # Deletion: Referenz überspringen
                else:
                    made_bases.extend(rng.choice("ACGT") for _ in range(size))  # Insertion
                    continue
            made_bases.append(piece[i])
            i += 1
        piece = "".join(made_bases[:length])
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
    extra = f":{indels}:{made}" if indel_rate > 0 else ""
    out.write(f"@{name}:{pos + 1}:{strand}:{errors}{extra}\n{''.join(bases)}\n+\n{'I' * length}\n")
    made += 1
