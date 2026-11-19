#!/usr/bin/env python3
"""Erzeugt aus einem Satz von Read-Paaren einen Satz mit Duplikaten und bekannter Wahrheit.

Aufruf: make_duplicates.py <r1.fq> <r2.fq> <aus_1.fq> <aus_2.fq> <wahrheit.tsv> [seed]

- Jedes Paar bekommt zufällige Basenqualitäten (damit die Wahl des besten Paares etwas zu
  entscheiden hat).
- 15 % der Paare bekommen 1 bis 3 Kopien. Eine Kopie hat denselben Ursprung (dieselben Basen),
  aber eigene Qualitäten; ein Viertel der Kopien behält die Qualitäten (Gleichstand), ein Viertel
  bekommt bis zu drei Lesefehler, und bei einem Zehntel ist der zweite Read durch Zufallsbasen
  ersetzt (Partner ohne Ort).
- Namen: Original `<name>`, Kopie `<name>_k<n>`. Die Wahrheit nennt je Read die Gruppe (Name des
  Originals) und ob es eine Kopie ist.
Die Reihenfolge der Paare wird gemischt, damit Kopien nicht nebeneinander stehen.
"""
import random, sys

r1, r2, o1, o2, truth_path = sys.argv[1:6]
rng = random.Random(int(sys.argv[6]) if len(sys.argv) > 6 else 7)


def read(path):
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().strip()
            f.readline()
            q = f.readline().strip()
            yield h.strip()[1:].split()[0].removesuffix("/1").removesuffix("/2"), s, q


def quals(n):
    # überwiegend gut, mit einigen schwachen Stellen und einem schwächeren Ende
    out = []
    for i in range(n):
        r = rng.random()
        q = 40 if r < 0.8 else rng.choice([37, 30, 25, 20, 14, 11, 2])
        if i > n - 20 and rng.random() < 0.3:
            q = rng.choice([14, 11, 2, 20])
        out.append(chr(33 + q))
    return "".join(out)


def errors(seq, k):
    s = list(seq)
    for _ in range(k):
        i = rng.randrange(len(s))
        s[i] = rng.choice([c for c in "ACGT" if c != s[i]])
    return "".join(s)


pairs = []
for (n1, s1, _), (n2, s2, _) in zip(read(r1), read(r2)):
    assert n1 == n2
    q1, q2 = quals(len(s1)), quals(len(s2))
    pairs.append((n1, n1, 0, s1, q1, s2, q2))
    if rng.random() < 0.15:
        for k in range(1, rng.randint(1, 3) + 1):
            c1, c2, d1, d2 = s1, s2, quals(len(s1)), quals(len(s2))
            kind = rng.random()
            if kind < 0.25:
                d1, d2 = q1, q2
            elif kind < 0.5:
                c1, c2 = errors(s1, rng.randint(0, 3)), errors(s2, rng.randint(0, 3))
            elif kind < 0.6:
                c2 = "".join(rng.choice("ACGT") for _ in s2)
            pairs.append((f"{n1}_k{k}", n1, 1, c1, d1, c2, d2))
rng.shuffle(pairs)
with open(o1, "w") as f1, open(o2, "w") as f2, open(truth_path, "w") as t:
    t.write("name\tgruppe\tkopie\n")
    for name, group, copy, s1, q1, s2, q2 in pairs:
        f1.write(f"@{name}/1\n{s1}\n+\n{q1}\n")
        f2.write(f"@{name}/2\n{s2}\n+\n{q2}\n")
        t.write(f"{name}\t{group}\t{copy}\n")
print(len(pairs), "Paare,", sum(p[2] for p in pairs), "Kopien")
