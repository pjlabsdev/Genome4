#!/usr/bin/env python3
"""Erzeugt Prüffälle für das lokale Alignment mit den Ergebnissen von bwa.

Ruft `ksw_align2` aus der bwa-Bibliothek, die pybwa mitbringt, direkt auf (wie bwa bei der
Partnersuche: zweitbester Score ab 19, mit Anfang). Jede Zeile: Read, Referenz (Ziffern 0 bis 4,
4 = N), Byte-Modus (1/0), dann score, te, qe, score2, te2, tb, qb.
Aufruf: make_local_cases.py anzahl saat > faelle.tsv
"""
import ctypes
import glob
import os
import random
import sys

import pybwa

n, seed = int(sys.argv[1]), int(sys.argv[2])
lib = ctypes.CDLL(glob.glob(os.path.join(os.path.dirname(pybwa.__file__), "libbwamem*.so"))[0])


class Result(ctypes.Structure):
    _fields_ = [(k, ctypes.c_int) for k in ("score", "te", "qe", "score2", "te2", "tb", "qb")]


lib.ksw_align2.restype = Result
lib.ksw_align2.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int,
                           ctypes.c_char_p] + [ctypes.c_int] * 5 + [ctypes.c_void_p]
MAT = bytes((1 if i == j else -4) & 0xFF if i < 4 and j < 4 else 0xFF for i in range(5) for j in range(5))
XSUBO, XSTART, XBYTE = 0x40000, 0x80000, 0x10000
rng = random.Random(seed)


def mutate(q):
    q = list(q)
    for _ in range(rng.choice([0, 0, 1, 2, 3, 6, 12])):
        if len(q) < 30:
            break
        p = rng.randrange(len(q))
        kind = rng.random()
        if kind < 0.6:
            q[p] = rng.choice([0, 1, 2, 3, 3, 4] if rng.random() < 0.1 else [0, 1, 2, 3])
        elif kind < 0.8:
            del q[p:p + rng.randint(1, 8)]
        else:
            q[p:p] = [rng.randrange(4) for _ in range(rng.randint(1, 8))]
    return q


for _ in range(n):
    tlen = rng.randint(40, 900)
    if rng.random() < 0.2:
        unit = [rng.randrange(4) for _ in range(rng.randint(1, 40))]  # Wiederholung
        t = (unit * (tlen // len(unit) + 1))[:tlen]
        t = [b if rng.random() > 0.03 else rng.randrange(4) for b in t]
    else:
        t = [rng.randrange(4) for _ in range(tlen)]
    if rng.random() < 0.25 and tlen > 300:  # zweite Kopie eines Abschnitts
        a, b, size = rng.randrange(tlen - 150), rng.randrange(tlen - 150), rng.randint(30, 150)
        t[b:b + size] = t[a:a + size]
    qlen = rng.choice([150, 150, 150, 100, 251, rng.randint(20, 300)])
    if rng.random() < 0.9:
        start = rng.randrange(-qlen // 2, tlen)
        q = [t[i] if 0 <= i < tlen else rng.randrange(4) for i in range(start, start + qlen)]
        q = mutate(q)
    else:
        q = [rng.randrange(4) for _ in range(qlen)]
    if rng.random() < 0.05:
        p = rng.randrange(len(t))
        t[p:p + rng.randint(1, 30)] = [4] * len(t[p:p + rng.randint(1, 30)])
    byte = len(q) < 250 and rng.random() < 0.8
    xtra = XSUBO | XSTART | (XBYTE if byte else 0) | 19
    r = lib.ksw_align2(len(q), bytes(q), len(t), bytes(t), 5, MAT, 6, 1, 6, 1, xtra, None)
    print("".join(map(str, q)), "".join(map(str, t)), int(byte), r.score, r.te, r.qe, r.score2, r.te2, r.tb, r.qb, sep="\t")
