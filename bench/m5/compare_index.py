#!/usr/bin/env python3
"""Vergleicht den Index von `genome bam-index` mit dem von `samtools index` (über pysam).

Aufruf: compare_index.py <genome> <arbeitsordner> <ergebnis.json> <sortiert.bam> [weitere.bam ...]

Je BAM:
 1. Inhalt beider `.bai`-Dateien zerlegt und verglichen: je Contig die Fächer mit ihren
    Dateistücken, der lineare Index, das Sonderfach mit den Zählern, die Zahl ohne Contig.
    (Die Dateien selbst müssen nicht Byte für Byte gleich sein: htslib schreibt die Fächer in der
    Reihenfolge seiner Hashtabelle, Genome 4.0 aufsteigend.)
 2. Abfragen: 2.000 zufällige Bereiche, einmal mit dem einen, einmal mit dem anderen Index über
    htslib gelesen. Die gelieferten Datensätze müssen gleich sein. Dazu `idxstats`.
"""
import json, os, random, shutil, struct, subprocess, sys
import pysam

genome, work, out_path, *bams = sys.argv[1:]
os.makedirs(work, exist_ok=True)


def parse(path):
    d = open(path, "rb").read()
    assert d[:4] == b"BAI\1"
    at = 4
    def take(fmt):
        nonlocal at
        v = struct.unpack_from("<" + fmt, d, at)
        at += struct.calcsize("<" + fmt)
        return v
    refs = []
    for _ in range(take("i")[0]):
        bins = {}
        for _ in range(take("i")[0]):
            b, n = take("Ii")
            bins[b] = [take("QQ") for _ in range(n)]
        (n,) = take("i")
        linear = list(take(f"{n}Q")) if n else []
        refs.append((bins, linear))
    unplaced = take("Q")[0] if at < len(d) else None
    return refs, unplaced


result = {"samtools": pysam.__samtools_version__, "dateien": []}
ok = True
for bam in bams:
    name = os.path.basename(bam)
    a, b = os.path.join(work, "g_" + name), os.path.join(work, "s_" + name)
    shutil.copy(bam, a)
    shutil.copy(bam, b)
    env = json.loads(subprocess.run([genome, "bam-index", a], check=True, capture_output=True, text=True).stdout)["result"]
    pysam.index(b)
    (ra, ua), (rb, ub) = parse(a + ".bai"), parse(b + ".bai")
    case = {
        "datei": name,
        "saetze": env["saetze"],
        "contigs": len(ra),
        "faecher_genome": sum(len(x[0]) for x in ra),
        "faecher_samtools": sum(len(x[0]) for x in rb),
        "faecher_gleich": all(x[0] == y[0] for x, y in zip(ra, rb)) and len(ra) == len(rb),
        "linearer_index_gleich": all(x[1] == y[1] for x, y in zip(ra, rb)),
        "ohne_contig_gleich": ua == ub,
        "dateien_bytegleich": open(a + ".bai", "rb").read() == open(b + ".bai", "rb").read(),
        "bytes_index": env["bytes_index"],
    }
    # Abfragen
    rng = random.Random(11)
    fa, fb = pysam.AlignmentFile(a), pysam.AlignmentFile(b)
    refs = [(r, l) for r, l in zip(fa.references, fa.lengths)]
    different = 0
    returned = 0
    for i in range(2000):
        ref, length = rng.choice(refs)
        span = rng.choice([1, 100, 5000, 200000, 5000000])
        start = rng.randrange(max(1, length))
        end = min(length, start + span)
        if i % 4 == 0 and fa.mapped:
            # ein Bereich um einen vorhandenen Datensatz, damit nicht alles leer ist
            pass
        xa = [r.to_string() for r in fa.fetch(ref, start, end)]
        xb = [r.to_string() for r in fb.fetch(ref, start, end)]
        returned += len(xa)
        different += xa != xb
    case["abfragen"] = 2000
    case["abfragen_verschieden"] = different
    case["datensaetze_geliefert"] = returned
    case["ohne_ort_gleich"] = [r.to_string() for r in fa.fetch(until_eof=True)] == [r.to_string() for r in fb.fetch(until_eof=True)]
    case["idxstats_gleich"] = pysam.idxstats(a) == pysam.idxstats(b)
    ok &= case["faecher_gleich"] and case["linearer_index_gleich"] and case["ohne_contig_gleich"] and different == 0 and case["idxstats_gleich"]
    result["dateien"].append(case)
result["alles_gleich"] = ok
json.dump(result, open(out_path, "w"), indent=1, ensure_ascii=False)
print(json.dumps(result, indent=1, ensure_ascii=False))
sys.exit(0 if ok else 1)
