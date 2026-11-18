#!/usr/bin/env python3
"""Erklärt Abweichungen der rohen Trefferzahl zwischen bwa und Genome.

Rohe Zahlen (Breite des Suffix-Array-Bereichs) zählen auch Fundstellen, die in eine N-Strecke
hineinreichen: dort stehen Zufallsbasen, bei bwa andere als bei uns. Die Wahrheit ist die direkte
Suche in der FASTA (ohne N). Geprüft wird: gefilterte Genome-Treffer == Wahrheit.
Aufruf: explain_diffs.py <genome-Programm> <index.gcfm> <referenz.fa> <anzahl> <ausgabe.json>
(gleiche Abfragen wie compare_bwa.py, gleicher Startwert)
"""
import json, os, random, re, subprocess, sys, tempfile
import pybwa, pysam

genome, index, fasta, n, out_path = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
fa = pysam.FastaFile(fasta)
seq = fa.fetch(fa.references[0]).upper()
comp = str.maketrans("ACGT", "TGCA")
rng = random.Random(20261118)
queries = []
while len(queries) < n:
    length = rng.randint(20, 100)
    start = rng.randrange(0, len(seq) - length)
    piece = seq[start:start + length]
    if "N" in piece:
        continue
    kind = len(queries) % 4
    if kind == 1:
        piece = piece.translate(comp)[::-1]
    elif kind == 3:
        at = rng.randrange(length)
        piece = piece[:at] + rng.choice([b for b in "ACGT" if b != piece[at]]) + piece[at + 1:]
    queries.append(piece)

with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
    f.write("\n".join(queries) + "\n")
    qfile = f.name
ours = [int(l.split("\t")[1]) for l in subprocess.run([genome, "index", "count", index, qfile], check=True, capture_output=True, text=True).stdout.splitlines()]
opt = pybwa.BwaAlnOptions()
opt.max_mismatches = 0; opt.max_gap_opens = 0; opt.seed_length = 1000; opt.threads = os.cpu_count() or 1
with pybwa.BwaAln(prefix=fasta) as aln:
    recs = aln.align(queries=queries, opt=opt)

def truth(q):
    total = 0
    for pat in (q, q.translate(comp)[::-1]):
        total += sum(1 for _ in re.finditer("(?=" + pat + ")", seq))
    return total

rows = []
for i, (q, rec, mine) in enumerate(zip(queries, recs, ours)):
    bwa = 0 if rec.is_unmapped else rec.get_tag("X0")
    if bwa == mine:
        continue
    res = json.loads(subprocess.run([genome, "index", "locate", index, q, "--limit", "100000000"], check=True, capture_output=True, text=True).stdout)
    rows.append({"abfrage": i, "laenge": len(q), "bwa_roh": bwa, "genome_roh": mine,
                 "genome_gefiltert": len(res["result"]["treffer"]), "wahrheit_fasta": truth(q)})
ok = sum(r["genome_gefiltert"] == r["wahrheit_fasta"] for r in rows)
result = {"abweichungen_roh": len(rows), "davon_genome_gefiltert_gleich_wahrheit": ok,
          "bwa_roh_ueber_wahrheit": sum(r["bwa_roh"] > r["wahrheit_fasta"] for r in rows),
          "genome_roh_ueber_wahrheit": sum(r["genome_roh"] > r["wahrheit_fasta"] for r in rows),
          "bwa_roh_unter_wahrheit": sum(r["bwa_roh"] < r["wahrheit_fasta"] for r in rows),
          "genome_roh_unter_wahrheit": sum(r["genome_roh"] < r["wahrheit_fasta"] for r in rows),
          "zeilen": rows}
json.dump(result, open(out_path, "w"), indent=1)
print(json.dumps({k: v for k, v in result.items() if k != "zeilen"}, indent=1)); print(rows[:6])
os.unlink(qfile)
