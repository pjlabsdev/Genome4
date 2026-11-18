#!/usr/bin/env python3
"""M2-Prüfung: exakte Trefferzahlen und Positionen des Genome-FM-Index gegen bwa.

bwa kommt hier über das Python-Paket pybwa (enthält den Original-C-Code von bwa).
Aufruf: compare_bwa.py <genome-Programm> <index.gcfm> <referenz.fa> <anzahl> <ausgabe.json>
Die bwa-Indexdateien liegen neben der FASTA (werden gebaut, wenn sie fehlen).
"""
import json, os, random, subprocess, sys, tempfile, time

import pybwa
import pysam

genome, index, fasta, n, out_path = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]

if not os.path.exists(fasta + ".bwt"):
    pybwa.BwaIndex.index(fasta=fasta)

fa = pysam.FastaFile(fasta)
contig = fa.references[0]
seq = fa.fetch(contig).upper()
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
        piece = piece.translate(comp)[::-1]          # Gegenstrang
    elif kind == 3:
        at = rng.randrange(length)                    # eine Base geändert: meist kein exakter Treffer
        piece = piece[:at] + rng.choice([b for b in "ACGT" if b != piece[at]]) + piece[at + 1:]
    queries.append(piece)

with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
    f.write("\n".join(queries) + "\n")
    qfile = f.name

t0 = time.time()
ours = subprocess.run([genome, "index", "count", index, qfile], check=True, capture_output=True, text=True).stdout
t_ours = time.time() - t0
our_counts = [int(line.split("\t")[1]) for line in ours.splitlines()]
assert len(our_counts) == n

opt = pybwa.BwaAlnOptions()
opt.max_mismatches = 0          # nur exakte Treffer
opt.max_gap_opens = 0
opt.seed_length = 1000          # ohne Seed-Abkürzung, ganze Lesung exakt
opt.threads = os.cpu_count() or 1
t0 = time.time()
with pybwa.BwaAln(prefix=fasta) as aln:
    recs = aln.align(queries=queries, opt=opt)
t_bwa = time.time() - t0

count_diff, examples = 0, []
with_hits = 0
pos_checked = pos_diff = 0
for i, (q, rec, mine) in enumerate(zip(queries, recs, our_counts)):
    bwa_count = 0 if rec.is_unmapped else rec.get_tag("X0")
    with_hits += bwa_count > 0
    if bwa_count != mine:
        count_diff += 1
        if len(examples) < 10:
            examples.append({"abfrage": i, "laenge": len(q), "bwa": bwa_count, "genome": mine})
    # Positionen: jede 200. Abfrage mit Treffer einzeln orten und mit bwa vergleichen
    if bwa_count > 0 and i % 200 == 0 and bwa_count <= 1000:
        res = json.loads(subprocess.run([genome, "index", "locate", index, q, "--limit", "100000"],
                                        check=True, capture_output=True, text=True).stdout)
        hits = {(h["pos"], h["strang"]) for h in res["result"]["treffer"]}
        pos_checked += 1
        if (rec.reference_start + 1, "-" if rec.is_reverse else "+") not in hits:
            pos_diff += 1

result = {
    "referenz": os.path.basename(fasta), "contig": contig, "abfragen": n,
    "abfragen_mit_treffern_bwa": with_hits,
    "abweichende_trefferzahlen": count_diff, "beispiele": examples,
    "positionen_geprueft": pos_checked, "positionen_abweichend": pos_diff,
    "bwa": "pybwa " + getattr(pybwa, "__version__", "?"),
    "sekunden_genome_count": round(t_ours, 2), "sekunden_bwa_aln": round(t_bwa, 2),
}
json.dump(result, open(out_path, "w"), indent=2, ensure_ascii=False)
print(json.dumps(result, indent=2, ensure_ascii=False))
os.unlink(qfile)
