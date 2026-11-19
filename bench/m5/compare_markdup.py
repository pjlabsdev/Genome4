#!/usr/bin/env python3
"""Vergleicht `genome markdup` mit `samtools markdup` (über pysam) und mit der Wahrheit.

Aufruf: compare_markdup.py <genome> <index.gcfm> <r1.fq> <r2.fq> <wahrheit.tsv> <arbeitsordner> <ergebnis.json>

Drei Prüfungen:
 A. Gleiche Eingabe: `genome index sam --partner-angaben`, `genome sort`, dann `genome markdup`
    und `samtools markdup` auf derselben sortierten BAM. Verglichen wird jede Zeile (alle Felder)
    und die Statistik.
 B. Ganze Kette: `genome index sam` (ohne Partner-Angaben), `samtools fixmate -m`,
    `samtools sort`, `samtools markdup` gegen die Kette von Genome 4.0. Verglichen wird je Zeile
    alles außer TLEN (das fixmate anders rechnet als bwa) und der Reihenfolge der Zusatzangaben.
 C. Wahrheit: Je Gruppe von Kopien (aus `make_duplicates.py`) soll genau ein Paar unmarkiert
    bleiben, und zwar eines mit der höchsten Qualitätssumme.
"""
import collections, json, os, subprocess, sys, time
import pysam

genome, index, r1, r2, truth_path, work, out_path = sys.argv[1:8]
os.makedirs(work, exist_ok=True)
p = lambda name: os.path.join(work, name)
RG = "@RG\\tID:l1\\tSM:p"


def run(*args):
    return subprocess.run([genome, *args], check=True, capture_output=True, text=True).stdout


def lines(path):
    return pysam.view(path).splitlines()


def canon(line, drop_tlen=False):
    f = line.split("\t")
    if drop_tlen:
        f[8] = "."
    return "\t".join(f[:11] + sorted(f[11:]))


# Kette Genome 4.0
run("index", "sam", index, "--fastq", r1, "--fastq2", r2, "--bam", p("g_roh.bam"), "--lesegruppe", RG, "--partner-angaben")
run("sort", p("g_roh.bam"), p("g_sortiert.bam"))
t = time.time()
ours_stats = json.loads(run("markdup", p("g_sortiert.bam"), p("g_markiert.bam"), "--kerne", "1"))["result"]
t_ours = time.time() - t

# A: samtools markdup auf derselben Eingabe
t = time.time()
pysam.markdup("-f", p("s_stats.json"), "--json", "-@", "0", p("g_sortiert.bam"), p("s_markiert_a.bam"), catch_stdout=False)
t_samtools = time.time() - t
ours = lines(p("g_markiert.bam"))
theirs = lines(p("s_markiert_a.bam"))
s = json.load(open(p("s_stats.json")))
stat_pairs = {
    "READ": "gelesen", "WRITTEN": "geschrieben", "EXCLUDED": "nicht_bewertet", "EXAMINED": "bewertet",
    "PAIRED": "mit_partner", "SINGLE": "ohne_partner", "DUPLICATE PAIR": "duplikate_mit_partner",
    "DUPLICATE SINGLE": "duplikate_ohne_partner", "DUPLICATE TOTAL": "duplikate",
    "ESTIMATED_LIBRARY_SIZE": "geschaetzte_bibliotheksgroesse",
}
a = {
    "zeilen": len(ours),
    "zeilen_gleich": sum(x == y for x, y in zip(ours, theirs)),
    "alle_zeilen_gleich": ours == theirs,
    "markiert_genome": sum(int(l.split("\t")[1]) & 1024 != 0 for l in ours),
    "markiert_samtools": sum(int(l.split("\t")[1]) & 1024 != 0 for l in theirs),
    "statistik": {k: {"samtools": s[k], "genome": ours_stats[v]} for k, v in stat_pairs.items()},
    "sekunden_genome": round(t_ours, 3),
    "sekunden_samtools": round(t_samtools, 3),
}
a["statistik_gleich"] = all(v["samtools"] == (v["genome"] or 0) for v in a["statistik"].values())

# B: ganze Kette mit samtools
run("index", "sam", index, "--fastq", r1, "--fastq2", r2, "--bam", p("s_roh.bam"), "--lesegruppe", RG)
pysam.fixmate("-m", p("s_roh.bam"), p("s_fix.bam"), catch_stdout=False)
pysam.sort("-o", p("s_sortiert.bam"), p("s_fix.bam"))
pysam.markdup(p("s_sortiert.bam"), p("s_markiert_b.bam"), catch_stdout=False)
chain = lines(p("s_markiert_b.bam"))
same = sum(canon(x, True) == canon(y, True) for x, y in zip(ours, chain))
b = {
    "zeilen": len(chain),
    "zeilen_gleich_ohne_tlen": same,
    "alle_gleich_ohne_tlen": same == len(ours) == len(chain),
    "zeilen_mit_anderem_tlen": sum(x.split("\t")[8] != y.split("\t")[8] for x, y in zip(ours, chain)),
    "markiert_samtools_kette": sum(int(l.split("\t")[1]) & 1024 != 0 for l in chain),
}

# C: Wahrheit
group_of = {}
with open(truth_path) as f:
    next(f)
    for line in f:
        name, group, _ = line.rstrip("\n").split("\t")
        group_of[name] = group
score = collections.defaultdict(int)  # Summe der Qualitäten ab 15 je Paar (Hauptzeilen mit Ort)
dup = collections.defaultdict(list)
placed = collections.defaultdict(int)
for l in ours:
    f = l.split("\t")
    flag = int(f[1])
    if flag & (256 | 2048):
        continue
    if not flag & 4:
        placed[f[0]] += 1
        dup[f[0]].append(bool(flag & 1024))
    score[f[0]] += sum(q - 33 for q in f[10].encode() if q - 33 >= 15)
groups = collections.defaultdict(list)
for name, group in group_of.items():
    if placed[name] == 2:  # nur Paare, bei denen beide Reads einen Ort haben
        groups[group].append(name)
c = collections.Counter()
for group, names in groups.items():
    if len(names) < 2:
        c["gruppen_mit_einem_paar"] += 1
        c["davon_faelschlich_markiert"] += any(dup[names[0]])
        continue
    c["gruppen_mit_kopien"] += 1
    unmarked = [n for n in names if not any(dup[n])]
    half = [n for n in names if any(dup[n]) and not all(dup[n])]
    if len(unmarked) == 1 and not half:
        c["genau_ein_paar_bleibt"] += 1
        c["das_beste_bleibt"] += score[unmarked[0]] == max(score[n] for n in names)
    elif len(unmarked) > 1:
        c["mehrere_bleiben"] += 1
    else:
        c["anderes"] += 1
result = {
    "samtools": pysam.__samtools_version__,
    "A_gleiche_eingabe": a,
    "B_ganze_kette": b,
    "C_wahrheit": dict(c),
}
json.dump(result, open(out_path, "w"), indent=1, ensure_ascii=False)
print(json.dumps(result, indent=1, ensure_ascii=False))
