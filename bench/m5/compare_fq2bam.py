#!/usr/bin/env python3
"""Prüft `genome fq2bam` als Ganzes.

Aufruf: compare_fq2bam.py <genome> <index.gcfm> <r1.fq> <r2.fq> <arbeitsordner> <ergebnis.json>

 A. Gegen die einzelnen Schritte von Genome 4.0 (`index sam --partner-angaben`, `sort`, `markdup`):
    dieselben Zeilen in derselben Reihenfolge, mit kleinem und großem Speicher und mit 1 und 4
    Kernen.
 B. Gegen die Kette mit samtools (`fixmate -m`, `sort`, `markdup`) auf den Zeilen von
    `genome index sam`: alles gleich außer TLEN und der Reihenfolge der Zusatzangaben.
 C. Index und Zahlen: `samtools quickcheck`, `idxstats` mit dem Index von fq2bam gleich wie mit
    dem von samtools, `flagstat` und die Statistik der Duplikate im Bericht gleich wie bei samtools.
"""
import json, os, shutil, subprocess, sys, time
import pysam

genome, index, r1, r2, work, out_path = sys.argv[1:7]
os.makedirs(work, exist_ok=True)
p = lambda name: os.path.join(work, name)
RG = "@RG\\tID:l1\\tSM:p"


def run(*args):
    return subprocess.run([genome, *args], check=True, capture_output=True, text=True).stdout


def lines(path):
    return pysam.view(path).splitlines()


def canon(line):
    f = line.split("\t")
    f[8] = "."
    return "\t".join(f[:11] + sorted(f[11:]))


runs = {}
for name, extra in {
    "kerne4_speicher2048": ["--kerne", "4"],
    "kerne1_speicher2048": ["--kerne", "1"],
    "kerne4_speicher8": ["--kerne", "4", "--speicher-mb", "8"],
}.items():
    t = time.time()
    rep = json.loads(run("fq2bam", index, "--fastq", r1, "--fastq2", r2, "--aus", p(name + ".bam"), "--lesegruppe", RG, *extra))["result"]
    runs[name] = (rep, round(time.time() - t, 2))
main_name = "kerne4_speicher2048"
main, report = lines(p(main_name + ".bam")), runs[main_name][0]

# A: einzelne Schritte
run("index", "sam", index, "--fastq", r1, "--fastq2", r2, "--bam", p("a_roh.bam"), "--lesegruppe", RG, "--partner-angaben")
run("sort", p("a_roh.bam"), p("a_sortiert.bam"))
run("markdup", p("a_sortiert.bam"), p("a_markiert.bam"))
steps = lines(p("a_markiert.bam"))
a = {
    "zeilen": len(main),
    "gleich_wie_einzelne_schritte": main == steps,
    "laeufe": {n: {"gleich_wie_hauptlauf": lines(p(n + ".bam")) == main, "teile": r["sortieren"]["teile"], "sekunden": s} for n, (r, s) in runs.items()},
}

# B: samtools-Kette
run("index", "sam", index, "--fastq", r1, "--fastq2", r2, "--bam", p("s_roh.bam"), "--lesegruppe", RG)
pysam.fixmate("-m", p("s_roh.bam"), p("s_fix.bam"), catch_stdout=False)
pysam.sort("-o", p("s_sortiert.bam"), p("s_fix.bam"))
pysam.markdup("-f", p("s_stats.json"), "--json", p("s_sortiert.bam"), p("s_markiert.bam"), catch_stdout=False)
chain = lines(p("s_markiert.bam"))
same = sum(canon(x) == canon(y) for x, y in zip(main, chain))
b = {
    "zeilen_samtools": len(chain),
    "zeilen_gleich_ohne_tlen": same,
    "alle_gleich_ohne_tlen": same == len(main) == len(chain),
    "zeilen_mit_anderem_tlen": sum(x.split("\t")[8] != y.split("\t")[8] for x, y in zip(main, chain)),
    "markiert_fq2bam": sum(int(l.split("\t")[1]) & 1024 != 0 for l in main),
    "markiert_samtools": sum(int(l.split("\t")[1]) & 1024 != 0 for l in chain),
}

# C: Index und Zahlen
bam = p(main_name + ".bam")
quick = subprocess.run([sys.executable, "-c", "import pysam,sys; pysam.quickcheck(sys.argv[1])", bam]).returncode == 0
ours_idx = pysam.idxstats(bam)
shutil.copy(bam, p("kopie.bam"))
pysam.index(p("kopie.bam"))
fs = json.loads(pysam.flagstat("-O", "json", bam))["QC-passed reads"]
md = json.load(open(p("s_stats.json")))
names = {"total": "gesamt", "mapped": "mit_ort", "duplicates": "duplikate", "properly paired": "richtiges_paar",
         "singletons": "nur_selbst_mit_ort", "paired in sequencing": "aus_paaren", "with itself and mate mapped": "beide_mit_ort"}
c = {
    "quickcheck": quick,
    "idxstats_gleich_wie_samtools_index": ours_idx == pysam.idxstats(p("kopie.bam")),
    "flagstat_im_bericht_gleich": all(fs[k] == report["flagstat"]["bestanden"][v] for k, v in names.items()),
    "duplikat_statistik_gleich": all(md[k] == report["duplikate"][v] for k, v in {
        "EXAMINED": "bewertet", "PAIRED": "mit_partner", "SINGLE": "ohne_partner", "DUPLICATE PAIR": "duplikate_mit_partner",
        "DUPLICATE SINGLE": "duplikate_ohne_partner", "ESTIMATED_LIBRARY_SIZE": "geschaetzte_bibliotheksgroesse"}.items()),
    "kopf": pysam.view("-H", "--no-PG", bam).splitlines()[:1] + [l[:60] for l in pysam.view("-H", "--no-PG", bam).splitlines()[-2:]],
}
result = {"samtools": pysam.__samtools_version__, "A_einzelne_schritte": a, "B_samtools_kette": b, "C_index_und_zahlen": c,
          "bericht_sekunden_schritte": report["sekunden_schritte"], "bericht_reads": report["reads"]}
ok = a["gleich_wie_einzelne_schritte"] and all(v["gleich_wie_hauptlauf"] for v in a["laeufe"].values()) and b["alle_gleich_ohne_tlen"] \
    and quick and c["idxstats_gleich_wie_samtools_index"] and c["flagstat_im_bericht_gleich"] and c["duplikat_statistik_gleich"]
result["alles_gleich"] = ok
json.dump(result, open(out_path, "w"), indent=1, ensure_ascii=False)
print(json.dumps(result, indent=1, ensure_ascii=False))
sys.exit(0 if ok else 1)
