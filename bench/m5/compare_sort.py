#!/usr/bin/env python3
"""Vergleicht `genome sort` mit `samtools sort` (über pysam) Satz für Satz.

Aufruf: compare_sort.py <genome> <eingabe.bam> <arbeitsordner> <ergebnis.json>

Geprüft wird für mehrere Speichergrößen (mit und ohne Teile auf der Platte) und für eine
gemischte Reihenfolge der Eingabe:
  - die Folge der Sätze ist Byte für Byte die von samtools (alle Felder, über `samtools view`),
  - kein Satz fehlt, keiner ist doppelt (Vergleich als Menge mit der Eingabe),
  - der Kopf nennt SO:coordinate, die übrigen Kopfzeilen sind unverändert.
"""
import hashlib, json, os, random, subprocess, sys, time
import pysam

genome, src, work, out_path = sys.argv[1:5]
os.makedirs(work, exist_ok=True)


def records(path):
    return pysam.view(path).splitlines()


def header(path):
    return pysam.view("-H", "--no-PG", path).splitlines()


def digest(lines):
    h = hashlib.sha256()
    for l in lines:
        h.update(l.encode() + b"\n")
    return h.hexdigest()


def shuffled(path, dst, seed):
    """Dieselben Sätze in zufälliger Reihenfolge, damit die Eingabe keine Vorordnung hat."""
    with pysam.AlignmentFile(path, check_sq=False) as f:
        recs = list(f)
        random.Random(seed).shuffle(recs)
        with pysam.AlignmentFile(dst, "wb", header=f.header) as o:
            for r in recs:
                o.write(r)


inputs = {"wie_ausgerichtet": src}
mixed = os.path.join(work, "gemischt.bam")
shuffled(src, mixed, 20261119)
inputs["gemischt"] = mixed


def synthetic(dst, n, seed):
    """Mehrere Contigs, viele gleiche Orte, beide Stränge, Sätze ohne Contig mit und ohne Position."""
    rng = random.Random(seed)
    names = [f"c{i}" for i in range(7)]
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": s, "LN": 5000} for s in names]}
    )
    with pysam.AlignmentFile(dst, "wb", header=hdr) as o:
        for i in range(n):
            r = pysam.AlignedSegment(hdr)
            r.query_name = f"r{i}"
            r.query_sequence = "ACGTACGTAC"
            r.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
            kind = rng.random()
            if kind < 0.05:
                r.reference_id, r.reference_start, r.flag = -1, -1, 4
            elif kind < 0.10:
                # ohne eigenen Ort, aber beim Partner einsortiert
                r.reference_id, r.reference_start = rng.randrange(7), rng.randrange(200)
                r.flag = 4 | 1 | (16 if rng.random() < 0.5 else 0)
            else:
                r.reference_id, r.reference_start = rng.randrange(7), rng.randrange(200)
                r.flag = 16 if rng.random() < 0.5 else 0
                r.cigarstring = "10M"
                r.mapping_quality = rng.randrange(61)
            o.write(r)


synth = os.path.join(work, "kuenstlich.bam")
synthetic(synth, 300000, 5)
inputs["kuenstlich_7_contigs"] = synth

result = {"samtools": pysam.__samtools_version__, "pysam": pysam.__version__, "faelle": []}
ok = True
for name, path in inputs.items():
    truth_path = os.path.join(work, f"samtools_{name}.bam")
    t = time.time()
    pysam.sort("-o", truth_path, "-@", "1", path)
    t_samtools = time.time() - t
    truth = records(truth_path)
    source = records(path)
    in_header = header(path)
    for mb in (2048, 16, 2):
        ours_path = os.path.join(work, f"genome_{name}_{mb}.bam")
        t = time.time()
        env = json.loads(
            subprocess.run(
                [genome, "sort", path, ours_path, "--speicher-mb", str(mb), "--kerne", "1"],
                check=True, capture_output=True, text=True,
            ).stdout
        )["result"]
        t_ours = time.time() - t
        ours = records(ours_path)
        h = header(ours_path)
        case = {
            "eingabe": name,
            "speicher_mb": mb,
            "saetze": len(ours),
            "teile": env["teile"],
            "gleich_wie_samtools": ours == truth,
            "erste_abweichung": next((i for i, (a, b) in enumerate(zip(ours, truth)) if a != b), None),
            "gleiche_menge_wie_eingabe": sorted(ours) == sorted(source),
            "kopf_so_coordinate": h[0].startswith("@HD") and "SO:coordinate" in h[0].split("\t"),
            "kopf_rest_unveraendert": [l for l in h if not l.startswith("@HD") and "ID:genome-sort" not in l]
            == [l for l in in_header if not l.startswith("@HD")],
            "quickcheck": subprocess.run(
                [sys.executable, "-c", "import pysam,sys; pysam.quickcheck(sys.argv[1])", ours_path]
            ).returncode == 0,
            "sha256_saetze": digest(ours)[:12],
            "sekunden_genome": round(t_ours, 3),
            "sekunden_samtools": round(t_samtools, 3),
        }
        ok &= all(case[k] for k in ("gleich_wie_samtools", "gleiche_menge_wie_eingabe", "kopf_so_coordinate", "kopf_rest_unveraendert", "quickcheck"))
        result["faelle"].append(case)
result["alles_gleich"] = ok
json.dump(result, open(out_path, "w"), indent=1, ensure_ascii=False)
print(json.dumps(result, indent=1, ensure_ascii=False))
sys.exit(0 if ok else 1)
