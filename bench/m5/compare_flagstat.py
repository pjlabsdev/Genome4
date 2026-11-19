#!/usr/bin/env python3
"""Vergleicht `genome flagstat` mit `samtools flagstat -O json` (über pysam).

Aufruf: compare_flagstat.py <genome> <arbeitsordner> <ergebnis.json> <a.bam> [b.bam ...]

Zusätzlich zu den übergebenen Dateien entsteht eine künstliche BAM mit zufälligen Flags (alle
zwölf Bits), Contigs, Partner-Contigs und MAPQ, damit jede Zeile der Zählung und beide Gruppen
(Qualitätsprüfung bestanden oder nicht) Werte bekommen.
"""
import json, os, random, subprocess, sys
import pysam

genome, work, out_path, *bams = sys.argv[1:]
os.makedirs(work, exist_ok=True)

synth = os.path.join(work, "zufallsflags.bam")
rng = random.Random(3)
hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": f"c{i}", "LN": 9000} for i in range(4)]})
with pysam.AlignmentFile(synth, "wb", header=hdr) as o:
    for i in range(400000):
        r = pysam.AlignedSegment(hdr)
        r.query_name = f"r{i}"
        r.query_sequence = "ACGT"
        r.flag = rng.randrange(4096)
        r.reference_id = rng.randrange(-1, 4)
        r.reference_start = rng.randrange(8000)
        r.next_reference_id = rng.randrange(-1, 4)
        r.next_reference_start = rng.randrange(8000)
        r.mapping_quality = rng.choice([0, 3, 4, 5, 6, 30, 60])
        if not r.flag & 4:
            r.cigarstring = "4M"
        o.write(r)
bams.append(synth)

NAMES = {
    "total": "gesamt", "primary": "hauptzeilen", "secondary": "sekundaer", "supplementary": "ergaenzend",
    "duplicates": "duplikate", "primary duplicates": "duplikate_hauptzeilen", "mapped": "mit_ort",
    "primary mapped": "hauptzeilen_mit_ort", "paired in sequencing": "aus_paaren", "read1": "read1",
    "read2": "read2", "properly paired": "richtiges_paar", "with itself and mate mapped": "beide_mit_ort",
    "singletons": "nur_selbst_mit_ort", "with mate mapped to a different chr": "partner_auf_anderem_contig",
    "with mate mapped to a different chr (mapQ >= 5)": "partner_auf_anderem_contig_mapq5",
}
result = {"samtools": pysam.__samtools_version__, "dateien": []}
ok = True
for bam in bams:
    ours = json.loads(subprocess.run([genome, "flagstat", bam], check=True, capture_output=True, text=True).stdout)["result"]
    theirs = json.loads(pysam.flagstat("-O", "json", bam))
    case = {"datei": os.path.basename(bam), "abweichungen": [], "zahlen": {}}
    for group_s, group_g in (("QC-passed reads", "bestanden"), ("QC-failed reads", "nicht_bestanden")):
        for name_s, name_g in NAMES.items():
            a, b = theirs[group_s][name_s], ours[group_g][name_g]
            case["zahlen"][f"{group_g}.{name_g}"] = b
            if a != b:
                case["abweichungen"].append({"zahl": f"{group_g}.{name_g}", "samtools": a, "genome": b})
    case["verglichene_zahlen"] = 2 * len(NAMES)
    ok &= not case["abweichungen"]
    result["dateien"].append(case)
result["alles_gleich"] = ok
json.dump(result, open(out_path, "w"), indent=1, ensure_ascii=False)
for c in result["dateien"]:
    print(c["datei"], c["verglichene_zahlen"], "Zahlen,", len(c["abweichungen"]), "Abweichungen", c["abweichungen"][:3], c["zahlen"]["bestanden.gesamt"], c["zahlen"]["nicht_bestanden.gesamt"], c["zahlen"]["bestanden.duplikate"])
sys.exit(0 if ok else 1)
