# Vorläufig: pybwa (eingebettetes bwa mem) -> BAM mit Flags 64/128, für den Probelauf.
import sys, pybwa, pysam
prefix, fq1, fq2, out = sys.argv[1:5]
q=[]
for a,b in zip(pysam.FastxFile(fq1), pysam.FastxFile(fq2)): q += [a,b]
opt = pybwa.BwaMemOptions(interleaved_paired_end=True, chunk_size=10_000_000)
with pybwa.BwaMem(prefix=prefix) as al:
    res = al.align(q, opt)
hdr = {"HD":{"VN":"1.6"},"SQ":[{"SN":"chr20","LN":64444167}]}
with pysam.AlignmentFile(out,"wb",header=hdr) as o:
    for i,recs in enumerate(res):
        for r in recs:
            r.flag |= 0x40 if i%2==0 else 0x80
            o.write(r)
