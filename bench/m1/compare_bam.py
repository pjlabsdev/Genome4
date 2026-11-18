"""Vergleicht zwei BAM-Dateien über samtools/htslib (pysam): Kopftext und jeder Datensatz.

Der Prüfer ist bewusst ein fremdes Werkzeug, damit er vom Geprüften unabhängig ist.
Aufruf: python compare_bam.py original.bam neu.bam   (Exit 0 = identisch)
"""
import pysam, hashlib, sys
def digest(p):
    h=hashlib.sha256(); n=0
    with pysam.AlignmentFile(p,"rb",threads=4) as f:
        h.update(str(f.header).encode())
        for r in f:
            h.update(r.to_string().encode()); h.update(b"\n"); n+=1
    return n,h.hexdigest()
a=digest(sys.argv[1]); b=digest(sys.argv[2]); print(a); print(b); print("IDENTICAL" if a==b else "DIFFERENT")
sys.exit(0 if a==b else 1)
