#!/usr/bin/env python3
"""Ruft die Paar-Funktionen von bwa einzeln auf (aus der Bibliothek, die pybwa mitbringt) und
schreibt Zwischenstände, um Abweichungen einer Stufe zuzuordnen.

Je Paar: Bereiche beider Reads aus mem_align1_core (Zeile SIG: Zahl, Summe der Anfänge, Summe
der Scores), dann mem_matesw, mem_mark_primary_se und mem_pair wie in mem_sam_pe (Zeile PAAR:
Score, zweitbester, Zahl fast gleich guter, Anfänge der gewählten Bereiche). Die Insertgröße
(Richtung FR) wird vorgegeben, die Paarnummer zählt ab 0.
Aufruf: trace_bwa_pairs.py bwa-prefix r1.fq r2.fq von bis mittel streuung > spur.txt
"""
import ctypes as C, glob, os, sys, pybwa, pysam
lib = C.CDLL(glob.glob(os.path.join(os.path.dirname(pybwa.__file__), "libbwamem*.so"))[0])
class Reg(C.Structure):
    _fields_=[("rb",C.c_int64),("re",C.c_int64),("qb",C.c_int),("qe",C.c_int),("rid",C.c_int),("score",C.c_int),("truesc",C.c_int),("sub",C.c_int),("alt_sc",C.c_int),("csub",C.c_int),("sub_n",C.c_int),("w",C.c_int),("seedcov",C.c_int),("secondary",C.c_int),("secondary_all",C.c_int),("seedlen0",C.c_int),("bits",C.c_int),("frac_rep",C.c_float),("hash",C.c_uint64)]
class RegV(C.Structure):
    _fields_=[("n",C.c_size_t),("m",C.c_size_t),("a",C.POINTER(Reg))]
class Pes(C.Structure):
    _fields_=[("low",C.c_int),("high",C.c_int),("failed",C.c_int),("avg",C.c_double),("std",C.c_double)]
class Idx(C.Structure):
    _fields_=[("bwt",C.c_void_p),("bns",C.c_void_p),("pac",C.c_void_p)]
lib.bwa_idx_load.restype=C.POINTER(Idx); lib.bwa_idx_load.argtypes=[C.c_char_p,C.c_int]
lib.mem_opt_init.restype=C.c_void_p
lib.mem_align1_core.restype=RegV; lib.mem_align1_core.argtypes=[C.c_void_p]*4+[C.c_int,C.c_char_p,C.c_void_p]
lib.mem_matesw.restype=C.c_int; lib.mem_matesw.argtypes=[C.c_void_p]*3+[C.POINTER(Pes),C.POINTER(Reg),C.c_int,C.c_char_p,C.POINTER(RegV)]
lib.mem_mark_primary_se.restype=C.c_int; lib.mem_mark_primary_se.argtypes=[C.c_void_p,C.c_int,C.POINTER(Reg),C.c_int64]
lib.mem_pair.restype=C.c_int; lib.mem_pair.argtypes=[C.c_void_p]*3+[C.POINTER(Pes),C.c_void_p,C.POINTER(RegV),C.c_int,C.POINTER(C.c_int),C.POINTER(C.c_int),C.POINTER(C.c_int),C.POINTER(C.c_int)]
prefix, fq1, fq2 = sys.argv[1:4]
low, high, avg, std = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
idx=lib.bwa_idx_load(prefix.encode(),7).contents
opt=lib.mem_opt_init()
nt={"A":0,"C":1,"G":2,"T":3}
pes=(Pes*4)()
for d in range(4): pes[d].failed=1
pes[1]=Pes(low,high,0,avg,std)
for pid,(q1,q2) in enumerate(zip(pysam.FastxFile(fq1),pysam.FastxFile(fq2))):
    qs=[q1.sequence,q2.sequence]; name=q1.name
    seqs=[C.create_string_buffer(bytes(nt.get(c,4) for c in q),len(q)) for q in qs]
    a=(RegV*2)()
    for i in range(2):
        a[i]=lib.mem_align1_core(opt,idx.bwt,idx.bns,idx.pac,len(qs[i]),seqs[i],None)
        rs=[a[i].a[j] for j in range(a[i].n)]
        print("SIG",name,i,len(rs),sum(r.rb for r in rs),sum(r.score for r in rs))
    b=[[Reg.from_buffer_copy(a[i].a[j]) for j in range(a[i].n) if a[i].a[j].score>=a[i].a[0].score-17] for i in range(2)]
    for i in range(2):
        for j in range(min(len(b[i]),50)):
            lib.mem_matesw(opt,idx.bns,idx.pac,pes,C.byref(b[i][j]),len(qs[1-i]),seqs[1-i],C.byref(a[1-i]))
    npri=(C.c_int*2)()
    for i in range(2):
        npri[i]=lib.mem_mark_primary_se(opt,a[i].n,a[i].a,pid<<1|i) if a[i].n else 0
    sub=C.c_int();nsub=C.c_int();z=(C.c_int*2)()
    o=lib.mem_pair(opt,idx.bns,idx.pac,pes,None,a,pid,C.byref(sub),C.byref(nsub),z,npri) if npri[0] and npri[1] else 0
    if o>0: print("PAAR",name,o,sub.value,nsub.value,a[0].a[z[0]].rb,a[1].a[z[1]].rb)
    else: print("PAAR",name,"-")
