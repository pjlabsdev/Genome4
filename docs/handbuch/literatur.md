# Literatur

Quellen der Verfahren, die Genome 4.0 umsetzt. Vor einer Veröffentlichung werden alle Angaben
noch einmal gegen die Originale geprüft.

| Kürzel | Quelle | Wofür |
|---|---|---|
| Burrows & Wheeler 1994 | Burrows M, Wheeler DJ. A block-sorting lossless data compression algorithm. Technical Report 124, Digital Equipment Corporation, 1994. | BWT |
| Ferragina & Manzini 2000 | Ferragina P, Manzini G. Opportunistic data structures with applications. FOCS 2000. | FM-Index |
| Nong, Zhang & Chan 2009 | Nong G, Zhang S, Chan WH. Linear suffix array construction by almost pure induced-sorting. Data Compression Conference 2009. | SA-IS |
| Burkhardt & Kärkkäinen 2003 | Burkhardt S, Kärkkäinen J. Fast lightweight suffix array construction and checking. CPM 2003. | Differenzüberdeckung, Prüfung des Suffix-Arrays |
| Colbourn & Ling 2000 | Colbourn CJ, Ling ACH. Quorums from difference covers. Information Processing Letters 75, 2000. | Konstruktion der Differenzüberdeckung |
| Kärkkäinen 2007 | Kärkkäinen J. Fast BWT in small space by blockwise suffix sorting. Theoretical Computer Science 387, 2007. | blockweiser Bau |
| Larsson & Sadakane 2007 | Larsson NJ, Sadakane K. Faster suffix sorting. Theoretical Computer Science 387, 2007. | Verdoppeln beim Sortieren der Stichprobe |
| Li & Durbin 2009 | Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25:1754-1760, 2009. | BWA, Index über beide Stränge, `bwa aln` als Vergleich |
| Gotoh 1982 | Gotoh O. An improved algorithm for matching biological sequences. Journal of Molecular Biology 162:705-708, 1982. | affine Lücken in quadratischer Zeit |
| Li 2012 | Li H. Exploring single-sample SNP and INDEL calling with whole-genome de novo assembly. Bioinformatics 28:1838-1844, 2012. | FMD-Index: Suche in beide Richtungen mit einer BWT |
| Farrar 2007 | Farrar M. Striped Smith-Waterman speeds database searches six times over other SIMD implementations. Bioinformatics 23:156-161, 2007. | lokales Alignment in gestreifter Anordnung (Partnersuche bei Paaren) |
| Li 2013 | Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997, 2013. | BWA-MEM (ab M3) |
| Vasimuddin et al. 2019 | Vasimuddin M, Misra S, Li H, Aluru S. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. IPDPS 2019. | BWA-MEM2 als Tempo-Vergleich |
| Knuth 1998 | Knuth DE. The Art of Computer Programming, Volume 3: Sorting and Searching. 2. Auflage, Addison-Wesley, 1998. Abschnitt 5.4 (externes Sortieren). | Sortieren in Teilen |
| Li et al. 2009 | Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25:2078-2079, 2009. | SAM/BAM |
| SAM-Spezifikation | The SAM/BAM Format Specification Working Group. Sequence Alignment/Map Format Specification. hts-specs. | BAM und BGZF im Einzelnen |

## Quelltexte, an denen der Nachbau geprüft wird
- **bwa** (Heng Li), `bwamem.c`, `bwamem_pair.c`, `ksw.c`, `ksort.h`, `kbtree.h`: MIT-Lizenz laut Dateikopf. Der
  Nachbau in `gc-align` ist eigener Rust-Code; wo er dem Ablauf von bwa Schritt für Schritt
  folgt, steht der Name der bwa-Funktion im Quelltext und auf der Verfahrensseite. Vor einer
  Weitergabe von Genome 4.0 ist der MIT-Hinweis von bwa und klib beizulegen.
