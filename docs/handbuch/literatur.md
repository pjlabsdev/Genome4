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
| Li 2012 | Li H. Exploring single-sample SNP and INDEL calling with whole-genome de novo assembly. Bioinformatics 28:1838-1844, 2012. | FMD-Index: Suche in beide Richtungen mit einer BWT |
| Li 2013 | Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997, 2013. | BWA-MEM (ab M3) |
| Vasimuddin et al. 2019 | Vasimuddin M, Misra S, Li H, Aluru S. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. IPDPS 2019. | BWA-MEM2 als Tempo-Vergleich |
| Li et al. 2009 | Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25:2078-2079, 2009. | SAM/BAM |
| SAM-Spezifikation | The SAM/BAM Format Specification Working Group. Sequence Alignment/Map Format Specification. hts-specs. | BAM und BGZF im Einzelnen |
