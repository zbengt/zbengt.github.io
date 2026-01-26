---
layout: post
title: "E5 Deep-Dive Expression homologous lncRNA sequences evaluation"
date: 2026-01-26
tags: [lncRNA, deep-dive]
categories: [coral, comparative-genomics]
---

Code repo link: [13.1-lncRNA-cross-species.Rmd](https://github.com/urol-e5/deep-dive-expression/blob/main/M-multi-species/code/13.1-lncRNA-cross-species.Rmd)

# Summary

Using resulting lncRNAs from reciprocal blasts under my relaxed criteria (≥50% identity, ≥30% coverage, ≤ 1e-5) while still enforcing 1:1 and 1:1:1 relationships, I created FASTAs of highly similar sequences across species and BLASTed them against NCBI RefSeq RNA.

-   **General Trends**: A handful of sequences in had no hits. Most sequences had reasonably high confidence hits with uncharaterized ncRNAs previously annotated as such in _Acropora_, _Montipora_, _Pocillopora_, and _Porites_ transcriptomes. This is a good sign that our pipeline did indeed identify putative lncRNAs. However, this isn't especially helpful in pointing toward potential functions of the lncRNAs we detected.

-   **Interesting Matches**: A small number of lncRNAs actually hit previously annotated mRNA sequences of known potential protein coding capabilities. Most notable of these was the hit for [**Nematostella vectensis transcriptional repressor NF-X1**](https://www.ncbi.nlm.nih.gov/nucleotide/XM_048723524.1) in our all 3 species match category, the only sequence showing stringent homology across the 3 species. This is especially interesting given NF-X1's role in transcriptional regulation, one of lncRNA's most supported functions. It could be that a lncRNA implicated in NF-X1 transcriptional represssion intersects with this gene and is involved in this repression pathway (see caution bullet point below for why this may not be the case). This same line of reasoning might be applied to the other small number of lncRNAs that hit functionally annotated mRNA sequences. Hits like [**Montipora capricornis cytidine and dCMP deaminase domain-containing protein 1-like**](https://www.ncbi.nlm.nih.gov/nucleotide/XM_068862579.1?report=genbank&log$=nuclalign&blast_rank=17&RID=RFRPVNCG014) for Apul-Peve and [**Pocillopora verrucosa protein fem-1 homolog B-like**](https://www.ncbi.nlm.nih.gov/nucleotide/XM_066169535.1?report=genbank&log$=nuclalign&blast_rank=3&RID=RFUZ4RCV016) for Peve-Ptuh are worth further exploration as well. **Generally, it seems like functionally annotated hits fit within gene regulation and cellular maintanence categories, which makes sense to me since we would expect lncRNAs in these basic/foundational categories to be the most conserved if conservation is possible.**

-   **Caution**: It's quite possible I am seeing hits classified as mRNAs because we did not successfully eliminate all mRNAs during our lncRNA discovery pipeline. I ran all putative homologs through CPC2 again to visually confirm whether or not any sequences predicted to code for proteins slipped through. None of the putative homologs were flagged as having coding potential, but I'm not sure I am convinced, particularly given the ORF integrity results from CPC2 indicating that at least some of the sequences begin with a valid start codon and end with an in-frame stop codon. The other components evaluated by CPC2 seem to overrule this coding potential metric, with putative peptide length, Fickett score, and isoeletric point being outside the normal distributions expected for protein coding sequences. This definitely merits further consideration.

Below is an example of BLAST refseq RNA results, specifically for the sequence found across all 3 species. As well as tables summarizing CPC2 analysis of coding potential for every proposed homolog.

# Next steps

-   Determine whether porposed lncRNA homologs match annotated mRNAs
-   Examine position and potential intersections of lncRNAs
-   Look at co-expression results for these lncRNAs

# RefSeq RNA BLAST Results

Ran using refseq_rna database in NCBI standard databases. Optimized for highly similar sequences (megablast).

### Apul–Peve–Ptuh

Apul sequences:

| Description                                                 | Scientific Name          | Max Score | Total Score | Query Cover |  E value | Per. ident | Acc. Len | Accession                                                                |
|---------|--------|-------:|-------:|-------:|-------:|-------:|-------:|-------------|
| PREDICTED: Acropora muricata uncharacterized lncRNA         | Acropora muricata        |     455.0 |       455.0 |        100% | 4.0e-124 |      96.38 |     1711 | [XR_010871610.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_010871610.1) |
| PREDICTED: Nematostella vectensis transcriptional regulator | Nematostella vectensis   |     224.0 |       224.0 |         95% |  1.0e-54 |      82.13 |     4373 | [XM_048723524.1](https://www.ncbi.nlm.nih.gov/nucleotide/XM_048723524.1) |
| PREDICTED: Acropora digitifera uncharacterized transcript   | Acropora digitifera      |     148.0 |       148.0 |         39% |  8.0e-32 |      91.59 |      220 | [XR_001561218.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_001561218.1) |
| PREDICTED: Saccoglossus kowalevskii uncharacterized gene    | Saccoglossus kowalevskii |     115.0 |       115.0 |         79% |  8.0e-22 |      76.47 |      280 | [XR_013422149.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_013422149.1) |
| PREDICTED: Acropora digitifera uncharacterized transcript   | Acropora digitifera      |      93.5 |        93.5 |         20% |  4.0e-15 |      96.49 |      102 | [XR_001564892.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_001564892.1) |
| PREDICTED: Pomacea canaliculata uncharacterized transcript  | Pomacea canaliculata     |      92.7 |        92.7 |         17% |  5.0e-14 |      74.36 |     6106 | [XM_025259772.1](https://www.ncbi.nlm.nih.gov/nucleotide/XM_025259772.1) |
| PREDICTED: Pomacea canaliculata uncharacterized transcript  | Pomacea canaliculata     |      92.7 |        92.7 |         17% |  5.0e-14 |      74.36 |     6410 | [XM_025259773.1](https://www.ncbi.nlm.nih.gov/nucleotide/XM_025259773.1) |
| PREDICTED: Saccoglossus kowalevskii uncharacterized gene    | Saccoglossus kowalevskii |      72.4 |        72.4 |         29% |  4.0e-10 |      91.07 |      249 | [XR_013422150.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_013422150.1) |

Peve sequences:

| Description                                                 | Scientific Name          | Max Score | Total Score | Query Cover | E value | Per. ident | Acc. Len | Accession                                                                |
|---------|-------|------:|------:|------:|------:|------:|------:|-------------|
| PREDICTED: Acropora muricata uncharacterized lncRNA         | Acropora muricata        |     346.0 |       346.0 |         98% | 3.0e-91 |      89.42 |     1711 | [XR_010871610.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_010871610.1) |
| PREDICTED: Nematostella vectensis transcriptional regulator | Nematostella vectensis   |     278.0 |       278.0 |         99% | 1.0e-70 |      85.00 |     4373 | [XM_048723524.1](https://www.ncbi.nlm.nih.gov/nucleotide/XM_048723524.1) |
| PREDICTED: Saccoglossus kowalevskii uncharacterized gene    | Saccoglossus kowalevskii |     139.0 |       139.0 |         85% | 5.0e-29 |      77.87 |      280 | [XR_013422149.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_013422149.1) |
| PREDICTED: Pomacea canaliculata uncharacterized transcript  | Pomacea canaliculata     |     134.0 |       134.0 |         91% | 2.0e-27 |      76.47 |     6106 | [XM_025259772.1](https://www.ncbi.nlm.nih.gov/nucleotide/XM_025259772.1) |
| PREDICTED: Pomacea canaliculata uncharacterized transcript  | Pomacea canaliculata     |     134.0 |       134.0 |         91% | 2.0e-27 |      76.47 |     6410 | [XM_025259773.1](https://www.ncbi.nlm.nih.gov/nucleotide/XM_025259773.1) |
| PREDICTED: Acropora digitifera uncharacterized transcript   | Acropora digitifera      |     130.0 |       130.0 |         37% | 3.0e-26 |      89.32 |      220 | [XR_001561218.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_001561218.1) |
| PREDICTED: Saccoglossus kowalevskii uncharacterized gene    | Saccoglossus kowalevskii |      82.4 |        82.4 |         41% | 8.0e-12 |      80.17 |      249 | [XR_013422150.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_013422150.1) |

Ptuh sequences:

| Description                                                 | Scientific Name          | Max Score | Total Score | Query Cover | E value | Per. ident | Acc. Len | Accession                                                                |
|---------|-------|------:|------:|------:|------:|------:|------:|-------------|
| PREDICTED: Acropora muricata uncharacterized lncRNA         | Acropora muricata        |     302.0 |       302.0 |        100% | 6.0e-78 |      86.99 |     1711 | [XR_010871610.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_010871610.1) |
| PREDICTED: Nematostella vectensis transcriptional regulator | Nematostella vectensis   |     259.0 |       259.0 |         94% | 3.0e-65 |      85.10 |     4373 | [XM_048723524.1](https://www.ncbi.nlm.nih.gov/nucleotide/XM_048723524.1) |
| PREDICTED: Saccoglossus kowalevskii uncharacterized gene    | Saccoglossus kowalevskii |     148.0 |       148.0 |         92% | 8.0e-32 |      77.69 |      280 | [XR_013422149.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_013422149.1) |
| PREDICTED: Pomacea canaliculata uncharacterized transcript  | Pomacea canaliculata     |     115.0 |       115.0 |         86% | 8.0e-22 |      76.07 |     6106 | [XM_025259772.1](https://www.ncbi.nlm.nih.gov/nucleotide/XM_025259772.1) |
| PREDICTED: Pomacea canaliculata uncharacterized transcript  | Pomacea canaliculata     |     115.0 |       115.0 |         86% | 8.0e-22 |      76.07 |     6410 | [XM_025259773.1](https://www.ncbi.nlm.nih.gov/nucleotide/XM_025259773.1) |
| PREDICTED: Acropora digitifera uncharacterized transcript   | Acropora digitifera      |      99.0 |        99.0 |         25% | 8.0e-17 |      92.65 |      220 | [XR_001561218.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_001561218.1) |
| PREDICTED: Saccoglossus kowalevskii uncharacterized gene    | Saccoglossus kowalevskii |      80.5 |        80.5 |         24% | 3.0e-11 |      89.23 |      249 | [XR_013422150.1](https://www.ncbi.nlm.nih.gov/nucleotide/XR_013422150.1) |

# Coding Potential Calculator 2 Summaries

Used CPC2 to double check coding potential and gather additional statistics on sequences.

### Apul–Peve–Ptuh (CPC2 results)

| ID                | peptide_length | Fickett_score |     pI | ORF_integrity | coding_probability | label     |
|-----------|----------:|----------:|----------:|----------:|----------:|-----------|
| Apul_lncRNA_16008 |             41 |       0.31522 | 4.4278 |             1 |          0.0225572 | noncoding |
| Peve_lncRNA_2714  |              0 |       0.34808 | 0.0000 |            -1 |       0.0000011485 | noncoding |
| Ptuh_lncRNA_1396  |             37 |       0.35813 | 9.6878 |             1 |          0.0146368 | noncoding |

### Peve–Ptuh (CPC2 results)

| ID                | peptide_length | Fickett_score |      pI | ORF_integrity | coding_probability | label     |
|-----------|----------:|----------:|----------:|----------:|----------:|-----------|
| Peve_lncRNA_2367  |             70 |       0.40302 |  6.4756 |            -1 |          0.0324823 | noncoding |
| Peve_lncRNA_2714  |              0 |       0.34808 |  0.0000 |            -1 |       0.0000011485 | noncoding |
| Peve_lncRNA_2976  |              7 |       0.40956 |  6.4910 |             1 |       0.0000125227 | noncoding |
| Peve_lncRNA_3463  |             58 |       0.37421 |  8.5521 |             1 |          0.0326739 | noncoding |
| Peve_lncRNA_4043  |            104 |       0.47052 |  9.4664 |            -1 |           0.122161 | noncoding |
| Peve_lncRNA_4250  |             49 |       0.32454 | 10.6385 |             1 |          0.0253991 | noncoding |
| Peve_lncRNA_6079  |             61 |       0.39530 | 10.7259 |             1 |          0.0644281 | noncoding |
| Peve_lncRNA_7236  |             56 |       0.40038 |  7.8256 |             1 |          0.0348187 | noncoding |
| Peve_lncRNA_7779  |             21 |       0.41468 |  5.7503 |             1 |         0.00591771 | noncoding |
| Peve_lncRNA_8456  |             52 |       0.40982 | 10.3061 |             1 |          0.0496757 | noncoding |
| Peve_lncRNA_8599  |             17 |       0.46386 | 11.0002 |             1 |          0.0239428 | noncoding |
| Peve_lncRNA_8863  |             18 |       0.35377 |  5.2191 |             1 |         0.00406401 | noncoding |
| Peve_lncRNA_9076  |             53 |       0.33811 |  9.4166 |             1 |          0.0239105 | noncoding |
| Peve_lncRNA_9286  |             15 |       0.34200 |  5.0606 |            -1 |          0.0170078 | noncoding |
| Ptuh_lncRNA_14108 |             73 |       0.40421 |  4.8149 |             1 |           0.205848 | noncoding |
| Ptuh_lncRNA_1396  |             37 |       0.35813 |  9.6878 |             1 |          0.0146368 | noncoding |
| Ptuh_lncRNA_921   |             59 |       0.40536 |  9.3326 |             1 |          0.0545046 | noncoding |
| Ptuh_lncRNA_6242  |             26 |       0.35867 |  8.9396 |             1 |         0.00711255 | noncoding |
| Ptuh_lncRNA_7905  |             26 |       0.33958 |  4.5296 |             1 |         0.00855639 | noncoding |
| Ptuh_lncRNA_11747 |             77 |       0.26715 |  8.4809 |             1 |          0.0438526 | noncoding |
| Ptuh_lncRNA_1654  |            177 |       0.50404 |  6.0023 |            -1 |           0.185281 | noncoding |
| Ptuh_lncRNA_1498  |             66 |       0.37567 |  8.1115 |             1 |          0.0467746 | noncoding |
| Ptuh_lncRNA_13467 |             21 |       0.41610 |  5.7503 |             1 |         0.00597448 | noncoding |
| Ptuh_lncRNA_3898  |             49 |       0.37818 | 11.0052 |             1 |          0.0339989 | noncoding |
| Ptuh_lncRNA_12651 |             33 |       0.46538 |  9.1005 |            -1 |          0.0442831 | noncoding |
| Ptuh_lncRNA_11868 |             28 |       0.37089 | 10.0285 |             1 |           0.011644 | noncoding |
| Ptuh_lncRNA_1486  |             53 |       0.32219 |  6.8935 |             1 |          0.0174918 | noncoding |
| Ptuh_lncRNA_4161  |             28 |       0.38118 |  9.1027 |             1 |         0.00905463 | noncoding |

### Apul–Ptuh (CPC2 results)

| ID                | peptide_length | Fickett_score |          pI | ORF_integrity | coding_probability | label     |
| ----------------- | -------------: | ------------: | ----------: | ------------: | -----------------: | --------- |
| Apul_lncRNA_11205 |            288 |       0.43082 | 9.456359863 |            -1 |           0.460064 | noncoding |
| Apul_lncRNA_12375 |             71 |       0.33092 | 4.601013184 |             1 |           0.103221 | noncoding |
| Apul_lncRNA_12954 |             22 |       0.38628 | 9.862487793 |             1 |          0.0090556 | noncoding |
| Apul_lncRNA_13651 |             31 |       0.34779 | 11.26373291 |             1 |          0.0158401 | noncoding |
| Apul_lncRNA_13665 |             67 |       0.37510 | 10.37664795 |             1 |          0.0649536 | noncoding |
| Apul_lncRNA_14079 |            113 |       0.49597 | 6.636901855 |            -1 |          0.0654743 | noncoding |
| Apul_lncRNA_14131 |             47 |       0.32934 | 7.837097168 |             1 |          0.0133536 | noncoding |
| Apul_lncRNA_14842 |             46 |       0.35722 |  10.1774292 |             1 |          0.0238125 | noncoding |
| Apul_lncRNA_16008 |             41 |       0.31522 |  4.42779541 |             1 |          0.0225572 | noncoding |
| Apul_lncRNA_16843 |             95 |       0.31484 | 9.738830566 |             1 |          0.0935923 | noncoding |
| Apul_lncRNA_17663 |             23 |       0.31636 | 5.274719238 |             1 |         0.00565893 | noncoding |
| Apul_lncRNA_18874 |             83 |       0.34204 | 9.805969238 |             1 |          0.0818223 | noncoding |
| Apul_lncRNA_19225 |             97 |       0.32626 | 9.181945801 |             1 |           0.111006 | noncoding |
| Apul_lncRNA_19808 |             71 |       0.37986 | 9.497619629 |             1 |          0.0731677 | noncoding |
| Apul_lncRNA_20865 |             54 |       0.42265 | 10.58270264 |            -1 |           0.122738 | noncoding |
| Apul_lncRNA_20866 |             23 |       0.31097 | 9.997497559 |            -1 |           0.178024 | noncoding |
| Apul_lncRNA_23028 |             30 |       0.37802 |  9.30279541 |             1 |          0.0104584 | noncoding |
| Apul_lncRNA_23865 |             82 |       0.38713 | 5.728210449 |             1 |           0.179713 | noncoding |
| Apul_lncRNA_24419 |             51 |       0.37101 |  4.82623291 |             1 |          0.0398573 | noncoding |
| Apul_lncRNA_25911 |             26 |       0.33362 | 8.938293457 |             1 |         0.00687062 | noncoding |
| Apul_lncRNA_26325 |             44 |       0.37237 |  7.97833252 |            -1 |          0.0506618 | noncoding |
| Apul_lncRNA_26468 |             77 |       0.37774 | 8.299987793 |             1 |          0.0828723 | noncoding |
| Apul_lncRNA_26864 |             56 |       0.34432 | 10.73809814 |             1 |          0.0348038 | noncoding |
| Apul_lncRNA_26928 |             62 |       0.32002 | 10.40802002 |             1 |          0.0363458 | noncoding |
| Apul_lncRNA_2721  |             29 |       0.45295 | 12.00006104 |             1 |           0.033624 | noncoding |
| Apul_lncRNA_27421 |             84 |       0.34200 | 9.692565918 |             1 |          0.0841464 | noncoding |
| Apul_lncRNA_27481 |             58 |       0.37752 | 5.460876465 |             1 |          0.0474245 | noncoding |
| Apul_lncRNA_27488 |             77 |       0.34629 | 10.50970459 |             1 |          0.0715875 | noncoding |
| Apul_lncRNA_27904 |             44 |       0.31880 | 9.097961426 |             1 |          0.0146653 | noncoding |
| Apul_lncRNA_28077 |             42 |       0.41225 | 5.188415527 |            -1 |          0.0177318 | noncoding |
| Apul_lncRNA_288   |             89 |       0.32577 | 5.410339355 |             1 |           0.162289 | noncoding |
| Apul_lncRNA_28932 |             67 |       0.45524 | 10.85906982 |             1 |           0.206445 | noncoding |
| Apul_lncRNA_29686 |            103 |       0.33026 | 9.061828613 |             1 |           0.143293 | noncoding |
| Apul_lncRNA_31231 |             13 |       0.37587 | 5.272766113 |             1 |         0.00324109 | noncoding |
| Apul_lncRNA_4002  |             41 |       0.37606 | 6.887878418 |             1 |          0.0121214 | noncoding |
| Apul_lncRNA_7519  |             39 |       0.37546 | 9.624694824 |             1 |          0.0172243 | noncoding |
| Apul_lncRNA_8146  |             98 |       0.35079 | 11.49920654 |             1 |           0.126951 | noncoding |
| Ptuh_lncRNA_3186  |             28 |       0.29928 | 8.254211426 |             1 |         0.00672497 | noncoding |
| Ptuh_lncRNA_3517  |             67 |       0.29244 | 9.724060059 |             1 |          0.0351263 | noncoding |
| Ptuh_lncRNA_16011 |             51 |       0.36620 | 5.871276855 |             1 |          0.0242798 | noncoding |
| Ptuh_lncRNA_12658 |             22 |       0.40225 | 3.998840332 |             1 |          0.0145898 | noncoding |
| Ptuh_lncRNA_11921 |             25 |       0.34918 | 9.496276855 |             1 |         0.00809292 | noncoding |
| Ptuh_lncRNA_1260  |             50 |       0.31328 | 5.737731934 |             1 |          0.0199501 | noncoding |
| Ptuh_lncRNA_14951 |             12 |       0.39832 | 12.01654053 |            -1 |           0.163744 | noncoding |
| Ptuh_lncRNA_11816 |             58 |       0.36569 | 8.988342285 |             1 |          0.0327129 | noncoding |
| Ptuh_lncRNA_1396  |             37 |       0.35813 | 9.687805176 |             1 |          0.0146368 | noncoding |
| Ptuh_lncRNA_1065  |            137 |       0.29594 | 10.05877686 |             1 |           0.231113 | noncoding |
| Ptuh_lncRNA_12620 |              6 |       0.30716 | 5.274963379 |             1 |          0.0025279 | noncoding |
| Ptuh_lncRNA_1279  |            136 |       0.31476 | 9.882873535 |             1 |           0.293585 | noncoding |
| Ptuh_lncRNA_5958  |             42 |       0.44463 | 7.944030762 |            -1 |          0.0303491 | noncoding |
| Ptuh_lncRNA_14029 |             52 |       0.44023 |  8.84197998 |             1 |          0.0519598 | noncoding |
| Ptuh_lncRNA_8435  |             41 |       0.44615 | 8.363464355 |             1 |           0.026238 | noncoding |
| Ptuh_lncRNA_7446  |             41 |       0.36506 | 8.363464355 |             1 |           0.012787 | noncoding |
| Ptuh_lncRNA_6767  |             43 |       0.32626 | 5.996154785 |             1 |           0.012574 | noncoding |
| Ptuh_lncRNA_7862  |             66 |       0.36763 | 12.01141357 |            -1 |            0.17804 | noncoding |
| Ptuh_lncRNA_11175 |             58 |       0.39670 | 6.554016113 |            -1 |           0.028832 | noncoding |
| Ptuh_lncRNA_9962  |              8 |       0.36590 | 6.490661621 |             1 |        9.79939e-06 | noncoding |
| Ptuh_lncRNA_1652  |             36 |       0.43633 | 10.06048584 |             1 |          0.0305934 | noncoding |
| Ptuh_lncRNA_5429  |             18 |       0.34526 | 9.700622559 |             1 |         0.00653488 | noncoding |
| Ptuh_lncRNA_10359 |             45 |       0.36998 |  8.80291748 |             1 |          0.0176839 | noncoding |
| Ptuh_lncRNA_5732  |            103 |       0.34787 | 10.01275635 |             1 |           0.175399 | noncoding |
| Ptuh_lncRNA_4953  |             39 |       0.34425 | 8.494445801 |             1 |          0.0108628 | noncoding |
| Ptuh_lncRNA_4630  |             48 |       0.43806 | 4.983825684 |            -1 |          0.0220361 | noncoding |
| Ptuh_lncRNA_1208  |             53 |       0.38010 |  9.77911377 |             1 |          0.0347193 | noncoding |
| Ptuh_lncRNA_10862 |             93 |       0.41551 | 10.05499268 |             1 |           0.294029 | noncoding |
| Ptuh_lncRNA_14402 |             70 |       0.39857 | 10.02850342 |            -1 |           0.118543 | noncoding |
| Ptuh_lncRNA_4981  |             14 |       0.25772 | 4.298156738 |             1 |          0.0235857 | noncoding |
| Ptuh_lncRNA_2723  |             61 |       0.29543 | 9.886901855 |             1 |          0.0304174 | noncoding |
| Ptuh_lncRNA_12289 |             83 |       0.38343 | 7.762512207 |             1 |           0.117721 | noncoding |
| Ptuh_lncRNA_1708  |             92 |       0.40564 | 8.638122559 |             1 |           0.232939 | noncoding |
| Ptuh_lncRNA_11007 |             10 |       0.37905 | 9.785217285 |             1 |         0.00493444 | noncoding |
| Ptuh_lncRNA_4599  |             93 |       0.32738 | 8.361877441 |             1 |          0.0967391 | noncoding |
| Ptuh_lncRNA_11883 |            113 |       0.32918 | 10.20269775 |             1 |           0.186662 | noncoding |
| Ptuh_lncRNA_7018  |             67 |       0.33099 | 9.995178223 |             1 |          0.0430005 | noncoding |

### Apul–Peve (CPC2 results)

| ID                | peptide_length | Fickett_score |          pI | ORF_integrity | coding_probability | label     |
| ----------------- | -------------: | ------------: | ----------: | ------------: | -----------------: | --------- |
| Apul_lncRNA_10150 |             35 |       0.36883 | 4.749206543 |             1 |          0.0152901 | noncoding |
| Apul_lncRNA_11023 |             52 |       0.33416 | 9.302307129 |             1 |          0.0219726 | noncoding |
| Apul_lncRNA_12170 |             52 |       0.35497 |  7.79510498 |             1 |          0.0190474 | noncoding |
| Apul_lncRNA_12365 |            100 |       0.39199 | 6.666564941 |             1 |           0.327169 | noncoding |
| Apul_lncRNA_12639 |             75 |       0.34107 |  6.88885498 |             1 |          0.0561701 | noncoding |
| Apul_lncRNA_13011 |             93 |       0.41674 | 6.268981934 |            -1 |          0.0429362 | noncoding |
| Apul_lncRNA_13934 |             19 |       0.32935 | 10.28900146 |             1 |         0.00864502 | noncoding |
| Apul_lncRNA_14577 |            103 |       0.33642 | 10.96063232 |             1 |           0.140515 | noncoding |
| Apul_lncRNA_14938 |             33 |       0.43317 | 10.83428955 |            -1 |           0.114709 | noncoding |
| Apul_lncRNA_15055 |             15 |       0.36090 | 12.00006104 |             1 |         0.00907571 | noncoding |
| Apul_lncRNA_16008 |             41 |       0.31522 |  4.42779541 |             1 |          0.0225572 | noncoding |
| Apul_lncRNA_16554 |             38 |       0.44024 |  10.3838501 |             1 |          0.0390687 | noncoding |
| Apul_lncRNA_17819 |             38 |       0.45336 | 10.04425049 |             1 |          0.0436952 | noncoding |
| Apul_lncRNA_17971 |              6 |       0.35409 | 5.274963379 |             1 |        1.15257e-05 | noncoding |
| Apul_lncRNA_20265 |             23 |       0.39772 | 9.186706543 |             1 |         0.00801587 | noncoding |
| Apul_lncRNA_21518 |             99 |       0.37418 | 6.892272949 |             1 |           0.247124 | noncoding |
| Apul_lncRNA_21555 |             81 |       0.37634 | 9.124572754 |             1 |           0.104222 | noncoding |
| Apul_lncRNA_22097 |             31 |       0.24527 | 10.27740479 |             1 |          0.0267043 | noncoding |
| Apul_lncRNA_25307 |            114 |       0.40949 | 11.48248291 |             1 |           0.429807 | noncoding |
| Apul_lncRNA_25972 |            114 |       0.32294 | 9.450134277 |             1 |           0.186137 | noncoding |
| Apul_lncRNA_26197 |             43 |       0.44402 | 10.84283447 |             1 |          0.0575382 | noncoding |
| Apul_lncRNA_26646 |            111 |       0.32475 | 5.200256348 |            -1 |          0.0717744 | noncoding |
| Apul_lncRNA_27113 |             48 |       0.42754 | 5.212097168 |            -1 |          0.0205513 | noncoding |
| Apul_lncRNA_27723 |            105 |       0.34839 | 8.487487793 |            -1 |           0.115908 | noncoding |
| Apul_lncRNA_30304 |             96 |       0.31136 | 6.375305176 |             1 |           0.138713 | noncoding |
| Apul_lncRNA_3622  |             77 |       0.40465 | 10.08538818 |             1 |           0.136416 | noncoding |
| Apul_lncRNA_3835  |             35 |       0.31898 | 9.224060059 |             1 |          0.0108196 | noncoding |
| Apul_lncRNA_5253  |            115 |       0.35665 | 8.667785645 |             1 |            0.30346 | noncoding |
| Apul_lncRNA_8970  |             21 |       0.43659 | 9.628356934 |             1 |          0.0120618 | noncoding |
| Peve_lncRNA_7892  |             28 |       0.34094 |  9.69720459 |             1 |         0.00971246 | noncoding |
| Peve_lncRNA_7749  |             10 |       0.36231 | 3.998840332 |             1 |         0.00459624 | noncoding |
| Peve_lncRNA_9439  |             25 |       0.32854 | 9.778991699 |             1 |         0.00893162 | noncoding |
| Peve_lncRNA_8466  |             51 |       0.34268 | 9.694885254 |             1 |          0.0241012 | noncoding |
| Peve_lncRNA_282   |            101 |       0.28950 | 9.409484863 |             1 |          0.0906865 | noncoding |
| Peve_lncRNA_2752  |             58 |       0.44168 | 9.781799316 |             1 |          0.0927713 | noncoding |
| Peve_lncRNA_6606  |             26 |       0.39585 | 9.524597168 |             1 |          0.0103367 | noncoding |
| Peve_lncRNA_5668  |             75 |       0.35911 | 11.43731689 |             1 |          0.0698126 | noncoding |
| Peve_lncRNA_5991  |              0 |       0.36037 |         0.0 |            -1 |         2.3051e-06 | noncoding |
| Peve_lncRNA_7422  |             19 |       0.37422 | 5.823059082 |            -1 |          0.0161845 | noncoding |
| Peve_lncRNA_2714  |              0 |       0.34808 |         0.0 |            -1 |         1.1485e-06 | noncoding |
| Peve_lncRNA_9893  |             25 |       0.35681 | 8.841491699 |             1 |         0.00654519 | noncoding |
| Peve_lncRNA_1809  |             54 |       0.44657 | 11.63092041 |             1 |          0.0965145 | noncoding |
| Peve_lncRNA_6695  |             24 |       0.47029 | 6.494934082 |             1 |         0.00763753 | noncoding |
| Peve_lncRNA_5732  |             27 |       0.35076 | 8.631286621 |             1 |         0.00662266 | noncoding |
| Peve_lncRNA_2842  |             52 |       0.32582 | 9.034851074 |             1 |          0.0200084 | noncoding |
| Peve_lncRNA_5392  |            104 |       0.38750 | 8.385681152 |             1 |           0.301176 | noncoding |
| Peve_lncRNA_6925  |             39 |       0.43292 |  10.5824585 |             1 |          0.0383648 | noncoding |
| Peve_lncRNA_6696  |            162 |       0.31435 | 10.01654053 |             1 |            0.48732 | noncoding |
| Peve_lncRNA_3490  |             47 |       0.27944 | 8.976867676 |             1 |           0.017202 | noncoding |
| Peve_lncRNA_8354  |             27 |       0.33239 | 7.979309082 |             1 |         0.00548449 | noncoding |
| Peve_lncRNA_996   |             55 |       0.43804 | 5.869567871 |             1 |          0.0582023 | noncoding |
| Peve_lncRNA_9305  |             16 |       0.43063 | 10.91607666 |             1 |          0.0131755 | noncoding |
| Peve_lncRNA_6872  |              7 |       0.36188 | 8.498474121 |             1 |        1.86288e-05 | noncoding |
| Peve_lncRNA_512   |             78 |       0.31973 | 8.706359863 |             1 |          0.0508621 | noncoding |
| Peve_lncRNA_2511  |             59 |       0.38720 | 9.489196777 |             1 |          0.0460588 | noncoding |
| Peve_lncRNA_2565  |             43 |       0.36227 | 10.00030518 |             1 |          0.0208519 | noncoding |
| Peve_lncRNA_2574  |             68 |       0.27517 | 9.819030762 |             1 |           0.037142 | noncoding |
| Peve_lncRNA_465   |             86 |       0.39172 |  9.18170166 |             1 |           0.156744 | noncoding |
