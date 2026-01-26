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

| ID                | peptide_length | Fickett_score |      pI | ORF_integrity | coding_probability | label     |
|-----------|----------:|----------:|----------:|----------:|----------:|-----------|
| Apul_lncRNA_11205 |            288 |       0.43082 |  9.4564 |            -1 |           0.460064 | noncoding |
| Apul_lncRNA_12375 |             71 |       0.33092 |  4.6010 |             1 |           0.103221 | noncoding |
| Apul_lncRNA_12954 |             22 |       0.38628 |  9.8625 |             1 |          0.0090556 | noncoding |
| Apul_lncRNA_13651 |             31 |       0.34779 | 11.2637 |             1 |          0.0158401 | noncoding |
| Apul_lncRNA_13665 |             67 |       0.37510 | 10.3766 |             1 |          0.0649536 | noncoding |
| Apul_lncRNA_14079 |            113 |       0.49597 |  6.6369 |            -1 |          0.0654743 | noncoding |
| Apul_lncRNA_14131 |             47 |       0.32934 |  7.8371 |             1 |          0.0133536 | noncoding |
| Apul_lncRNA_14842 |             46 |       0.35722 | 10.1774 |             1 |          0.0238125 | noncoding |
| Apul_lncRNA_16008 |             41 |       0.31522 |  4.4278 |             1 |          0.0225572 | noncoding |
| …                 |              … |             … |       … |             … |                  … | …         |
| Ptuh_lncRNA_11883 |            113 |       0.32918 | 10.2027 |             1 |           0.186662 | noncoding |
| Ptuh_lncRNA_7018  |             67 |       0.33099 |  9.9952 |             1 |          0.0430005 | noncoding |

### Apul–Peve (CPC2 results)

| ID                | peptide_length | Fickett_score |      pI | ORF_integrity | coding_probability | label     |
|-----------|----------:|----------:|----------:|----------:|----------:|-----------|
| Apul_lncRNA_10150 |             35 |       0.36883 |  4.7492 |             1 |          0.0152901 | noncoding |
| Apul_lncRNA_11023 |             52 |       0.33416 |  9.3023 |             1 |          0.0219726 | noncoding |
| Apul_lncRNA_12170 |             52 |       0.35497 |  7.7951 |             1 |          0.0190474 | noncoding |
| Apul_lncRNA_12365 |            100 |       0.39199 |  6.6666 |             1 |           0.327169 | noncoding |
| Apul_lncRNA_12639 |             75 |       0.34107 |  6.8889 |             1 |          0.0561701 | noncoding |
| Apul_lncRNA_13011 |             93 |       0.41674 |  6.2690 |            -1 |          0.0429362 | noncoding |
| Apul_lncRNA_13934 |             19 |       0.32935 | 10.2890 |             1 |         0.00864502 | noncoding |
| Apul_lncRNA_14577 |            103 |       0.33642 | 10.9606 |             1 |           0.140515 | noncoding |
| Apul_lncRNA_14938 |             33 |       0.43317 | 10.8343 |            -1 |           0.114709 | noncoding |
| …                 |              … |             … |       … |             … |                  … | …         |
| Peve_lncRNA_2565  |             43 |       0.36227 | 10.0003 |             1 |          0.0208519 | noncoding |
| Peve_lncRNA_2574  |             68 |       0.27517 |  9.8190 |             1 |           0.037142 | noncoding |
| Peve_lncRNA_465   |             86 |       0.39172 |  9.1817 |             1 |           0.156744 | noncoding |
