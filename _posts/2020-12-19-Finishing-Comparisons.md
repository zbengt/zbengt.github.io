---
layout: post
title:  Lists of Differentially-Expressed Genes (plus graphs)
date: '2020-12-19'
categories: hematodinium
tags: hematodinium, DESeq2, Kallisto
---

**Outline**

I put a whole bunch of libraries through the same protocol as before - use Kallisto to create pseudoalignments for libraries, use a Trinity pipeline script to merge Kallisto counts into a matrix, analyze matrix with DESeq2, produce list of significantly different genes (padj <= 0.005)

At this point, I now have a list of significantly-different genes for each difference in condition

Next steps: 
- Create newline-separated text file containing significantly-different genes for each comparison. Is input for next step.
- Use [this script](https://github.com/RobertsLab/code/blob/master/script-box/uniprot2go.sh) to get GO terms for significantly-different genes for each comparison.
- Use GO assignments to perform enrichment analysis 

As before, the comparisons are the following:
- Day 0 Ambient vs. Day 17 Ambient (118, 132, 178 vs 463, 481, 485)
- Day 0+2 Elevated vs Low (127/173/72/272/280/294 vs 151/254)
- Day 0+2 Elevated vs Ambient (127/173/72/272/280/294 vs 118/132/178/334/349/359)
- Day 0+2 Ambient vs Low (118/132/178/334/349/359 vs 151/254)

Alright, here's a link to my results from the DESeq analyses. Again, the most important is the table of genes with significantly different expressions.

[Here's the script used to generate all results](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/1_2_kallisto_to_deseq_to_accessionIDs.Rmd)

- [Results](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/DESeq2_output/cbai_transcriptomev2.0), which contain the following:
- Table of genes with significantly different expressions (adjusted pval <= 0.005), both with and without column headers
- A variety of MA plots with the following conditions:
    - All results
    - All results, LFC estimates shrunk using apeglm
    - Results with a p-value <= 0.05
    - All results, with p-values <= 0.005 highlighted
- Dispersion estimates



