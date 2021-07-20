---
layout: post
title:  Day 0 Ambient vs. Day 17 Ambient
date: '2020-12-18'
categories: hematodinium
tags: hematodinium, DESeq2, Kallisto
---

**Comparison**

Yesterday, I described how I planned to start with a comparison of Day 0 and Day 17 libraries at ambient temperatures. I was able to complete that! Methods followed the same protocol outlined previously - use Kallisto to create pseudoalignments for libraries, use a script from Trinity pipeline to create matrix of counts, put matrix of counts into DESeq2 and analyze. Here are links to the following results:

- [My DESeq2 script](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/1_2_kallisto_to_deseq_to_accessionIDs.Rmd)

- [Results](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/DESeq2_output/cbai_transcriptomev2.0), which contain the following:
- Table of genes with significantly different expressions (adjusted pval <= 0.005), both with and without column headers
- A variety of MA plots with the following conditions:
    - All results
    - All results, LFC estimates shrunk using apeglm
    - Results with a p-value <= 0.05
    - All results, with p-values <= 0.005 highlighted
- Dispersion estimates

The most important of those is likely the table of genes with significantly different expressions, but before I do much with that, I'm going to do the next 3 comparisons, as follow:
- Day 0+2 Elevated vs Low (127/173/72/272/280/294 vs 151/254)
- Day 0+2 Elevated vs Ambient (127/173/72/272/280/294 vs 118/132/178/334/349/359)
- Day 0+2 Ambient vs Low (118/132/178/334/349/359 vs 151/254)

Hopefully I've got enough storage space for all those libraries!






