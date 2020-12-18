---
layout: post
title:  Day 0 Ambient vs. Day 17 Ambient
date: '2020-12-18'
categories: hematodinium
tags: hematodinium, DESeq2, Kallisto
---

**Comparison**

Yesterday, I described how I planned to start with a comparison of Day 0 and Day 17 libraries at ambient temperatures. I was able to complete that! Methods followed the same protocol outlined previously - use Kallisto to create pseudoalignments for libraries, use a script from Trinity pipeline to create matrix of counts, put matrix of counts into DESeq2 and analyze. Here are links to the following results:

- [My DESeq2 script](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/Scripts/day0_v_day17_ambient_comparison.R)
-[Table of genes with significantly different expressions (adjusted pval <= 0.005)](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/0vs17_DEGlist.txt), [same table but with column headers](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/0vs17_DEGlist_wcols.txt)
- A variety of MA plots with the following conditions:
    - [All results](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/allres_MAplot.png)
    - [All results, LFC estimates shrunk using apeglm](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/allres_shrunken_MAplot.png)
    - [Results with a p-value <= 0.05](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/res05_MAplot.png)
    - [All results, with p-values <= 0.005 highlighted](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/normalizedcts_v_log2foldchange.png)
- [PCA plot](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/PCA_plot.png)
- [Dispersion estimates](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/dispersion_estimates.png)

The most important of those is likely the table of genes with significantly different expressions, but before I do much with that, I'm going to do the next 3 comparisons, as follow:
- Day 0+2 Elevated vs Low (127/173/72/272/280/294 vs 151/254)
- Day 0+2 Elevated vs Ambient (127/173/72/272/280/294 vs 118/132/178/334/349/359)
- Day 0+2 Ambient vs Low (118/132/178/334/349/359 vs 151/254)

Hopefully I've got enough storage space for all those libraries!






