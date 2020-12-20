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

Alright, here's a link to my results from the DESeq analyses. Again, the most important is the table of genes with significantly different expressions - :

**Day 0 Ambient vs. Day 17 Ambient**

(also posted Day 0 Amb vs. Day 17 Amb results  yesterday, posting here for completion)
- [My DESeq2 script](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/Scripts/day0_v_day17_ambient_comparison.R)
-[Table of genes with significantly different expressions (adjusted pval <= 0.005)](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/0vs17_DEGlist.txt), [same table but with column headers](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/0vs17_DEGlist_wcols.txt)
- A variety of MA plots with the following conditions:
    - [All results](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/allres_MAplot.png)
    - [All results, LFC estimates shrunk using apeglm](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/allres_shrunken_MAplot.png)
    - [Results with a p-value <= 0.05](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/res05_MAplot.png)
    - [All results, with p-values <= 0.005 highlighted](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/normalizedcts_v_log2foldchange.png)
- [PCA plot](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/PCA_plot.png)
- [Dispersion estimates](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/day0_day17_ambient/dispersion_estimates.png)

**Day 0+2 Elevated vs. Day 0+2 Low**
- [My DESeq2 script](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/Scripts/day02_elev_v_day02_low_comparison.R)
-[Table of genes with significantly different expressions (adjusted pval <= 0.005)](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_low_day02/Elev_vsLow_DEGlist.txt), [same table but with column headers](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_low_day02/Elev_vsLow_DEGlist_wcols.txt)
- A variety of MA plots with the following conditions:
    - [All results](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_low_day02/allres_MAplot.png)
    - [All results, LFC estimates shrunk using apeglm](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_low_day02/allres_shrunken_MAplot.png)
    - [Results with a p-value <= 0.05](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_low_day02/res05_MAplot.png)
    - [All results, with p-values <= 0.005 highlighted](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_low_day02/normalizedcts_v_log2foldchange.png)
- [PCA plot](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_low_day02/PCA_plot.png)
- [Dispersion estimates](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_low_day02/dispersion_estimates.png)

**Day 0+2 Elevated vs Day 0+2 Ambient**
- [My DESeq2 script](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/Scripts/day02_elev_v_day02_amb_comparison.R)
-[Table of genes with significantly different expressions (adjusted pval <= 0.005)](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_amb_day02/Elev_vsAmb_DEGlist.txt), [same table but with column headers](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_amb_day02/Elev_vsAmb_DEGlist_wcols.txt)
- A variety of MA plots with the following conditions:
    - [All results](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_amb_day02/allres_MAplot.png)
    - [All results, LFC estimates shrunk using apeglm](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_amb_day02/allres_shrunken_MAplot.png)
    - [Results with a p-value <= 0.05](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_amb_day02/res05_MAplot.png)
    - [All results, with p-values <= 0.005 highlighted](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_amb_day02/normalizedcts_v_log2foldchange.png)
- [PCA plot](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_amb_day02/PCA_plot.png)
- [Dispersion estimates](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/elev_v_amb_day02/dispersion_estimates.png)

**Day 0+2 Ambient vs Day 0+2 Low**
- [My DESeq2 script](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/Scripts/day02_amb_v_day02_low_comparison.R)
-[Table of genes with significantly different expressions (adjusted pval <= 0.005)](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/amb_v_low_day02/Amb_vsLow_DEGlist.txt), [same table but with column headers](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/amb_v_low_day02/Amb_vsLow_DEGlist_wcols.txt)
- A variety of MA plots with the following conditions:
    - [All results](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/amb_v_low_day02/allres_MAplot.png)
    - [All results, LFC estimates shrunk using apeglm](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/amb_v_low_day02/allres_shrunken_MAplot.png)
    - [Results with a p-value <= 0.05](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/amb_v_low_day02/res05_MAplot.png)
    - [All results, with p-values <= 0.005 highlighted](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/amb_v_low_day02/normalizedcts_v_log2foldchange.png)
- [PCA plot](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/amb_v_low_day02/PCA_plot.png)
- [Dispersion estimates](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/amb_v_low_day02/dispersion_estimates.png)




