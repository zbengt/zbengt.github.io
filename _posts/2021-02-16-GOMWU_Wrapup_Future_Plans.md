---
layout: post
title:  GO-MWU Wrapup + Future Plans
date: '2021-02-16'
categories: hematodinium
tags: overview, GO-MWU
---

## Wrapup of GO-MWU

Last post, I had run GO-MWU on my two analyses - Elevated Day 2 vs. Ambient Day 0+2 (individual libraries only), and Ambient Day 0+2+17 + Elevated Day 0 + Lowered Day 0 vs. Elevated Day 2 (individual and pooled libraries). 

I've been continuing on a bunch of different analyses and running them through my pipeline. Reminder: pipeline is kallisto matrix -> DESeq2 -> GO-MWU (with steps in between for formatting). I've ran the following pairwise comparisons all the way through both:

- Ambient Day 0 vs. Ambient Day 2 (individual libraries only)
- Ambient Day 0 vs. Ambient Day 17 (individual libraries only)
- Ambient Day 2 vs. Ambient Day 17 (individual libraries only)
- Ambient Day 2 vs. Elevated Day 2 (individual libraries only)
- Elevated Day 0 vs. Elevated Day 2 (individual libraries only)

[Here's a link to all visualizations of significant GO categories](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/GOMWU_output/cbai_transcriptomev2.0), and [here's a link to the output from GO-MWU](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/output/GO-MWU_output/cbai_transcriptomev2.0).


Findings from GO-MWU analyses:

-  Tons of significant GO categories for Ambient Day 0 or 2 vs. 17!
-  Not many significant GO categories when looking at different crabs with different temperature groups (ex: Day 2 ambient vs. Day 2 elevated). However, when we track the same crabs over time (ex: Day 0 elevated vs. Day 2 elevated), we see quite a few significant GO categories. Reminder: Day 0 elevated samples were taken when all crabs were held at ambient water temperature.
- To me, the most promising areas for further analysis are Ambient Day 0 vs. Ambient Day 17, Ambient Day 2 vs. Ambient Day 17, and Elevated Day 0 vs. Elevated Day 2

## What direction to go next

At this point, the project is splitting along two general lines - DEG analysis and correlation analysis. I'll break down the plans for each of those.

### DEG analysis:

The goal of this is to investigate the role of temperature on _Hematodinium_ and _C. bairdi_ expression. Lists of DEGs will be produced for two particularly relevant comparisons - [Ambient Day 2 vs. Elevated Day 2 (individual libraries only)](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/amb2_vs_elev2_indiv/DEGlist_wcols.txt) and [Elevated Day 0 vs. Elevated Day 2 (individual libraries only)](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/elev0_vs_elev2_indiv/DEGlist_wcols.txt). Both contrast samples taken at ambient temperatures vs. samples taken at elevated temperatures. 

We will take these DEGs and BLAST them against an NCBI database of all _Alveolata_ sequences (_Alveolata_ is the superphylum containing _Hematodinium_, as there are relatively few _Hematodinium_ sequences). From there, we can specifically investigate differentially-expressed _Hematodinium_ genes on an individual level.

### Correlation Analysis
Essentially, this tries to determine the roles of both temperature and time. This involves taking a table of all genes and their transcripts per million (TPM) from the kallisto output files, and examining the correlation.

I'm still in the really early stages of figuring this stuff out, but based on the lab meeting today, there are two possibilities for a correlation analysis. 

Option 1: [ANOVA-simultaneous component analysis (ASCA)](https://academic.oup.com/bioinformatics/article/21/13/3043/197836). Shelly's done quite a bit of work on this, including a paper that examines [the effect of both temperature and time](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07127-3), which seems extremely relevant to what I'm doing.

Option 2: [Weighted correlation network analysis (WGCNA)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559). Yaamini's the expert on this one - [here's a script she wrote](https://github.com/eimd-2019/project-EWD-transcriptomics/blob/master/analyses/WGCNA/WGCNA.md). 

Both seem like good options, but I'm not really working on this immediately. My plan is to [try out the WGCNA tutorial Yaamini linked in her script](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html), do some reading, and see which one I prefer/fits project goals better.

