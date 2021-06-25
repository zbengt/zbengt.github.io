---
layout: post
title:  Workflow Standardization
date: '2021-01-20'
categories: hematodinium
tags: deseq2, GO-MWU
---

**Research Update - _Hematodinium_ Investigation**

Yesterday's notebook update focused on my plans to investigate future research questions. This one focuses on my current project - analyzing differential expression in _Hematodinium_ in Grace's samples. 

I'm roughly at the same stage as before - I've got all my differentialy-expressed GO terms, but realized a few things. The most important of those: for analysis with GO-MWU/GOseq/some other gene enrichment analysis, I need to provide **all** my genes, not just the differentially-expressed ones. Because of that, I had to start over on my R script that produces a newline-separated file of UniProt accessions. 

When going through my R script, I realized I've learned a lot in the last few months. The script was badly written, and I had to either sub out specific parts each time, or have differnt scripts for different files. Because of that, I decided to rewrite my scripts, and turned the whole DESeq2 analysis into a single function. I did the same for the process that creates newline-separated files. I now have two important R scripts - [one containing all my functions](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/hematodinium_analysis_functions.R), and another [that calls them](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/1_2_kallisto_to_deseq_to_accessionIDs.Rmd). I also create a Venn diagram in the latter R script, since I didn't think such a specific task was worth turning into a function.

Anyway, that was most of what I was working on recently. I now have newline-separated files of my DEG and all genes, and my next steps are as follows (I was incorrect in the order for next steps in my 12/21 post):

- Get input for GO-MWU (a two-column table with genes + GO terms, and a two-column table of genes + unadjusted p-value0)
- Perform enrichment analysis using GO-MWU (or if I can't, use topGO or something else)

The process of getting GO terms for all genes (rather than just DEGs) might take a long, long time. But hey, nothing I can do about that except start on it soon!



