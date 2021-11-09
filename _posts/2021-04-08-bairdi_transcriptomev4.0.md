---
layout: post
title: Using C. bairdi transcriptome v4.0
date: '2021-04-08'
categories: hematodinium
tags: cbai_transcriptomev4.0, kallisto, DESeq2, GO-MWU
---

## Recap

Last time, I described aligning all libraries to the _Hematodinium_-specific transcriptome hemat_transcriptomev1.6 and running them through my pipeline (kallisto -> DESeq2 -> GO-MWU). I made two comparisons - Ambient Day 2 vs. Elevated Day 2 and Elevated Day 0 vs. Elevated Day 2 (reminder: Day 0 samples were taken when all crabs were still at ambient temp). 

Since then, I repeated the same process, but with a _C. bairdi_ - specific transcriptome, cbai_transcriptomev4.0! I made the same two comparisons as before. Scripts are noted with a [prefix in the 40s](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/scripts). Prior to aligning libraries, I first built an index - that process is shown within the scripts.

Scripts available as follows:

[Creating BLASTx index for cbai_transcriptomev4.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/4_0_cbai4.0_indexcreation.ipynb)

[Running kallisto](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/4_1_download_libraries_run_kallisto.ipynb)

[Running DESeq2](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/4_2_kallisto_to_deseq_to_accessionIDs.Rmd)

[Obtaining GO terms](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/4_3_uniprot_to_GO.Rmd)

[Preparing one of the inputs for GO-MWU](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/4_4_eliminate_duplicates.ipynb)

[Preparing the other input for GO-MWU](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/4_5_GO-MWU_prep.Rmd)

[Running GO-MWU](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/4_6_running_GO-MWU/4_6_running_GO-MWU.R)

Again, two comparisons were made - Elevated Day 0 vs. Elevated Day 2, individual libraries only (remember, Elevated Day 0 samples were taken prior to exposure to elevated-temperature water), and Ambient Day 2 vs. Elevated Day 2, individual libraries only. All 3 groups (Elev. Day 0, Elev. Day 2, Amb. Day 2) have 3 libraries.

## Results of Analysis 3

## [DESeq2 results](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/DESeq2_output/cbai_transcriptomev4.0)

Alright, next steps on this are as follows: 
- Learn to run WGCNA
- Run WGCNA on a single crab over time, and then scale up from there, potentially running GO-MWU within each module

I also have some next steps when it comes to planning future Hematodinium-culturing projects! Generally, I'm planning to explore the logistical feasibility of working out of either Juneau or Kodiak. Chatted with Andy the other day, emailing his old colleage April Rebert at Juneau ADFG right now!
