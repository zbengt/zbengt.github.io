---
layout: post
title: Manual Clustering
date: '2021-05-19'
categories: hematodinium
tags: cbai_transcriptomev2.0, WGCNA, clustering
---

## Moving Forward?

Good news: I made a whole bunch of progress

Bad news: Doesn't seem like WGCNA will work out well

So I figured out how to run RStudio on Mox, and then used that to run WGCNA. I ran every group of crab over time, and with all possible transcriptomes. In case you haven't been following along too closely, that means I ran WGCNA on ambient-temperature crabs (that's crabs A, B, and C), aligned to a transcriptome unfiltered by taxa (cbai_transcriptomev2.0). I then repeated for those same libraries, but aligned to a transcriptome BLASTed against a _Chionoecetes_ genome (cbai_transcriptomev4.0) and another transcriptome filtered to contain only _Alveolata_ sequences (hemat_transcriptomev1.6). I then repeated the process for elevated-temperature crabs (or crabs G, H, and I). 

All those scripts are available in [my scripts directory](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/scripts), with the notation 5#\_WGCNA\_[transcriptome of choice]_[temperature_treatment]Crabs.Rmd

Results are available within [output/WGCNA_output](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/output/WGCNA_output).

WGCNA produces modules, which cluster the genes together based on expression patterns. I graphed each module separately for each run of WGCNA. However, my initial graphs had some problems. First, they didn't show differences in expression by crab, and second, many graphs had extremely distorted scales from a few anomalous counts (one or two genes would have, say, 90,000 counts at a given time point, and nothing else was particularly visible). As a result, I graphed again, producing a total of 5 graphs per module. They are as follows (example for a hypothetical module named blue)

blue_module_TPM_graph.png: TPM counts over time

blue_module_TPM_graph_crab_colors.png: TPM counts over time, with line colors corresponding to crab IDs

blue_module_TPM_graph_crab_colors_short_yaxis.png: TPM counts over time, with line colors corresponding to crab IDs and the y-axis height limited at 100 TPM

blue_module_log_TPM_graph.png: TPM counts over time on the log scale

blue_module_log_TPM_graph_crab_colors.png: TPM counts over time on the log scale, with line colors corresponding to crab ID.

There is a LOT of data here. Nothing looks particularly illuminating at first, but I really do need to go back over this to see what (if anything) pops out. In the meantime, I've moved on to working on creating heat maps of gene expression within each individual crab. Laura posted a [very neat way of accomplishing this](https://github.com/RobertsLab/resources/discussions/1206), which I have adapted for larger numbers of crabs.

For the sake of completeness, I downloaded and aligned the libraries for crabs D and E (previously excluded from this analysis as they were uninfected) and aligned them to cbai_transcriptomev4.0. As they were uninfected, it didn't make much sense to align them to a Hematodinium-only or an unfiltered transcriptome. 

Regardless, scripts are available for [cbai_transcriptomev2.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/7_1_manual_clustering_cbaiv2.0.Rmd), [hemat_transcriptomev1.6](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/7_2_manual_clustering_hematv1.6.Rmd), and for [cbai_transcriptomev4.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/7_3_manual_clustering_cbaiv4.0.Rmd).

There is an element of subjectivity that I'm not a huge fan of when it comes to this method, as it requires me to determine how many modules are apparent in each heat map. However, it is an interesting way to track individual crab over time, and I'm intrigued by what it might reveal.

As an important note, each cluster and accompanying heatmap is labeled with a two- or three-letter string. These correspond to the subjective levels of expression present. A module that goes low expression on Day 0 -> high expression on Day 2 -> medium expression on Day 17 will be labeled cluster_LHM.txt and cluster_LHM_heatmap.png. Elevated-temperature crab only have 2 timepoints, and therefore only 2 clusters (LH and HL).


[NOTE: These labels were later changed to only describe overall expression patterns - see next blog post for details]

I ran the full script for hemat_transcriptomev1.6 and cbai_transcriptomev4.0. Output from my manual clustering is available [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/output/manual_clustering), with each transcriptome in a separate directory, and then each crab in a separate sub-directory. However, when I try this with cbai_transcriptomev2.0, Rstudio just refuses to run it - not enough memory on my local machine. 

To solve this, it's back to RStudio on Mox. But first, I need to figure out how to install the `pheatmap` package - shouldn't be too tricky for, say, Sam, but I did spend a few hours today trying and failing to figure it out. Hopefully I can get some help and keep moving along! 

Next steps:
- Examine the output of my manual clustering
    - Compare modules to see what the main expression patterns were
    - Compare genes between crabs to see if there's overlap
- Install `pheatmap` and continue by examining cbai_transcriptomev2.0
- Try repeating, but using relative gene counts (to Day 0) rather than TPM
- Go through WGCNA output in detail (and get a better idea of what I'm looking for)


