---
layout: post
title: CEABiGR lncRNA - WGCNA-01
date: '2023-09-28'
categories: bioinformatics CEABiGR
tags: CEABiGR WGCNA
---

Running Ariana's version which has much more stringent filtering steps before performing the WGCNA:
* Removing genes that don't show any expression (done in previous WGCNA)
* Choosing how many 0s to allow per gene
* Choosing expression level thresholds

Current status of code is [here](https://github.com/zbengt/oyster-lnc/blob/main/code/06-WGCNA-TreeCut-distance-matrix.Rmd).

After completing these additional filtration steps, 4190 of the 4750 lncRNAs of interest were kept. Currently working through an issue at the TOM dissimilarity step, coming back to it and will post an issue tagging Ariana if I can't get it within 15 min of coming back to it.

Here is the issue...

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/WGCNA-01-tom-dissim.png?raw=true)

Clustering at this stage should represent all genes, this clearly isn't that. So figuring out where the error is. Seems like it might be a data loading issue.

UPDATE: I made some stupid mistakes, no surprise there.
Next steps discussed with Ariana:
* Go back and normalize data with DESeq2 steps
* Change the filtering (5 counts, 1 sample out of all of them)
* Re-select the soft-threshold after normalization
* Make expression histograms and do some of the outlier identification we talked about (dendrogram - module, histograms - specific lncRNAs)



