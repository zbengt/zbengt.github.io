---
layout: post
title: Manual Clustering 2
date: '2021-06-09'
categories: hematodinium
tags: manual clustering
---

## Manual Clustering, Continued

Alright, some updates from last time:

- Rather than setting the number of modules independently for each individual crab/transcriptome, I specified a single cut height, which was used on the dendrogram for each crab to separate into modules. For crabs A-F (where three time points are present), this cut height was 1.8, and was chosen because it generally separated crabs into 5-8 modules, all of which appeared to contain single expression patterns. Crabs G, H, and I were part of the elevated-temperature treatment group, which only had 2 time points, and so that cut height just wasn't adequate. As a result, for these crabs, I bumped up the cut height to 10. These changes can be seen in scripts 71-73. All are the same except for initial inputs and module naming but [here's one as an example](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/7_1_manual_clustering_cbaiv2.0.Rmd)

- I changed my name scheme for the manual clustering method. Previously, each timepoint was described with HML (high, medium, or low), so a module with three timepoints could be MMH (if expression was medium on Day 0, medium on Day 2, and high on Day 17). However, we decided to get a bit more general. Now, expression is grouped into overall expression patterns. There are four possibilities - high to low (HTL), low to high (LTH), high-low-high (HLH), and low-high-low (LHL). You'll notice this groups together, say, modules that previously would have been called HLL and HHL together as HTL. As a result, there's some duplication in our module names, so if there are two LTH modules, the first will be called LTH and the second is LTH2.

- I was able to run this manual clustering method for cbai_transcriptomev2.0 on Mox after installing the `pheatmap`  package! It required me to build a new Singularity container - description of how [here](https://github.com/RobertsLab/resources/discussions/1218) and [an example here](https://github.com/RobertsLab/project-oyster-oa/blob/master/code/Haws/04-methylKit-RStudio.sh)

Beyond that, not too much progress was made on the manual clustering (most of my time was spent just working on stuff for classes, although I found some interesting stuff - coming to you in a blog post soon!) Next steps on this are still fairly similar to what they used to be:

- Examine output of my manual clustering
    - Compare modules to see what the main expression patterns were
    - Compare genes between crabs to see if there's overlap
- Try using the methods from a paper (link broken, here's the DOI: https://doi.org/10.1101/2021.04.14.439692) - GitHub repo [here](https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries) - to try to do some visualization and analyze expression with GO terms
- Again, go back through WGCNA output and see if there's anything good to glean from there

