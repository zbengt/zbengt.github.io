---
layout: post
title: Manual Modules 3
date: '2021-07-12'
categories: hematodinium
tags: manual clustering
---

UPDATE: Changed the bar for gene expression by 2021-08-09. Protocol was the same otherwise For more info, refer to the lab notebook post for Aug. 9th, 2021

This hasn't exactly been the most productive few weeks I've had, but for some extremely happy reasons! Matthew and I took a trip up to Alaska and got engaged!! As a result, I haven't exactly been too focused on crabs for the past week or two, but you know what - I think it's worth it.

Alright, now we're back in Seattle, so now it's back to business. A quick and broad recap:

Goal (primary): Analyze gene expression in _Hematodinium_-infected _Chionoecetes_ 

At this point, rather than examining differentially-expressed genes, I'm clustering genes into modules and then examining the composition of each module over time/temperature treatment There are two ways I'm trying to accomplish this. 

- WGCNA. This takes a complete set of expression data and clusters genes into modules based on similar expression patterns. It then examines modules to see if any show a significant correlation with either time or temperature treatment (or, if including the few uninfected crabs we have, infection status)

- Manual clustering (the focus of the rest of this post). "Manual" is a bit of a misnomer - all I'm doing is setting the cut height and naming the modules. In this method I take a set of expression data from **one crab**. This means 3 libraries if examining ambient- or lowered-temp crab, and just 2 if examining elevated-temp crab. I then use a script to cluster genes based on similar expression patterns. A cut height is then chosen for all crab (well, for each treatment group - crab with 2 timepoints need a higher cut height) and applied to the gene cluster to cut it into modules. I chose a cut height that produced ~5-7 modules for most crab. Modules are then assigned one of the following names

    - 3 timepoints (ambient- and lowered- temperature crab))
        - High to low (HTL): Expression decreases over time (regardless of whether the decrease took place on Day 2 or Day 17)

        - Low to high (LTH): Expression increases over time (regardless of whether the increase took place on Day 2 or Day 17)

        - Low High Low (LHL): Expression increases from Day 0 to Day 2, and then drops back on Day 17

        - High Low High (HLH): Expression drops from Day 0 to Day 2 and then increases on Day 17

        - Mixed (MIX): Expression within the module follows no clear pattern

    - 2 timepoints (elevated-temperature crab):
        - Low Low (LL): expression stays low

        - High High (HH): expression stays high

        - Low High (LH): expression goes from low to high

        - High Low (HL): expression goes from high to low

        - Mixed (MIX): no clear pattern of expression within the module

Note that this is just a bit different from my previous naming scheme, as described in my post titled Manual Modules 2 (I added a MIX category)!

For many crabs, multiple modules had the same designation. For these, a number was added to the end of the label if > 1. For instance, let's say Crab B had two LHL modules - these would be labeled cluster_LHL and cluster_LHL2.

However, we don't really care how many modules we have with a designation - we care about how many genes we have with a designation! So as a result, I merged all modules with the same designation for each crab/transcriptome combo. Here are scripts for the following:

- [Unfiltered transcriptome](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/7_4_merging_manual_clusters_cbaiv2.0.ipynb)

- [Transcriptome filtered to only include presumed C. bairdi genes](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/7_6_merging_manual_clusters_cbaiv4.0.ipynb)

- [Transcriptome filtered to only include presumed _Hematodinium_ genes](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/7_5_merging_manual_clusters_hematv1.6.ipynb)

I then produced tables counting the number and percentage of genes belonging to each category for each crab/transcriptome combination. To avoid posting 6 lengthy tables below, I'll just say to follow the link - two tables (one for raw counts, the other for percentage) is available at the bottom of each script for the three transcriptomes above.

For all three transcriptomes, the percentages in each category varied pretty wildly. This held true both between temperature treatments, and within each treatment. This indicates one of three things. 
- 1: Simply, my manual clustering may not have been stringent enough in its filtering. I was quite generous with what got included (>= 5 counts for that crab), and should likely bump that up a bit. Suggestion was >= 30 counts in at least 1 sample (counting ALL crabs, not just the selected crab)

- 2: There's just a ton of natural variation in crab expression patterns

- 3: Another factor - perhaps the level of _Hematodinium_ infection - could be responsible! 

We really want to figure out whether infection level plays a role in any of this, so top priorities are as follows:

1. Track down data on infection level of crabs over time, and plot 

2. For each of our three transcriptomes, make a PCA with one point for each library. See if points cluster together by day, temperature, infection level (from qPCR data), and crab. This will be done using DESeq2

3. Repeat manual clustering, but with a higher bar for filtering out low-count genes (as described above)

We also have a few lower-tier priorities:

1. Align crabs D and F to hemat_transcriptomev1.6 and cbai_transcriptomev2.0. This should just check whether our alignment works well - since these are uninfected crab, we shouldn't get any alignments to hemat_v1.6, and the alignment to cbai_v2.0 (which is unfiltered) should give an idea of the accuracy of our filtering for cbai_v4.0 (which is filtered to exclude all non-bairdi sequences).

2. Write a blog post on the bittercrab website!

3. For WGCNA, try reducing minimum module size, particularly for hemat_transcriptomev1.6