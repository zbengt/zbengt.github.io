---
layout: post
title: Manually Clustering Immune Genes
date: '2021-08-23'
categories: hematodinium
tags: immune genes, clustering
---

Previously, we created files of transcripts that matched GO terms related to immune response. However, to make our steps more logically consistent, we recently changed those files to only include transcripts that matched the specific GO term (GO:0006955) for Immune Response. We created files for each of our three transcriptomes - cbai_transcriptomev2.0 (unfiltered), cbai_transcriptomev4.0 (bairdi only), and hemat_transcriptomev1.6 (hematodinium only). That was done using [this script](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/8_1_examining_immune_genes.Rmd)

Once we had those files, we pulled both their raw counts and TPM counts for all immune-related transcripts for the two cbai transcriptomes. The hemat transcriptome was skipped, as there were only 5 transcripts that fit the immune response GO term. The raw counts were used to create PCAs of immune gene expression by each of our five possible variables - crab, day, hematodinium level, infeciton level, and temperature. Those PCAs are available for [cbai_v2.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/DESeq2_output/cbai_transcriptomev2.0/immune_genes_all_libs) and [cbai_v4.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/DESeq2_output/cbai_transcriptomev4.0/immune_genes_all_libs). Both PCAs were created with [this script](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/8_2_immune_PCAs.Rmd)

Disappointingly, when the PCAs were examined, there seems to be absolutely no clustering around any of our five variables.

Following this, I took the TPM counts and put them through the manual clustering script we previously ran on all transcripts for each transcriptome. Again, this was only done for cbai_transcriptomev2.0 and cbai_transcriptomev4.0. This involved clustering transcripts based on expression patterns, and then manually naming each cluster based on its overall pattern. The following notation was used - same as previous manual clustering work done on this project:

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

This manual clustering was done using [this script](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/8_3_manual_clustering_cbaiv2.0_immune_genes.Rmd) for cbai_v2.0, and [this script](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/8_4_manual_clustering_cbaiv4.0_immune_genes.Rmd) for cbai_v4.0

Often, there was overlap in our clusters (i.e. Crab A may have had several HTL clusters). When this occurred, a number was added to the end of the cluster name (e.g. HTL2). However, we wanted to produce tables of the raw numbers, and thus had to merge clusters with the same base name. This was done for both [cbai_v2.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/8_5_merging_immune_manual_clusters_cbaiv2.0.ipynb) and [cbai_v4.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/8_6_merging_immune_manual_clusters_cbaiv4.0.ipynb). Following this merging, we took the counts for each cluster name and turned them into a table for each transcriptome. This was performed [using this script](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/8_7_merged_immune_cluster_counts.Rmd)

Cbai_transcriptomev2.0:

- Raw counts [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/cbai_transcriptomev2.0/immune_genes/merged_modules_counts_table.csv)

- Percentages [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/cbai_transcriptomev2.0/immune_genes/merged_modules_percentages_table.csv)

Cbai_transcriptomev4.0:

- Raw counts [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/cbai_transcriptomev4.0/immune_genes/merged_modules_counts_table.csv)

- Percentages [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/cbai_transcriptomev4.0/immune_genes/merged_modules_percentages_table.csv)

Importantly, we should note that all of these counts are extremely small - in the single to low double digits. Therefore, there really might not be much beyond random noise here. With that said, here's some tables of the average percentages for each category for each of our five factors.

Note that the numbers shown are the average of the AVERAGES, not of the raw counts. For instance, if one crab had 50/100 HTL genes and the only other in its group had 0/10 HTL genes, the table below would show 25% HTL 

For variables (like infection status) where some crabs use the two-timepoints notation (HH, HL, etc) and others use the three-timepoints notation (HTL, LTH, etc), averages are from all crabs in the groups with the notation of choice. Ex: if we have 3 uninfected crabs of 10 total, and 1 is from our elevated-temp treatment, the values for the uninfected row in columns HH/HL/etc. will consist solely of the percentages in that single crab in the elevated-temp treatment group. We will count MIX as a two-timepoints column, since none of our three-timepoints crab had any mixed modules

Also if you're looking for a refresher on which crabs map to which variables, [head over here!](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/data/sample_ids.csv). Infection status and hematodinium level are treated as constants for each crab.

## Temperature Treatment

Note that this doesn't *really* map to temperature treatment, as on Day 0, all crabs were at ambient temperature. But it does map to treatment group, which is what we care more about anyways.

#### cbai_transcriptomev2.0

|              | HLH  | HTL   | LHL   | LTH   | HH | HL    | LH   | LL    | MIX   |
|--------------|------|-------|-------|-------|----|-------|------|-------|-------|
| Amb (n = 3)  | 5.6% | 28.5% | 28.8% | 37.1% | NA | NA    | NA   | NA    | 0%    |
| Low (n = 3)  | 0%   | 67.2% | 24.6% | 19.9% | NA | NA    | NA   | NA    | 0%    |
| Elev(n = 3)  | NA   | NA    | NA    | NA    | 0% | 40.4% | 7.9% | 11.1% | 40.5% |

#### cbai_transcriptomev4.0

|              | HLH  | HTL   | LHL   | LTH   | HH | HL    | LH    | LL | MIX   |
|--------------|------|-------|-------|-------|----|-------|-------|----|-------|
| Amb (n = 3)  | 7.4% | 28.5% | 32.6% | 31.5% | NA | NA    | NA    | NA | NA    |
| Low (n = 3)  | 6.3% | 57.6% | 26.4% | 9.7%  | NA | NA    | NA    | NA | NA    |
| Elev (n = 3) | NA   | NA    | NA    | NA    | 0% | 51.2% | 17.5% | 0% | 31.3% |

## Infection status


#### cbai_transcriptomev2.0

|                    | HLH | HTL   | LHL   | LTH   | HH | HL    | LH   | LL    | MIX   |
|--------------------|-----|-------|-------|-------|----|-------|------|-------|-------|
| Infected (n = 7)   | 0%  | 39.7% | 26.4% | 29.7% | 0% | 40.5% | 7.9% | 11.1% | 40.5% |
| Uninfected (n = 2) | 0%  | 64.3% | 27.4% | 8.3%  | NA | NA    | NA   | NA    | NA    |

#### cbai_transcriptomev4.0

|                     | HLH  | HTL   | LHL   | LTH   | HH | HL    | LH    | LL | MIX   |
|---------------------|------|-------|-------|-------|----|-------|-------|----|-------|
| Infected (n = 7)    | 5.6% | 40.2% | 27.6% | 26.7% | 0% | 51.2% | 17.5% | 0% | 31.3% |
| Uninfected (n = 2)  | 9.6% | 49.0% | 33.3% | 8.4%  | NA | NA    | NA    | NA | NA    |


## Hematodinium infection level

#### cbai_transcriptomev2.0

|              | HLH | HTL   | LHL   | LTH   | HH | HL    | LH    | LL    | MIX   |
|--------------|-----|-------|-------|-------|----|-------|-------|-------|-------|
| High (n = 5) | 0%  | 39.6% | 30.8% | 24.0% | 0% | 36.6% | 11.9% | 16.7% | 40.5% |
| Low (n = 4)  | 0%  | 56.1% | 22.7% | 21.2% | 0% | 48.3% | 51.7% | 0%    | 0%    |

#### cbai_transcriptomev4.0

|              | HLH  | HTL   | LHL   | LTH   | HH | HL    | LH    | LL | MIX   |
|--------------|------|-------|-------|-------|----|-------|-------|----|-------|
| High (n = 5) | 7.4% | 20.4% | 35.2% | 37.0% | 0% | 53.8% | 26.3% | 0% | 20.0% |
| Low (n = 4)  | 6.3% | 65.8% | 23.8% | 4.2%  | 0% | 46.2% | 0%    | 0% | 53.8% |

### Day

Since we're clustering based on expression over time, this really can't be done

### Crab

THis is already accomplished in the original tables linked above!

First impressions: wow it is tricky to tease apart effects here. Definitely some differences, but there's so much overlap in our variables that it's hard to determine what causes what. For instance, both our uninfected crab are also in the lowered-temp treatment group. Is the large number of HTL-clustering transcripts due to the effect of temperature or the effect of infection? The lack of repetition impacts this further. All our two-day cluster numbers for the effect of hematodinium infection level come from a single crab (Crab G). Is this an anomalous crab or is there something real happening?

Second, it's also tricky to see what's important with these extremely low numbers of genes. I'll talk to people in the lab tomorrow (first real in-person day, we're all going in to clean up the place!) to hopefully get some insights. I should also recreate these tables with the clusterings made earlier of all genes. Previously, only the effect of crab was summarized in this way - something interesting might pop out with gene numbers that are in the thousands rather than the single digits!