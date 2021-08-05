---
layout: post
title: Crab Plasticity
date: '2021-08-04'
categories: hematodinium
tags: PCAs, plasticity
---

In our previous post, we created PCAs for each transcriptome, with each point as one library. This raised some interesting questions, notably that we didn't necessarily see libraries from the same crab clustered together (though for some crab, they did clearly cluster). This could be an indication of crab plasticity, as they adjust to both tank effects and the temperature treatment. 

We initially considered looking at plasticity qualitatively, but figured it'd be better to go the extra mile and take a quantitative approach. Therefore, we modified our DESeq2 analysis scripts for our three transcriptomes -  [cbai_v2.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/1_2_kallisto_to_deseq_to_accessionIDs.Rmd), [cbai_v4.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/4_2_kallisto_to_deseq_to_accessionIDs.Rmd), and [hemat_v1.6](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/3_2_kallisto_to_deseq_to_accessionIDs.Rmd). Yet another reminder - cbai_v2.0 is unfiltered, cbai_v4.0 is presumably _Chionoecetes_ genes only, and hemat_v1.6 is presumably _Hematodinium_ genes only.

These modifications created a PCA for each, extracted the coordinates for PC1 and PC2 for each library, used that to get the centroid for each crab, and then calculated the mean distance between each library and the centroid for that particular crab. That gave us a number to describe how closely each crab's libraries clustered for each transcriptome. Honestly, it was some pretty neat work and I'm proud of it! 

The raw data is available for [cbai_v2.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/all_indiv_libraries/PCA_plasticity.txt), [cbai_v4.0](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev4.0/all_indiv_libraries/PCA_plasticity.txt), and [hemat_v1.6](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/hemat_transcriptomev1.6/all_indiv_libraries/PCA_plasticity.txt). 

But in case you want something a bit simpler, here's the rankings for the crabs for each transcriptome! 1 = least plastic (lowest mean distance to centroid), 9 = most plastic (highest mean distance to centroid).

Additional reminder: We have 3 libraries for crabs A-F (days 0, 2, and 17) and 2 for G-I (days 0 and 2)

I also added a few columns at the end:

- Total Plasticity Points: The sum of each crab's rankings. This is a bit redundant, since cbai_v4.0 and hemat_v1.6 are both derived from cbai_v2.0 (though there's plenty of cbai_v2.0 genes that aren't in either), but is still an interesting measure of how much the libraries for each crab differ for all transcriptomes

- Total Plasticity Ranking: Same as Total Plasticity Points, but ranked. Again, 1 = least plastic, 9 = most plastic.

| Crab | Treatment | cbai_v2.0 | cbai_v4.0 | hemat_v1.6 | Total_Plasticity_Points | Total_Plasticity_Ranking |
|------|-----------|-----------|-----------|------------|-------------------------|--------------------------|
| A    | Ambient   | 6         | 9         | 7          | 22                      | 7                        |
| B    | Ambient   | 8         | 8         | 9          | 25                      | 9                        |
| C    | Ambient   | 1         | 2         | 6          | 9                       | T2                       |
| D    | Lowered   | 4         | 7         | 5          | 16                      | T5                       |
| E    | Lowered   | 2         | 1         | 1          | 4                       | 1                        |
| F    | Lowered   | 7         | 5         | 4          | 16                      | T5                       |
| G    | Elevated  | 9         | 6         | 8          | 23                      | 8                        |
| H    | Elevated  | 3         | 3         | 3          | 9                       | T2                       |
| I    | Elevated  | 5         | 4         | 2          | 11                      | 4                        |

A few things of note here:

First, there doesn't seem to be a clear, explicable relationship between temperature treatment and plasticity. Ambient-temperature crabs have extremely high plasticity points for all 3 transcriptomes. This indicates that we could be looking at either natural variability or tank effects here

Second, crabs with high plasticity for cbai_v4.0 also generally seem to have a high plasticity for hemat_v1.6. This is pretty interesting! Could indicate expression in both _Hematodinium_ and _C. bairdi_ is shifting simultaneously as each tries to find weaknesses in the defense of the other. Definitely worth examining more!

Third, both uninfected crabs (D and F) have a pretty middle-of-the-pack plasticity for all samples. Note that these crab might not actually be uninfected (judging by our qPCR results).

