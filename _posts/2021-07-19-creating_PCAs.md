---
layout: post
title: Creating PCAs
date: '2021-07-19'
categories: hematodinium
tags: PCAs, kallisto, DESeq2
---

Alright, as of last time, our general goal was to do the following (along with a few other priorities):

1. Track down data on infection level of crabs over time, and plot 

2. For each of our three transcriptomes, make a PCA with one point for each library. See if points cluster together by day, temperature, infection level (from qPCR data), and crab. This will be done using DESeq2

3. Align crabs D and F to hemat_transcriptomev1.6 and cbai_transcriptomev2.0. This should just check whether our alignment works well - since these are uninfected crab, we shouldn't get any alignments to hemat_v1.6, and the alignment to cbai_v2.0 (which is unfiltered) should give an idea of the accuracy of our filtering for cbai_v4.0 (which is filtered to exclude all non-bairdi sequences).

4. Write a blog post on the bittercrab website

We can check off number 4 - take a peek [here!](https://bittercrab.wordpress.com/2021/07/13/the-value-of-crab/) The rest of the entry will focus on a combination of Numbers 1-3

First, we were able to track down the qPCR data (we think)! It appears that each samples were only taken at a single timepoint. Importantly, that timepoint is unknown as of yet - hopefully they were taken on Day 0, so that there isn't any experimental effect! Anyway, that's all available [here](https://github.com/crmateo/crab-capstone/blob/main/data/hematqpcr_crabRNA.csv). A quick note on decoding: relevant columns are Uniq_ID and sq_all.runs_mean. Uniq_ID is UniqueSampleNumber_GenewizID_Day. We've been using the Genewiz ID (those 3-digit numbers) to ID our samples, so to match each up, just look for that. Ex: our sample 123 might be, say, their 9999_123_12.

The sheet's notation of Day counts Day 0 as the day of collection, which was nine days prior to the start of the experiment (our Day 0). So essentially, our Days 0, 2, and 17 correspond to their 9, 12, and 26. You might think this is confusing, but luckily [differences in date notation have never caused issues before](https://www.si.com/extra-mustard/2013/12/30/the-extra-mustard-trivia-hour-when-a-calendar-defeated-russia-in-the-1908-olympics)

Anyway, we looked at the qPCR data and noticed a few things. First, some of the crabs that were counted as uninfected within the experiment (Crabs D and F) were, in fact, infected according to the qPCR data. And second, generally, there were 2 groups of crab - the lightly-infected (below 5000 - units presumed to be ng/ul) and heavily-infected (above 100,000). As a result, we added hematodinium infection level as a binary variable - Low or High. Here's a table of this:

| Crab | Treatment | Assigned_Infection_Status | qPCR_Infection_Level | Infection_Level_Assigned |
|------|-----------|---------------------------|----------------------|--------------------------|
| A    | Ambient   | Infected                  | 1119                 | Low                      |
| B    | Ambient   | Infected                  | 105477               | High                     |
| C    | Ambient   | Infected                  | 333000               | High                     |
| D    | Lowered   | Uninfected                | 380502               | High                     |
| E    | Lowered   | Infected                  | 1550                 | Low                      |
| F    | Lowered   | Uninfected                | 3069                 | Low                      |
| G    | Elevated  | Infected                  | 4515                 | Low                      |
| H    | Elevated  | Infected                  | 119500               | High                     |
| I    | Elevated  | Infected                  | 211000               | High                     |


Alright, so we've got a new variable in play here! Now, our next goal was to see the effect of each of our variables - crab, day, assigned infection status, qPCR infection level, and temperature. However, first we had to fix a few things. Our previous assumption was that we didn't need to align libraries for crabs D and F to cbai_transcriptomev2.0 (our unfiltered transcriptome) or hemat_transcriptomev1.6 (our hematodinium transcriptome). Instead, since those crabs were uninfected, we only had to align them to cbai_transcriptomev4.0 (our transcriptome filtered to be just bairdi sequences). However, as shown above, that...might not be true. So instead, we moved on to Step 3. This was a bit complicated and took a while, but we were able to figure it out (had to run a Trinity script on Roadrunner). So now that we had our kallisto libraries, we moved on to...

### Producing PCAs!

We used DESeq2 to produce five PCAs for each of our three transcriptomes. Each of the PCA differs only in appearance, highlighting a different variable

Reminder: our variables are
- Crab (A, B, C...)
- Day (0, 2, 17)
- Infection Status (Infected/Uninfected)
- qPCR infection level (Low/High)
- Temperature at time of sampling (Low/Ambient/Elevated)

So here are all those PCAs!

#### cbai_transcriptomev2.0
- [Crab](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/all_indiv_libraries/crab_as_var/PCA_plot.png)
- [Day](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/all_indiv_libraries/day_as_var/PCA_plot.png)
- [Infection Status](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/all_indiv_libraries/infection_as_var/PCA_plot.png)
- [Infection Level](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/all_indiv_libraries/hemat_level_as_var/PCA_plot.png)
- [Temperature](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/all_indiv_libraries/temp_as_var/PCA_plot.png)

#### cbai_transcriptomev4.0
- [Crab](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev4.0/all_indiv_libraries/crab_as_var/PCA_plot.png)
- [Day](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/DESeq2_output/cbai_transcriptomev4.0/all_indiv_libraries/day_as_var)
- [Infection Status](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev4.0/all_indiv_libraries/infection_as_var/PCA_plot.png)
- [Infection Level](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev4.0/all_indiv_libraries/hemat_level_as_var/PCA_plot.png)
- [Temperature](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev4.0/all_indiv_libraries/temp_as_var/PCA_plot.png)

#### hemat_transcriptomev1.6
- [Crab](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/DESeq2_output/hemat_transcriptomev1.6/all_indiv_libraries/crab_as_var)
- [Day](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/hemat_transcriptomev1.6/all_indiv_libraries/day_as_var/PCA_plot.png)
- [Infection Status](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/hemat_transcriptomev1.6/all_indiv_libraries/infection_as_var/PCA_plot.png)
- [Infection Level](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/hemat_transcriptomev1.6/all_indiv_libraries/hemat_level_as_var/PCA_plot.png)
- [Temperature](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/hemat_transcriptomev1.6/all_indiv_libraries/temp_as_var/PCA_plot.png)

Of these, most don't really show any clear correlations going on. However, there are a few ones of interest. 

hemat_transcriptomev1.6 seems to show some kind of correlation between infection status and temperature. However, the two are pretty correlated (all uninfected crab were part of the low-temperature treatment group). And if those 2 crabs (D and F) were, in fact, uninfected, then this certainly wouldn't be surprising at all.

Anyways, we have a few possibilities to look at when analyzing these PCAs.

First, there do seem to be clear clusters. They just don't correspond to our variables. As a result, it's possible that some batch effect is present. To check this, we'll look at the days the libraries were prepped and see if that corresponds to anything. (UPDATE 2021-08-04: All libraries were sent in on the same day)

It's also possible that another variable (such as carapace width) could have an influence here. It would be worth re-running the PCAs, just to eliminate that possibility.

The crab PCAs for cbai_v2.0 and cbai_v4.0 seem to show that some crab samples are tightly clustered, while others are all over the place. This indicates that plasticity may vary between crabs. It may be worth grouping more and less plastic crabs - either qualitatively, or based on the mean distance from the centroid of that crab.

We also have some other tasks that aren't PCA-related that should still be worked on.

- Update our tables of merged module counts to include crabs D and F for cbai_v2.0 and hemat_v1.6

- Make our PCA plots interactive, using plotly (Laura sent over some code)

- Set up a meeting with Grace/Steven/Claudia to discuss qPCR data and the reliability 

- Look specifically at immune-related genes. [This discussion post Steven made](https://github.com/RobertsLab/resources/discussions/1250) earlier would be helpful here.

- On a very different note, follow up with ADFG about obtaining survey data!

Alright, that's a whole lot of work, so I'd better get started on it! See you next time!




