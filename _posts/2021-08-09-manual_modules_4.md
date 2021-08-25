---
layout: post
title: Manual Modules 4
date: '2021-08-09'
categories: hematodinium
tags: manual clustering
---

Back in mid-July, I manually created modules to cluster genes by expression patterns. Again, "manual" oversells it a bit - I really just provided a cut height, looked at the resulting clusters, and then named them by overall expression patterns. A recap of that naming system can be found in my lab notebook entry on 2021-07-12 (named Manual Modules 3)

Anyway, clusters were created for each individual crab. However, after some discussion, we concluded that the previous bar for genes was a bit low. Previously, our bar was taking the libraries for the individual crab (say, the 3 libraries for Crab A), and only including genes with 5+ counts among all libraries. So for instance, a gene with 4 counts on Day 0, 0 on Day 2, and 1 on Day 17 would be included. Of course, the downside of this is that a) that's a really low bar, so we were looking at expression that probably wasn't meaningful, and b) we were ignoring genes that might be highly-expressed in other crab

As a result, we changed the bar - the new one is only including genes with 30+ counts in at least one library for ANY crab, and at least 1 count in that particular crab. So a Crab A gene with 1 count on Day 0 and 0 on Day 2 and 17 (but 35 counts on Day 2 for Crab B) would be included for Crab A.

Here are our results! I'll post the tables for percentages, along with a link to the actual files for both percentages and raw counts.

#### cbai_transcriptomev2.0

[Percentages](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/cbai_transcriptomev2.0/all_genes/merged_modules_percentages_table.csv), [Raw Counts](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/cbai_transcriptomev2.0/all_genes/merged_modules_raw_counts.txt)

| Crab   | HLH   | HTL   | LHL   | LTH   | HL    | LH    | MIX   |
|--------|-------|-------|-------|-------|-------|-------|-------|
| Crab_A | 0.008 | 0.242 | 0.27  | 0.481 | NA    | NA    | NA    |
| Crab_B | 0.096 | 0.363 | 0.336 | 0.205 | NA    | NA    | NA    |
| Crab_C | 0.23  | 0.412 | 0.12  | 0.238 | NA    | NA    | NA    |
| Crab_D | 0.122 | 0.318 | 0.32  | 0.24  | NA    | NA    | NA    |
| Crab_E | 0.071 | 0.383 | 0.267 | 0.279 | NA    | NA    | NA    |
| Crab_F | 0.007 | 0.585 | 0.103 | 0.305 | NA    | NA    | NA    |
| Crab_G | NA    | NA    | NA    | NA    | 0.043 | 0.233 | 0.724 |
| Crab_H | NA    | NA    | NA    | NA    | 0.086 | 0.138 | 0.776 |
| Crab_I | NA    | NA    | NA    | NA    | 0.31  | 0.187 | 0.503 |

#### cbai_transcriptomev4.0

[Percentages](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/cbai_transcriptomev4.0/all_genes/merged_modules_percentages_table.csv), [Raw Counts](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/cbai_transcriptomev4.0/all_genes/merged_modules_counts_table.csv)

| Crab   | HLH   | HTL   | LHL   | LTH   | MIX   | HL    | LH    |
|--------|-------|-------|-------|-------|-------|-------|-------|
| Crab_A | 0.006 | 0.329 | 0.355 | 0.309 | NA    | NA    | NA    |
| Crab_B | 0.219 | 0.24  | 0.342 | 0.199 | NA    | NA    | NA    |
| Crab_C | 0.125 | 0.384 | 0.133 | 0.358 | NA    | NA    | NA    |
| Crab_D | 0.057 | 0.051 | 0.388 | 0.505 | NA    | NA    | NA    |
| Crab_E | 0.086 | 0.331 | 0.39  | 0.192 | NA    | NA    | NA    |
| Crab_F | 0     | 0.592 | 0.246 | 0.162 | 0     | NA    | NA    |
| Crab_G | NA    | NA    | NA    | NA    | 0.519 | 0.039 | 0.442 |
| Crab_H | NA    | NA    | NA    | NA    | 0.978 | NA    | 0.022 |
| Crab_I | NA    | NA    | NA    | NA    | 0.959 | 0.013 | 0.028 |

#### hemat_transcriptomev1.6

[Percentages](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/hemat_transcriptomev1.6/all_genes/merged_modules_percentages_table.csv), [Raw Counts](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/manual_clustering/hemat_transcriptomev1.6/all_genes/merged_modules_counts_table.csv)

| Crab   | HTL   | LHL   | LTH   | HLH   | HL    | MIX   | LH    |
|--------|-------|-------|-------|-------|-------|-------|-------|
| Crab_A | 0.304 | 0.313 | 0.383 | NA    | NA    | NA    | NA    |
| Crab_B | 0.453 | 0.383 | 0.125 | 0.039 | NA    | NA    | NA    |
| Crab_C | 0.128 | 0.208 | 0.46  | 0.204 | NA    | NA    | NA    |
| Crab_D | 0.164 | 0.288 | 0.377 | 0.171 | NA    | NA    | NA    |
| Crab_E | 0.312 | 0.249 | 0.439 | 0     | NA    | NA    | NA    |
| Crab_F | 0.397 | 0.081 | 0.485 | 0.037 | NA    | NA    | NA    |
| Crab_G | NA    | NA    | NA    | NA    | 0.009 | 0.991 | NA    |
| Crab_H | NA    | NA    | NA    | NA    | NA    | 1     | NA    |
| Crab_I | NA    | NA    | NA    | NA    | NA    | 0.869 | 0.131 |

### Observations

- For cbai_v2.0, Crabs D-F (lowered-temp treatment group) seem to have more HTL and LHL genes than Crabs A-C (ambient-temperature treatment group). This also seems to be true for cbai_v4.0, but does NOT seem to be true for hemat_v1.6

- Overall, it doesn't look like the variation is quite as extreme as it was with our previous method. 

- Crabs G-I really don't have much useful data - too many MIX to indicate much of anything

- Generally, not really seeing too much of an established pattern here