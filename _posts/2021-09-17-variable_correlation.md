---
layout: post
title: Correlating Crab Variables
date: '2021-09-17'
categories: hematodinium
tags: Variable Correlation
---

NOTE: As of 9/27, this post is outdated. Check the next one!

A pretty basic question that I haven't really answered yet is how different variables are actually correlated with each other! So to figure this out, I decided to make a correlation matrix. Each variable applies to one individual crab. So we're leaving the time dimension out of it, and are just examining covariance between crabs. 

I made [this table](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/data/indiv_crab_summary.csv) a while back. Each row is a value for each individual crab. 

As a reminder, the three transcriptomes used for alignment in this experiment are:
- cbai_transcriptomev2.0: unfiltered by taxa
- cbai_transcriptomev4.0: filtered, presumably only _C. bairdi_ seqs
- hemat_transcriptomev1.6: filtered, presumably only _Hematodinium_ seqs

The following columns are present:

- Crab ID (A-I)

- Uniq_ID: This links the notation used in Grace's analysis (and mine) with some of the earlier Jensen tables

- Treatment: The temperature treatment the crab was exposed to

- Timepoints:The number of timepoints (3 for ambient- and lowered-temperature treatment groups, 2 for elevated)

- cPCR_Positive: Did the crab test positive for _Hematodinium_ with conventional PCR?

- qPCR_SQ_mean (one column for each of 2 timepoints). The qPCR results when tested for _Hematodinium_. Higher = more _Hematodinium_ present

- Plasticity (three columns, one for each transcriptome): I created PCAs with each aligned library as a single plot. This value is the mean distance between each individual library for that crab and the overall centroid for that crab.

- Pct Aligned (three columns, one for each transcriptome): Using kallisto, I aligned libraries to the transcriptome. This is the mean percent aligned for libraries belonging to that crab.

- Transcripts (three columns, one for each transcriptome): The total number of expressed transcripts in all samples for that crab. (e.g. if Transcript 1 is expressed in all samples for Crab A, it's still only counted as one transcript).

- Immune Transcripts (three columns, one for each transcriptome): Same as above, but ONLY counting transcripts associated with the GO term for immune response

- Carapace Width: The width of the crab carapace (not inclusive of spines) 

- Shell Condition: The level of wear on the shell. Roughly corresponds to time since molt

- Chela Height: The height of the right claw (chela). The ratio of carapace width to chela height indicates whether the crab is mature

- Imm_Mature: Whether the crab is immature or mature

- Death: The date of death, if applicable.

I took all these columns and created a correlation matrix, removing the ones that didn't make sense in that context (such as Uniq_ID). The script I used can be found [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/10_2_summary_corr_matrix.Rmd). I created two separate visualizations - a dot plot, available [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/correlation/all_variables_corr_dot_plot.png), and a chart, available [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/correlation/all_variables_corr_chart.png). The image is pretty dang small, and in the chart not all variables are easily readable, but they're in the same order as in the dot plot (which has its variables in nearly the same order as the original table at the top of this post).

I'll do a quick dive into the results here!

### Correlations Found:

I'll do these in order from highest to lowest by absolute value (i.e. negative correlations will also be noted). Ones in **bold** are the interesting ones!

- hemat_v1.6 transcripts and hemat_v1.6 immune transcripts (0.98). This makes a good deal of sense - the more transcripts a crab has, the higher chance of immune transcripts. But not very meaningful. Furthermore, probably worth disregarding regardless - the max for an individual crab's hemat. immune transcripts was 5

- cbai_v2.0 transcripts and cbai_v4.0 transcripts (0.88). The cbai_v4.0 transcriptome is derived from the cbai_v2.0 transcriptome, and cbai_v2.0 is predominantly _C. bairdi_ genes. Not interesting.

- **cbai_v2.0 plasticity** and **hemat_v1.6 pct aligned** (-0.85). Well well well, we have an interesting one! The higher the percentage of the libraries that was aligned to the _Hematodinium_ transcriptome, the lower the overall plasticity was. However, note that this does **NOT** apply to cbai_v4.0 plasticity and hemat. percent aligned (which was -0.50). Ultimately interesting. But unless we over-filtered to create our _C. bairdi_ transcriptome, no clearly biologically meaningful link I see.

- **hemat_v1.6 plasticity** and **cbai_v2.0 transcripts** (0.83). A higher number of unfiltered transcripts corresponds to a higher _Hematodinium_ plasticity. This is pretty dang neat, and could be something biologically real - higher gene expression means that _Hematodinium_ might be varying things up more. Of course, cbai_v2.0 is partially made up of _Hematodinium_ genes (although only a very small percentage), so a higher plasticity does directly contribute to a higher number of transcripts. However, it's definitely a cool find!

- cbai_v2.0 plasticity and cbai_v4.0 plasticity (0.82). For the same reason as discussed with cbai_v2.0 and cbai_v4.0 transcripts (cbai_v4.0 is derived from v2.0), not really interesting.

- cbai_v2.0 immune transcripts and cbai_v4.0 immune transcripts (0.81). Again, boring - we know cbai_v2.0 and v4.0 track pretty closely

- Chela height and maturity (0.75). We only have 1 immature crab, so really we can discard everything related to maturity

- Chela height and carapace width (0.74). Front page of tomorrow's newspaper: "BIGGER CRABS HAVE BIGGER CLAWS"

- **hemat_v1.6 transcripts** and **cbai_v4.0 immune transcripts** (0.73). Now we're talkin'! Not anything particularly novel - the more _Hematodinium_ expression, the more immune activity from _C. bairdi_. But neat to see nonetheless!

- hemat_v1.6 plasticity and hemat_v1.6 immune transcripts (0.72): Makes sense - a higher plasticity means more variance in infection, which means more overall immune-related transcripts expressed.

- **cbai_v4.0 immune transcripts** and **hemat_v1.6 immune transcripts** (0.71): Makes sense - more activity on the part of the parasite (or host) means more activity on the part of the other! However, take this with a grain of salt - the number of hemat_v1.6 immune transcripts maxes out at 5.

- hemat_v1.6 transcripts and cbai_v2.0 immune transcripts (0.70). This isn't necessarily _boring_, but is already covered by hemat_v1.6 transcripts vs. cbai_v4.0 immune transcripts (3 points above this entry).

- cbai_v2.0 pct aligned and hemat_v1.6 pct aligned (0.68): Since hemat_v1.6 is derived from cbai_v2.0 (although only constitutes a very small part of it), this isn't really that interesting.

- Chela height and treatment group (-0.68): Hmm, sounds like the distribution of sampled crab size wasn't uniform among teratments! Was lower for carapace width (-0.53) so probably nothing to worry about.

- cbai_v4.0 pct aligned and cbai_v2.0 transcripts (-0.67): Again, not really interested in v2.0 vs v4.0 measurements

Alright, we'll stop there - that's everything with a correlation greater than 0.66. A bit arbitrary, but I checked everything down to 0.6, and there's nothing below that bar that's interesting.

### Correlations NOT found

These are our problem children - where we SHOULD be finding correlations but are not.

- **cPCR_positive** and **qPCR_SQ_mean**. This is the big one. According to cPCR, we have two uninfected crabs. According to qPCR, we have zero uninfected crabs - just 4 with low-level infections and 5 with heavy infections on Day 0. However, some of those crabs with heavy infections suddenly have low-level infections on Day 17. And to make matters far worse, one of the crab deemed uninfected by cPCR has the heaviest infection of all on Day 0, followed by the lightest on Day 17. Our correlation here is a measly 0.33

- **cPCR_positive** and **ANYTHING RELATED TO GENE EXPRESSION**. Ideally, we'd see some differences between uninfected crab and infected crab in one of the many ways we've examined gene expression. But...nope. Nothing.

- **qPCR_SQ_mean** and **ANYTHING RELATED TO GENE EXPRESSION**. Well, qPCR isn't any better, sadly. Zero correlation to any gene expression.

- **Treatment** and **ANYTHING RELATED TO GENE EXPRESSION**. Again, a big swing and a miss here. No correlations between temperature treatment and anything else here.

Gotta be honest, these aren't exactly the most encouraging results in the world. But we do have some neat-ish stuff in here! 