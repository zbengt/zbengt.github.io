---
layout: post
title:  Adding Individual Libraries
date: '2020-12-17'
categories: hematodinium
tags: hematodinium, DESeq2, Kallisto
---

**DESeq2**

First, I finished up my first quarter! Congratulations to me! 

Alright, I've made a lot of progress since last time. For the moment, building an index for transcriptome 3.0 is on the backburner - the top priority is just getting a list of differentially-expressed genes. I was able to run Kallisto and produce abundance files for all of the 5 libraries I'm examining. Then I used [this script from Trinity - under "Build Transcript and Gene Expression Matrices -"](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification) to turn those abundance files into a single matrix. [Link to my script](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/Scripts/01_cbai_hemat_pooled_libraries-Copy1.ipynb)

I then read my matrix into R, and analyzed it with DESeq2- [link to my script](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/scripts/01_cbai_hemat_pooled_libraries.ipynb). Now, here's where I ran into problems. As a reminder, here are the libraries I was examining: 
- Library 2 (Day 2 low-temp)
- Library 4 (Day 2 high-temp)
- Library 6 (Day 0 ambient)
- Library 8 (Day 2, all temps)
- Library 10 (Day 17, all temps)

You'll notice that there's practically no replication there. And that means that DESeq2 (or edgeR or really any method of differential analysis) can't look to see what genes are statistically significant in their difference, because they can't look at the variance. And I'd really like to look at which genes are statistically different between temperature/time treatments

Solution is to go forward with repeating the protocol. It should be much, much, much quicker this time, for several reasons. First, I actually know what I'm doing. Second, I've already got my transcriptome matrix built for Kallisto. It'll still be a day or two, but that's a far cry from...however long this has been. I'll be examining the following samples, which can be downloaded either from Grace's google sheet or [Gannet](https://gannet.fish.washington.edu/Atumefaciens/20200318_cbai_RNAseq_fastp_trimming).

Just as before, all samples came from infected crab. The only individual infected crab sample that will not be part of this analysis is 445 (temp = lowered, day = 17), as there are no replicates.

| ID  | Temp  | Day  |
|---|---|---|
| 127 | Elevated  | 0  |
| 173 | Elevated  | 0  |
| 72  | Elevated  | 0  |
| 272 | Elevated  | 2  |
| 280 | Elevated  | 2  |
| 294 | Elevated  | 2  |
| 118 | Ambient   | 0  |
| 132 | Ambient   | 0  |
| 178 | Ambient   | 0  |
| 334 | Ambient   | 2  |
| 349 | Ambient   | 2  |
| 359 | Ambient   | 2  |
| 463 | Ambient   | 17 |
| 481 | Ambient   | 17 |
| 485 | Ambient   | 17 |
| 151 | Lowered   | 2  |
| 254 | Lowered   | 2  |

4 total comparisons will be made 
- Day 0 vs Day 17 Ambient (118/132/178 vs 463/481/485)
- Day 0+2 Elevated vs Low (127/173/72/272/280/294 vs 151/254)
- Day 0+2 Elevated vs Ambient (127/173/72/272/280/294 vs 118/132/178/334/349/359)
- Day 0+2 Ambient vs Low (118/132/178/334/349/359 vs 151/254)

Well, there's good news and bad news. The good news is that I don't actually need to run the TRINITY pipeline to create a reference transcriptome - we already have some reference transcriptomes assembled! The bad news is...dang, that's a lot of work down the drain. As a result, I've shifted towards the downstream analyses found in the sidebar of [the Trinity manual](https://github.com/trinityrnaseq/trinityrnaseq/wiki). 

General overview of next steps:
- Quantify transcript abundance for Day 0/Day 17 Ambient with Kallisto
- Merge abundances with Trinity script linked above
- Analyze with DESeq2 as linked above
- Repeat with temperature comparisons (if local storage allows for all to be done simultaneously)





