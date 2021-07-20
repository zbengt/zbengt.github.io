---
layout: post
title:  Back to Square One?
date: '2021-02-09'
categories: hematodinium
tags: kallisto, DESeq2, GO-MWU
---

**My analysis outline? Incorrect? It's more likely than you think**

Last time we checked in, I had finally ran GO-MWU on my four analyses, and came up with more or less nothing. However, after running GO-MWU, I realized I had made an error when determining what analyses to perform. My analyses largely involved comparing different temperature treatments (elevated vs. ambient vs. low) by pooling together samples from Day 0 and Day 2. For instance, to compare elevated and low, I pooled all elevated crab from Day 0+2, and compared gene expression against all lowered-temp crab from Day 0+2. 

However, as it turns out, Day 0 samples were taken prior to any exposure to different temperatures! 

As a result, I decided to redo my whole analysis from scratch. This wasn't as big of a downside as you may think - I started this project back when I was a tiny little grad student who knew nothing (as opposed to current me, a grad student with 1.5 whole quarters of knowledge). Therefore, I made some fairly simple mistakes that would've caused me to eventually need to redo my analysis. For instance, I didn't use checksums to verify my files downloaded properly - at the time, I didn't even know what checksums were. 

So I ditched my 4 examinations of differential expression and replaced them with the following:

1: Ambient Day 0+2 vs Elevated Day 2 (individual libraries only)

| Crab ID | Library ID | Day | Temperature at time |
|---------|------------|-----|---------------------|
| G       | 272        | 2   | Elevated            |
| H       | 294        | 2   | Elevated            |
| I       | 280        | 2   | Elevated            |
| A       | 178        | 0   | Ambient             |
| A       | 359        | 2   | Ambient             |
| B       | 118        | 0   | Ambient             |
| B       | 349        | 2   | Ambient             |
| C       | 132        | 0   | Ambient             |
| C       | 334        | 2   | Ambient             |

2: All samples taken from ambient-temperature crabs vs. all samples taken from elevated-temperature crabs. 
**IMPORTANT: Day 0 samples from crabs in the elevated/lowered-treatment group are counted as ambient-temperature, since they had not yet been exposed to different water temperatures. Those samples are shown in _italics_ **

| Crab ID | Library ID | Day | Temperature at time |
|---------|------------|-----|---------------------|
| G       | 272        | 2   | Elevated            |
| H       | 294        | 2   | Elevated            |
| I       | 280        | 2   | Elevated            |
| pooled  | 380825     | 2   | Elevated            |
| *G*       | *173*        | *0*   |*Ambient*             |
| *H*       | *72*         | *0*   |*Ambient*            |
| *I*       | *127*        | *0*   | *Ambient*             |
| A       | 178        | 0   | Ambient             |
| A       | 359        | 2   | Ambient             |
| A       | 463        | 17  | Ambient             |
| B       | 118        | 0   | Ambient             |
| B       | 349        | 2   | Ambient             |
| B       | 481        | 17  | Ambient             |
| C       | 132        | 0   | Ambient             |
| C       | 334        | 2   | Ambient             |
| C       | 485        | 17  | Ambient             |
| *E*       | *151*        | *0*   | *Ambient*             |
| pooled  | 380820     | 0   | Ambient             |


**Reference: Here is a table of all libraries of infected crabs that have been sequenced**

| Crab ID | Library ID | Day | Temperature at time |
|---------|------------|-----|---------------------|
| G       | 272        | 2   | Elevated            |
| H       | 294        | 2   | Elevated            |
| I       | 280        | 2   | Elevated            |
| pooled  | 380825     | 2   | Elevated            |
| *G*       | *173*        | *0*   |*Ambient*             |
| *H*       | *72*         | *0*   |*Ambient*            |
| *I*       | *127*        | *0*   | *Ambient*             |
| A       | 178        | 0   | Ambient             |
| A       | 359        | 2   | Ambient             |
| A       | 463        | 17  | Ambient             |
| B       | 118        | 0   | Ambient             |
| B       | 349        | 2   | Ambient             |
| B       | 481        | 17  | Ambient             |
| C       | 132        | 0   | Ambient             |
| C       | 334        | 2   | Ambient             |
| C       | 485        | 17  | Ambient             |
| *E*       | *151*        | *0*   | *Ambient*             |
| pooled  | 380820     | 0   | Ambient             |
| E       | 254        | 2   | Lowered             |
| E       | 445        | 17  | Lowered             |
| pooled  | 380823     | 2   | Lowered             |



I then started running these 2 analyses through the pipeline. I got onto Roadrunner, downloaded all files, checked they downloaded correctly, and then used kallisto to get the counts of each library. Interestingly, while Transcriptome 3.0 (which was used previously) seemed to produce more matches and longer fragment lengths, I chose to use Transcriptome 2.0, as it was built from all libraries, while 3.0 was built from a relatively small subset of libraries). We previously used 3.0, as it proved impossible to build a Kallisto index for 2.0 on a local machine. But now that I'm on Roadrunner, it's not an issue!

I then used the Trinity script to create a matrix of counts, and then ran that matrix through DESeq2. [Here's a link to my DESeq2 output](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/DESeq2_output/cbai_transcriptomev2.0). I then ran my R script to get a newline-separated file of accession IDs.

The next step is to get all my GO IDs. I might try to find a way past that - potentially using my recently-created tables of accession IDs and GO IDs. Failing that, I may see if I can use Roadrunner to run Sam's shell script instead of my local machine. By my estimation, it'll take about 4 days to run each analysis on my local machine, so I'd like to avoid that if possible!
