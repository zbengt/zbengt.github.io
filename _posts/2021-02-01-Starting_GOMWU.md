---
layout: post
title:  Starting GO-MWU
date: '2021-02-01'
categories: hematodinium
tags: GO-MWU
---

**Getting GO-MWU Input**
Alright, I finally progressed to the next stage - running [GO-MWU](https://github.com/z0on/GO_MWU)!

In order to run GO-MWU, I needed two files. First, a two-column CSV table with gene IDs and a measure of significance (I used unadjusted p-value). Second, a two-column tab-separated file with gene IDs and GO terms.

Obtaining the first was fairly straightforward - my R script can be found [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/1_5_GO-MWU_prep.Rmd). In that script, I pulled all genes - including those with an unadjusted p-value of NA. That might be relevant later - may want to redo with only non-NA values.

Obtaining the second was also straightforward, but took a long, long time. I had previously created a newline-separated file of accession IDs for each of my 4 comparisons (reminder: elevated vs ambient, elev. vs low, amb. vs low, day 0 vs day 17 amb.). I plugged that into a [script from Sam](https://github.com/RobertsLab/code/blob/master/script-box/uniprot2go.sh), which got the GO terms for each accession ID. However, it took an extremely long time to run - 1-2 days per file. I tried running it on Roadrunner, but failed - believe the problem was that I didn't have the permissions to connect to a URL. Might've been fixable, but I only had one file left to run, so I just ran it on my local machine. In summation, this produced a two-column tab-separated file with gene IDs and GO terms - almost exactly what GO-MWU needs as input!

But not quite! GO-MWU requires just one line per gene, and my files had many genes on multiple lines. Luckily, there was an easy solution. I used the perl script nrify_GOtable.pl included with GO-MWU to reduce to one line per gene. This gave me all the inputs necessary for GO-MWU.

When I finally ran GO-MWU, the results were somewhat underwhelming. Of my 4 comparisons, 3 had no GO terms pass 10% FDR, and the fourth (Elevated vs. Ambient) had only 3. Of those 3 GO terms, no categories had a p-value below 0.05. Maybe the problem was my inputs - I might've screwed up by including all genes with an unadjusted p-value of NA. Worth checking later to fix. Also need to do a lot more analysis - [Yaamini's script](https://github.com/epigeneticstoocean/paper-gonad-meth/blob/master/code/14-Gene-Enrichment-with-GO-MWU.Rmd) has some helpful analyses of post-GO-MWU results. I think that'll be the next step, but stay tuned!
