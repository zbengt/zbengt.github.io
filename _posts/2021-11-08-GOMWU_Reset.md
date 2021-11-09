---
layout: post
title:  Fixing GO-MWU
date: '2021-11-08'
categories: hematodinium
tags: GO-MWU
---

## GO-MWU Repair

Alright, so I'll keep this relatively quick, since it's quite late and I just spent a bunch of time fixing old mistakes.

So back in February when I was a bright-eyed and bushy-tailed young graduate student, I ran GO-MWU on a long series of pairwise comparisons for crab libraries aligned to each of the three transcriptomes. I've recently been going through my old work to put it together into paper format. In doing so, I ended up poking through the [GO-MWU README](https://github.com/z0on/GO_MWU/blob/master/README.md), and realized that GO-MWU made a really, **really** important update to their instructions.

In June, they added a short section to their FAQ saying the following: 

 `NOTE: In read-based gene expression analysis (RNA-seq, TagSeq) p-values may be biased towards highly abundant genes, especially when the read depth is low. This may result in the corresponding GO bias. Use log2-fold changes to avoid this.` 

 I had previously used unadjusted p-value, not log2-fold change. And this describes my data quite well. As a result, I spent most of tonight switching all my GO-MWU analyses over to log2-fold change. This resulted in the following:

 - Some substantial changes in results (following notation is treatment-day, so amb0 is ambient-group, day 0)
    - cbai_v2.0: More significant modules for all comparisons, including our amb0/amb2 and elev0/elev2
    - cbai_v4.0: Less significant modules for elev0/elev2, more for amb2/elev2 and amb0/amb2
    - hemat_v1.6: Lost all significant modules for all comparisons

- Adjusted some of the graphing parameters, as they're supposed to change when you switch over to using log2-fold change (described in each graphing chunk in both my files and the original repo)

- Revised and renamed the function that gets the p-values for GO-MWU. Since, after all, it's now getting log2 fold change

- While I was doing that, I also added some lines of code to the bottom of each chunk to auto-move files into the proper folders. Speeds up the process and minimizes error by just a bit

Alright, that's all I'm reporting. But really exceptionally glad that I went back and read the README again!