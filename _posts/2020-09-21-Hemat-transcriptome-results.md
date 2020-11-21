---
layout: post
title: He-Mat and the Masters of the Database
date: '2020-09-21'
categories: onboarding Ubuntu hematodinium
tags: hematodinium Ubuntu
---

**Success!**

For some context, this post is a continuation of my entries from the past few days - I've been working to take hematodinium transcriptome data and locally blast against the Swiss-prot protein database

Well, I figured out a way to bring my files from WSL to Windows and link em in my notebook! My initial results [can be seen here](https://github.com/afcoyle/afcoyle-misc/blob/master/hematodinium/hemat-uniprot_blastx.tab)

**Important caveat**: The results above were produced following the same outline as the tutorial. Matches were made using max_target_seqs (set to 1). However, while investigating how to use blastx, I found info indicating that max_target_seqs just grabs the first sequence that meets the evalue standards (set to 1E-20 here). To me, that seemed sub-optimal, and so I reran blastx with a line added setting max_hsps to 1. [Those results are linked here](https://github.com/afcoyle/afcoyle-misc/blob/master/hematodinium/hemat_uniprot_maxhsps_blastx.tab)


After linking my results, I wrote my biography for the lab website. I then started on my next significant task - getting Jupyter Notebooks set up and linked to GitHub and re-doing the hematodinium transcriptome blast there! Made some progress, but am having some problems getting Jupyter Notebooks set up. Will update more tomorrow!

Since I was hitting a bit of a wall, I decided to put down that particular task for the night, and worked on a few hours of UW (non-SAFS) onboarding material that I'd just received. Their lessons weren't exactly the most inspired (one memorable slogan: "Theft - just don't do it"), but it's now complete! 



