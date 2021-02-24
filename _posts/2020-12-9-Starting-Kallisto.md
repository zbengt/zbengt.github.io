---
layout: post
title:  Starting Kallisto
date: '2020-12-9'
categories: hematodinium
tags: hematodinium, Trinity, Kallisto
---

**TRINITY Pipeline**

Well, there's good news and bad news. The good news is that I don't actually need to run the TRINITY pipeline to create a reference transcriptome - we already have some reference transcriptomes assembled! The bad news is...dang, that's a lot of work down the drain. As a result, I've shifted towards the downstream analyses found in the sidebar of [the Trinity manual](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification). 

General overview of next steps:
- Select which transcriptomes to examine - not choosing crab- or hemat-only (DONE)
    - cbai_hemat_transcriptome_v2.0: contains all pooled libraries & individual crab samples.
        - Encountering problems with running it through Kallisto (elaborated upon below)
    - cbai_hemat_transcriptome_v3.0: contains all pooled libraries
- Quantify transcript abundance with kallisto
    - Build index
    - Quantify abundance in each library
- Join with BLAST annotation tables in R or Bash

**Kallisto**
Well, Kallisto must've run into Hera, because it's definitely having some problems. 

I successfully built an index for transcriptome 3.0 (although it took around 36 hours), but my computer refused to build an index for transcriptome 2.0 - not enough RAM. I'm looking into getting Mox access (which I set aside after figuring out I wouldn't need it to create a reference transcriptome), but for the time being, I'm planning to move forward with my analysis using transcriptome 3.0, and come back to 2.0 when I can.

Again, I'm performing 4 total comparisons between different groups of infected crab
- Library 2 vs 4 (Day 2 low-temp vs. Day 2 high-temp)
- Library 2 vs 6 (Day 2 low-temp vs. Day 0 ambient)
- Library 4 vs 6 (Day 2 high-temp vs Day 0 ambient)
- Library 8 vs 10 (Day 2 all temps vs. Day 17 all temps)



