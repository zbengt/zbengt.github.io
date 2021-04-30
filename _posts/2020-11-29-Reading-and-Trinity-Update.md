---
layout: post
title:  TRINITY Pipeline Update
date: '2020-11-29'
categories: hematodinium, background
tags: hematodinium, Trinity, notes
---

**Reading Update**

Continued to dive into the reading, and have some additional potential research ideas! For easier navigation, those have been added to my previous update, which also focuses on research ideas.

**TRINITY Update**

For the past few weeks, I've been trying to examine differential gene expression in infected samples taken during Grace's experiment, with the goal of finding all differentially-expressed genes for both hematodinium and the crab. This will be accomplished by putting the crabs through the TRINITY pipeline - for a thorough grounding in what the pipeline does, check out [Haas et al. 2013](https://www.researchgate.net/publication/279835240_De_novo_transcript_sequence_reconstruction_from_RNA-Seq_reference_generation_and_analysis_with_Trinity)

So far my main difficulties have been in actually getting the pipeline to run - working through Ubuntu is tricky! But I got it running eventually by following the [official instructions](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity), and when that failed or was unclear, using Bernadette Johnson's walkthrough instructions, available [here](https://bernadettebiology.weebly.com/protocols--tutorials.html). 

I'm performing 4 total comparisons between different groups of infected crab
- Library 2 vs 4 (Day 2 low-temp vs. Day 2 high-temp)
- Library 2 vs 6 (Day 2 low-temp vs. Day 0 ambient)
- Library 4 vs 6 (Day 2 high-temp vs Day 0 ambient)
- Library 8 vs 10 (Day 2 all temps vs. Day 17 all temps)

Ideally would have had a Day 2 ambient library, but no such library existed. And since the temperature wasn't actually changed for the ambient group, transcriptomes from Day 0 should be more or less the same as those from Day 2.

Library info found [here](https://github.com/RobertsLab/paper-tanner-crab/blob/master/supplementary-information/S1-RNAsequ-libraries.csv), sequences gathered from Grace's Google Sheets master doc (title is C_bairdi_RNAseq)

The Trinity pipeline takes an extremely, extremely long time to run (upwards of 80 hours per sample). I'm still waiting for my first run to finish, but subsequent comparisons will be made using the lab's computing resources. First I need to figure out how to remote into them - that's the next task ahead!

