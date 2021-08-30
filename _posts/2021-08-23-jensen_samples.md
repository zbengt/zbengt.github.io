---
layout: post
title: Jensen Samples
date: '2021-08-30'
categories: hematodinium
tags: archived samples, ADF&G
---

Alright, at this point, I've got a few different projects going. First, I'm examining gene expression within crabs infected with Hematodinium (using Grace Crandall's samples). Second, I'm working on modeling a few decades of survey data on Hematodinium infection rates from the Alaska Department of Fish and Game (I've been meaning to write up a lab notebook post on that, will do it soon). 

However, recently a third project got added to the pile! During Hack Week, we went into the lab to tidy things up. While there, I found a HUGE number of archived samples that were sent over by Pam Jensen after her retirement. They generally consist of hemolymph samples in deep 96-well plates, and appear to contain samples from both infected and uninfected crabs. Samples date back to the mid-2000s, and there's literally four large totes full of them.

This is an absolute goldmine of data, and I really really want to do something exciting with it. However, the first step of that is actually determining what the samples are. [Here's a list of all boxes](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/data/jensen_archived_samples/jensen_archived_samples.csv). However, it's missing a lot of info - which wells map to which crabs? When and where was each crab sampled? Hopefully, those questions can be answered by [the data that Pam sent over](https://gannet.fish.washington.edu/hematodinium/), particularly the Access databases. My plan is to take the following steps:

- Merge tables in Access databases to map each well to a specific crab with collection and sampling time info. This should give us data through 2013.

- Merge all Excel tables (from 2014-2019) to also map each well to a specific crab, as above.

- See if the two merged tables can also be merged without losing substantial information

- Look at the resulting table(s) to see what interesting questions can be answered with this! I'm particularly angling towards coming up with a GRFP topic - since I'm really the only one with access to these samples, it seems like a pretty appealing project.

