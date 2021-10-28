---
layout: post
title:  Background Info on Grace's Research
date: '2020-10-11'
categories: hematodinium background
tags: interview
---

**First, Some Project Background**

In order to get some background and understanding of the hematodinium work done in the Roberts lab, I chatted with Grace Crandall via Zoom for an hour or so. 

As background, Grace is a member of the Roberts lab who just graduated. She's been investigating how infected and uninfected Tanner crab respond to various temperature regimes by collecting and analyzing transcriptomes. Here's a (very very) rough outline of the experimental protocol, which was conducted in partnership with Pam Jensen at the Alaska Department of Fish and Game (ADF&G):

- 400 male Tanner crab were collected from Stephen's Passage, near Juneau AK, and brought back to the ADF&G lab in Juneau

- All were tested for hematodinium infection via cPCR, 180 were selected for temperature treatments (90 infected, 90 uninfected)

- Crabs were divided into 9 tanks, with 10 infected and 10 uninfected in each

- Tanks were at 3 different temperatures - 4°C (cool), 7.5°C (ambient), and 10°C (elevated)

- Hemolymph samples were taken 3 times - initial sample, Day 2, and Day 17, with Day 0 being the start of the differential temp treatments

- RNA was extracted from hemolymph samples

- Gene expression differences over time and between treatments were analyzed

**Notes on Conversation with Grace**

The following notes are in no particular order


- Grant money for project was released in fall 2017

- Collections targeted immatures. Since crabs are more prone to infection while molting, and Tanner crabs have a terminal molt, immature crabs have a higher chance of having BCS

- There's a new version of the maturity ratio (chela height vs carapace width, used to determine maturity in Tanner crab) that's specific to SE Alaska

- Crabs were left in tanks for ~9 days to acclimate

- Crabs weren't fed for two reasons. First, it would've been an extra variable to manage. Second, because water quality was an issue

- Prevalence of hematodinium in Kodiak/Cook Inlet (interested me because it's between the high-infection areas of SE Alaska and the low-infection-rate Bering Sea): Grace isn't sure, will send over papers

- All tanks had header tanks that were fed into main tanks for temperature treatment. No water mixing between tanks

- Some initially-uninfected crabs had extremely low infection levels at end of experiment. Likely because cPCR was used initially, qPCR was used at the end

- There's a ton of good qPCR data that's largely unused - that'll be great for looking at the hematodinium side of things. How much hemat is present in which crabs? What temp treatment were they?

- No hemolymph samples were taken from dead crab

- Heater likely broke during 9-day acclimation period, but doesn't remember exactly

- Several samples were taken from each crab at Day 17 as a safety net. 3 samples each from cold and ambient crab, 6 each from warm crab (since 95% of warm-treatment crab died by Day 10)

- There are still a ton of hemolymph samples in the -80 freezer

- Hematodinium/crab comparisons done with MEGAN (which assigns phylogenetic stuff to sequence data) and Busco. Steps on github wiki

- Sequenced the crabs that fit with libraries we wanted. As am example: wanted pool of samples with Day 2 uninfected crabs. Picked ones that matched it with good RNA yields

- Morado and Jeff Fields have a lot of good hematodinium research

- Manifestations of hematodinium infection differ some between hosts. BCD might be unique for its bitter taste

- Samples in the -80 freezer don't have RNALater anymore - they were pelleted and the supernatant was extracted

- In Grace's paper draft, Supplemental I is all samples, Supplemental II is libraries

- Later, figured how to do differential gene expression using package DESeq2

- Certain stages of hematodinium may be more infectious than others. Dinospore may be more infectious

- Libraries were limited by time, since so much time was spent figuring out extraction protocol


**Biggest Obstacles**

- The biggest problem was finding a good method of RNA extraction that wouldn't result in degraded or impure RNA. This problem took a year to solve! Eventually Pam Jensen found a method that worked. 

- GeneWiz was used for sequencing, since they could sequence smaller amounts of RNA and were much quicker

- The large die-off (95% of warm-treatment crabs dying by Day 10) limited what could be analyzed

- The NW Genomic Center (UW facility in Figgy Hall) took over 6 months to sequence samples

- Since crab hemolymph is clear for uninfected crabs, it was sometimes difficult to determine whether you were getting seawater or hemolymph

- When trying to pellet samples, was getting more of a slush. This might have been why the RNA wasn't extracting well. After switching to a hemolymph extraction protocol, this was largely resolved


**Next Steps**

- Looking at hematodinium. Have only looked at host response to temperature. What does hematodinium do inside host? Is it impacted by temperature?

- Hemolymph was placed in RNALater, spun down, supernatant was pulled off and kept. Might be something to be extracted in there

- Getting more of a time series. Won't be super robust for warm treatment. Already partially done, but a more statistic-focused approach (on either crab or hematodinium) would be great