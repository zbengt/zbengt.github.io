---
layout: post
title: Adaptive Immune Gene Dive
date: '2021-09-09'
categories: hematodinium
tags: immune genes
---

Previously, I identified [several genes with GO terms associated with an adaptive immune response](https://afcoyle.github.io/2021-08-18-immune_genes/). A writeup to each of these is described in the link, but here's the suspects:

- Q92956 (TNR14_HUMAN): Tumor necrosis factor receptor, interacts with CD160
- O95971 (BY55_HUMAN): CD160 antigen. Likely plays role in both anti-viral innate immune response along with adaptive immune response
- P60033 (CD81_HUMAN): CD81 antigen. Whole bunch of roles, inc. structural component, enables receptor complex assembly upon encounter with pathogens, and more.
    - NOTE: P35762 (CD81_MOUSE) was also present in the crab transcriptome.

Alright, so we've got genes linked to adaptive immunity! This raises two key questions. 
1. Are these genes actually related to adaptive immunity? 
2. Are these genes correctly identified?
3. Are these genes at a significant expression level?

1 and 2 are the big ones - the presence alone of adaptive immune-related genes would be a neat little finding. But 3 is also important to this system specifically. Let's go through each.

## Are the genes related to adaptive immunity?

TNR14_HUMAN: Did some searching, looks like there are a few papers describing it as potentially having a role in innate immune response (I think, wow it's a tough read), such as [this](https://www.jimmunol.org/content/191/2/828). Plus it also has the GO term associated with innate immune response. 
- Conclusion: NOT definitively adaptive-response

BY55_HUMAN: That same paper also examined CD160, and found that it is also linked to innate immune response. Plus a few others - it's pretty well-established. Plus again, also has the GO term associated with innate immune response.
- Conclusion: NOT definitively adaptive-response

CD81_HUMAN: No clear demonstrations I could find of its involvement in innate immune response. However, it also has [about a billion other  roles]https://www.uniprot.org/uniprot/P60033), including protein stability regulation, receptor internalization, protein localization, and a bunch more. 
- Conclusion: NOT definitively adaptive-response

Alright, well that's a bummer. For the sake of completeness and because it's quick, let's look at whether the genes were correctly IDed. Here's the e-values of each transcript matching a gene described above (plus CD81_MOUSE).

## Are the genes correctly identified?

| Transcript_ID             | Gene_ID     | e-value  |
|---------------------------|-------------|----------|
| TRINITY_DN277121_c0_g1_i2 | TNR14_HUMAN | 1.2e-171 |
| TRINITY_DN277121_c0_g1_i1 | TNR14_HUMAN | 2.3e-103 |
| TRINITY_DN472921_c0_g1_i1 | BY55_HUMAN  | 1.5e-18  |
| TRINITY_DN633565_c0_g1_i1 | CD81_HUMAN  | 1.5e-65  |
| TRINITY_DN39038_c3_g1_i2  | CD81_MOUSE  | 2.6e-131 |
| TRINITY_DN39038_c3_g1_i1  | CD81_MOUSE  | 1.7e-130 |

The BY55_HUMAN e-value looks a liiitle bit questionable, but all others are pretty dang solid. I'd feel confident in saying that genes were, in fact, IDed correctly

Alright, we just won't worry about Question 3 at this point - in light of Question 2, it's probably irrelevant. Ah well, onto the next thing!