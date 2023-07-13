---
layout: post
title: E5 3 Species 03 - Pocillopora Mapping
date: '2023-07-10'
categories: bioinformatics 3species
tags: E5 3species
---

I ran HISAT2 to get alignments for the E5 Pocillopora data using the _P. meandrina_ reference we had decided on [here](https://github.com/urol-e5/deep-dive/wiki/Species-Characteristics-and-Genomic-Resources) in the wiki. Overall alignment rates look kind of low. Are we sure about using this?

## Alignment Summary


| Reads        | Aligned concordantly 0 times | Aligned concordantly exactly 1 time | Aligned concordantly >1 times | Overall Alignment Rate |
|--------------|------------------------------|-------------------------------------|-------------------------------|------------------------|
| 54,219,233   | 24,674,401 (45.51%)         | 24,699,780 (45.56%)                 | 4,845,052 (8.94%)             | 64.15%                 |
| 51,579,450   | 24,586,231 (47.67%)         | 22,163,871 (42.97%)                 | 4,829,348 (9.36%)             | 60.90%                 |
| 55,302,459   | 26,420,392 (47.77%)         | 26,136,411 (47.26%)                 | 2,745,656 (4.96%)             | 62.04%                 |
| 53,353,608   | 26,586,350 (49.83%)         | 24,622,354 (46.15%)                 | 2,144,904 (4.02%)             | 60.58%                 |
| 42,513,488   | 25,177,347 (59.22%)         | 16,194,243 (38.09%)                 | 1,141,898 (2.69%)             | 50.59%                 |


