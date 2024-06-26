---
layout: post
title: E5 3 Species 04 - Pocillopora Mapping Comparison
date: '2023-07-13'
categories: bioinformatics 3species
tags: E5 3species
---

Ran HISAT2 alignments for Pocillopora E5 data using P. effusa and P. verrucosa references. P. effusa had the worst alignment. P. verrucosa and P. meandrina are very similar in alignment rates, so we will go with P. meandrina since it has a better annotation.


# Alignment Results


## _P. effusa_

| Sample     | Total Reads | Paired Reads | Concordantly Aligned 0 Times | Concordantly Aligned 1 Time | Concordantly Aligned >1 Times | Discordantly Aligned 1 Time | Overall Alignment Rate |
|------------|-------------|--------------|-----------------------------|----------------------------|-------------------------------|-----------------------------|------------------------|
| Sample 1   | 54,219,233  | 54,219,233   | 25,580,382 (47.18%)        | 19,792,662 (36.50%)       | 8,846,189 (16.32%)           | 1,398,649 (5.47%)          | 62.33%                 |
| Sample 2   | 51,579,450  | 51,579,450   | 25,482,782 (49.40%)        | 17,573,943 (34.07%)       | 8,522,725 (16.52%)           | 1,095,747 (4.30%)          | 59.19%                 |
| Sample 3   | 55,302,459  | 55,302,459   | 27,995,152 (50.62%)        | 22,903,282 (41.41%)       | 4,404,025 (7.96%)            | 1,378,919 (4.93%)          | 59.72%                 |
| Sample 4   | 53,353,608  | 53,353,608   | 27,900,563 (52.29%)        | 22,025,395 (41.28%)       | 3,427,650 (6.42%)            | 1,504,147 (5.39%)          | 58.50%                 |
| Sample 5   | 42,513,488  | 42,513,488   | 26,242,576 (61.73%)        | 14,491,467 (34.09%)       | 1,779,445 (4.19%)            | 1,242,202 (4.73%)          | 48.40%                 |


## _P. verrucosa_

| Sample     | Total Reads | Paired Reads | Concordantly Aligned 0 Times | Concordantly Aligned 1 Time | Concordantly Aligned >1 Times | Discordantly Aligned 1 Time | Overall Alignment Rate |
|------------|-------------|--------------|-----------------------------|----------------------------|-------------------------------|-----------------------------|------------------------|
| Sample 1   | 54,219,233  | 54,219,233   | 26,646,141 (49.15%)        | 26,364,668 (48.63%)       | 1,208,424 (2.23%)           | 2,398,499 (9.00%)          | 63.73%                 |
| Sample 2   | 51,579,450  | 51,579,450   | 26,909,421 (52.17%)        | 23,565,675 (45.69%)       | 1,104,354 (2.14%)           | 1,931,846 (7.18%)          | 59.94%                 |
| Sample 3   | 55,302,459  | 55,302,459   | 26,695,037 (48.27%)        | 27,523,792 (49.77%)       | 1,083,630 (1.96%)           | 2,019,789 (7.57%)          | 62.96%                 |
| Sample 4   | 53,353,608  | 53,353,608   | 26,520,728 (49.71%)        | 25,805,042 (48.37%)       | 1,027,838 (1.93%)           | 2,110,685 (7.96%)          | 61.79%                 |
| Sample 5   | 42,513,488  | 42,513,488   | 24,714,378 (58.13%)        | 17,202,227 (40.46%)       | 596,883 (1.40%)             | 1,649,983 (6.68%)          | 52.22%                 |


## _P. meandrina_

| Sample     | Total Reads | Paired Reads | Concordantly Aligned 0 Times | Concordantly Aligned 1 Time | Concordantly Aligned >1 Times | Discordantly Aligned 1 Time | Overall Alignment Rate |
|------------|-------------|--------------|-----------------------------|----------------------------|-------------------------------|-----------------------------|------------------------|
| Sample 1   | 54,219,233  | 54,219,233   | 24,674,401 (45.51%)        | 24,699,780 (45.56%)       | 4,845,052 (8.94%)            | 1,657,119 (6.72%)          | 64.15%                 |
| Sample 2   | 51,579,450  | 51,579,450   | 24,586,231 (47.67%)        | 22,163,871 (42.97%)       | 4,829,348 (9.36%)            | 1,280,768 (5.21%)          | 60.90%                 |
| Sample 3   | 55,302,459  | 55,302,459   | 26,420,392 (47.77%)        | 26,136,411 (47.26%)       | 2,745,656 (4.96%)            | 1,564,617 (5.92%)          | 62.04%                 |
| Sample 4   | 53,353,608  | 53,353,608   | 26,586,350 (49.83%)        | 24,622,354 (46.15%)       | 2,144,904 (4.02%)            | 1,693,268 (6.37%)          | 60.58%                 |
| Sample 5   | 42,513,488  | 42,513,488   | 25,177,347 (59.22%)        | 16,194,243 (38.09%)       | 1,141,898 (2.69%)            | 1,400,533 (5.56%)          | 50.59%                 |








