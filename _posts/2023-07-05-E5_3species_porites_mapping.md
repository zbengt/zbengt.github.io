---
layout: post
title: E5 3 Species 02 - Porites Mapping
date: '2023-01-05'
categories: bioinformatics 3species
tags: E5 3species
---

I ran HISAT2 to get alignments for the E5 Porites data using _P. lobata_ and _P. evermanni_ references from [here](https://www.genoscope.cns.fr/corals/genomes.html). Looks like the _P. evermanni_ reference wins out with higher overall alignments.

## _Porites lobata_

| Run | Paired Reads | Aligned Concordantly 0 Times | Aligned Concordantly Exactly 1 Time | Aligned Concordantly >1 Times | Aligned Discordantly 1 Time | Overall Alignment Rate |
|-----|--------------|------------------------------|-------------------------------------|-------------------------------|----------------------------|-----------------------|
| 1   | 50,831,351   | 16,973,703                   | 8,204,327                          | 25,653,321                    | 510,835                    | 73.58%                |
| 2   | 51,385,213   | 12,840,566                   | 5,019,290                          | 33,525,357                    | 327,247                    | 79.12%                |
| 3   | 49,828,147   | 15,597,984                   | 7,821,811                          | 26,408,352                    | 440,009                    | 74.77%                |
| 4   | 49,976,568   | 10,802,392                   | 5,426,523                          | 33,747,653                    | 349,808                    | 82.63%                |
| 5   | 48,908,730   | 20,325,695                   | 9,895,038                          | 18,687,997                    | 609,205                    | 67.15%                |


## _Porites evermanni_

| Run | Paired Reads | Aligned Concordantly 0 Times | Aligned Concordantly Exactly 1 Time | Aligned Concordantly >1 Times | Aligned Discordantly 1 Time | Overall Alignment Rate |
|-----|--------------|------------------------------|-------------------------------------|-------------------------------|----------------------------|-----------------------|
| 1   | 50,831,351   | 11,712,763                   | 30,209,320                          | 8,909,268                     | 2,145,413                  | 86.40%                |
| 2   | 51,385,213   | 11,914,068                   | 28,680,866                          | 10,790,279                    | 2,306,668                  | 87.10%                |
| 3   | 49,828,147   | 11,096,891                   | 29,566,240                          | 9,165,016                     | 1,929,423                  | 85.96%                |
| 4   | 49,976,568   | 10,514,435                   | 27,932,980                          | 11,529,153                    | 2,532,283                  | 89.07%                |
| 5   | 48,908,730   | 12,403,366                   | 28,498,490                          | 8,006,874                     | 2,215,956                  | 84.18%                |




