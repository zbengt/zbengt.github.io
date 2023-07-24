---
layout: post
title: CEABiGR lncRNA 01 - DESeq2 Runs
date: '2023-07-24'
categories: bioinformatics CEABiGR
tags: CEABiGR
---

Ran DESeq2 for females only, males only, and combined not accounting for sex. Volcano plots, PCAs, and heatmaps of the top 30 most differentially expressed included for each.

### Number of lncRNAs with < 0.05 p-values
* Females: 1
* Males: 5
* Combined: 2

### Females - OA vs Control

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/Female-OA-volcanoplot.png?raw=true)

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/Female-OA-PCA.png?raw=true)

No obvious clustering as a whole.

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/Female-OA-Top30.png?raw=true)

Control grouped with control and exposed grouped with exposed for the top 30.


### Males - OA vs Control

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/Males-OA-volcano-plot.png?raw=true)

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/Males-OA-PCA.png?raw=true)

Clustering of most males.

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/Males-OA-Top30.png?raw=true)

Control grouped with control and exposed grouped with exposed for the top 30.


### Combined (Males and Females) - OA vs Control

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/Combined-OA-volcanoplot.png?raw=true)

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/Combined-OA-PCA.png?raw=true)

Clustering consistent for males.

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/Combined-OA-heatmap.png?raw=true)

Grouping based on sex.



