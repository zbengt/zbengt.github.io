---
layout: post
title: 00 - Size Distribution - CEABiGR lncRNA
date: '2024-01-11'
categories: lncRNA results status
tags: CEABiGR
---

Size distribution figures of CEAGBiGR lncRNAs. [Link to code](https://github.com/zbengt/oyster-lnc/blob/main/code/09-size-distribution-and-location.Rmd).

## Full size distribution

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR_full_size_distribution.png?raw=true)

## Binned size distribution

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR_binned_size%20distribution.png?raw=true)

This is suspicious. The vast majority are greater than 1000nt and I don't buy it. The FASTA also includes sequences below a 200nt threshold. Seems like this FASTA from NCBI might just be  uncharacterized regions as opposed to lncRNAs. This would make sense since most of the RNAcentral hits from this blast run are "uncharacterized regions" from the same annotation associated with this FASTA.





