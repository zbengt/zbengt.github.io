---
layout: post
title: 02 - Differential Expression - CEABiGR mRNA
date: '2024-01-17'
categories: mRNA results status
tags: CEABiGR
---

Differential expression runs with DESeq2 for coding sequences. [Link to code](https://github.com/zbengt/oyster-lnc/blob/main/code/11-mRNA-counts-DESeq2.Rmd).

This wasn't the most interesting story. There are many DEGs, but appears to be mostly sex based. Super small number of differentially expressed DEGs when males and females are run separately. Continues to fit the story that most differences are sex-based. Can look at annotation to get a better idea of what the DEGs found actually are.

## Combined Run

### Volcano Plot
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-combined-volcano.png?raw=true)

* 36242 significantly differentially expressed coding sequences

### PCA - treatment
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-combined-PCA-treatment.png?raw=true)

### PCA - sex
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-combined-PCA-sex.png?raw=true)

### PCA - treatment - labelled
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-combined-PCA-treatment-labelled.png?raw=true)

### PCA - sex - labelled
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-combined-PCA-sex-labelled.png?raw=true)

### Top 50 Heatmap
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-combined-top50-heatmap.png?raw=true)

* Clustered by sex

## Female Run

### Volcano Plot
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-female-volcano.png?raw=true)

* 25 significantly differentially expressed lncRNA

### PCA - treatment
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-female-PCA.png?raw=true)

### PCA - treatment - labelled
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-female-PCA-labelled.png?raw=true)

### Top 50 Heatmap
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-female-top50-heatmap.png?raw=true)

* Clustered by treatment group

## Male Run

### Volcano Plot
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-male-volcano.png?raw=true)

* 107 significantly differentially expressed lncRNAs

### PCA - treatment
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-male-PCA.png?raw=true)

### PCA - treatment - labelled
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-male-PCA-labelled.png?raw=true)

### Top 50 Heatmap
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-male-heatmap.png?raw=true)

* Clustered by treatment group except for one, view below:

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/02-DE-male-heatmap-treatment.png?raw=true)
