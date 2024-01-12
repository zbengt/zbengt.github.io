---
layout: post
title: 01 - Differential Expression - CEABiGR lncRNA
date: '2024-01-12'
categories: lncRNA results status
tags: CEABiGR
---

Differential expression runs with DESeq2. [Link to code](https://github.com/zbengt/oyster-lnc/blob/main/code/10-count-matrices-DESeq2-final.Rmd).

I think what ended up being the most interesting is the combined run accouning for sex and treatment, All of the observed differences in expression appear to be sex-based, but this in and of itself is interesting and it gives us the IDs of 2370 differentially expressed lncRNAs between female and male sample groups.

## Combined Run

### Volcano Plot
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-combined_volcanoplot.png?raw=true)

* 2370 significantly differentially expressed lncRNA

### PCA - treatment
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-combined_PCA-treatment.png?raw=true)

### PCA - sex
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-combined_PCA-sex.png?raw=true)

### PCA - treatment - labelled
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-combined_PCA_treatment-labelled.png?raw=true)

### PCA - sex - labelled
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-combined_PCA_sex-labelled.png?raw=true)

### Top 50 Heatmap
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-combined_heatmap50.png?raw=true)

* Clustered by sex

## Female Run

### Volcano Plot
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-female_volcanoplot.png?raw=true)

* 1 significantly differentially expressed lncRNA

### PCA - treatment
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-female_PCA.png?raw=true)

### PCA - treatment - labelled
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-female_PCA-labelled.png?raw=true)

### Top 50 Heatmap
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-female_heatmap50.png?raw=true)

* Clustered by treatment group

## Male Run

### Volcano Plot
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-male_volcanoplot.png?raw=true)

* 5 significantly differentially expressed lncRNAs

### PCA - treatment
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-male_PCA.png?raw=true)

### PCA - treatment - labelled
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-male_PCA-labelled.png?raw=true)

### Top 50 Heatmap
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/CEABiGR-male_heatmap50.png?raw=true)

* Clustered by treatment group
