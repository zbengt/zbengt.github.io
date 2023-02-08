---
layout: post
title: lncRNA in Express Compare Data
date: '2022-06-07'
categories: lncRNA
tags: workflow
---

Identifying long non-coding RNA in RNA-seq data from the Putnam Lab [Express-Compare repo](https://github.com/hputnam/Express_Compare).

## Species
* Montipora capitata
* Pocillopora acuta

## Reference Genomes
* [Montipora capitata](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_006542545.1/)
* [Pocillopora damicornis](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_003704095.1/), used for P. acuta

## RNA-seq Data hosted on NCBI
![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/lncRNA-EC_SRcodeTable.png?raw=true)

## General Workflow
![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/LncPipeWorkflow.png)

## Progress
* File download fixed with [SRA Toolkit configuration](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration) (run in Mac Terminal)
* FastQC run
* Run through HISAT for alignment