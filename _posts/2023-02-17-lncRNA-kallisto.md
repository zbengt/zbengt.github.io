---
layout: post
title: CEABiGR Kallisto Run for lncRNA Count Matrix
date: '2023-02-17'
categories: genomics expression
tags: analysis lncRNA
---

**Completed**

* Run trimmed RNA-seq data in kallisto using the lncRNA fasta file as reference

* Code can be found [here](https://github.com/sr320/ceabigr/blob/main/code/13-lncRNA-kallisto.Rmd)

* Follow up at pub-a-thon to merge abundance.tsv files for each library.

* Create new repo for this using the CEABiGR code as a reference


Using kallisto to get count data for lncRNA...

Need fasta and fq

Grab .fq files

https://gannet.fish.washington.edu/Atumefaciens/20220224_cvir_gonad_RNAseq_fastp_trimming/S12M_R1.fastp-trim.20bp-5prime.20220224.fq.gz



```{bash}
wget -r \
--no-check-certificate \
--quiet \
--no-directories --no-parent \
-P /home/shared/8TB_HDD_02/zbengt/data/oyster-lnc/kallisto \
-A *fq.gz \
https://gannet.fish.washington.edu/Atumefaciens/20220224_cvir_gonad_RNAseq_fastp_trimming/
```


Downloading lncRNA fasta

```{bash}
wget --no-check-certificate \
-P /home/shared/8TB_HDD_02/zbengt/data/oyster-lnc/01-lncRNA-kallisto/fasta \
https://gannet.fish.washington.edu/Atumefaciens/20220217-cvir-lncRNA_subsetting/GCF_002022765.2_C_virginica-3.0_lncRNA.fa
```


Creating an index
```{bash}
/home/shared/kallisto/kallisto \
index -i /home/shared/8TB_HDD_02/zbengt/data/oyster-lnc/01-lncRNA-kallisto/index/GCF_002022765.2_C_virginica-3.0_lncRNA.index \
/home/shared/8TB_HDD_02/zbengt/data/oyster-lnc/01-lncRNA-kallisto/fasta/GCF_002022765.2_C_virginica-3.0_lncRNA.fa
```

```{bash}
find /home/shared/8TB_HDD_02/zbengt/data/oyster-lnc/01-lncRNA-kallisto/fq/*_R1.fastp* \
| xargs basename -s _R1.fastp-trim.20bp-5prime.20220224.fq.gz | xargs -I{} /home/shared/kallisto/kallisto \
quant -i /home/shared/8TB_HDD_02/zbengt/data/oyster-lnc/01-lncRNA-kallisto/index/GCF_002022765.2_C_virginica-3.0_lncRNA.index \
-o /home/shared/8TB_HDD_02/zbengt/github/oyster-lnc/output/01-lncRNA-kallisto/{} \
-t 40 \
/home/shared/8TB_HDD_02/zbengt/data/oyster-lnc/01-lncRNA-kallisto/fq/{}_R1.fastp-trim.20bp-5prime.20220224.fq.gz /home/shared/8TB_HDD_02/zbengt/data/oyster-lnc/01-lncRNA-kallisto/fq/{}_R2.fastp-trim.20bp-5prime.20220224.fq.gz \
2>&1 | tee ../output/01-lncRNA-kallisto/{}.out
```
