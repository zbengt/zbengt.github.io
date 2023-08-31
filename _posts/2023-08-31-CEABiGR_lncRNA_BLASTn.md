---
layout: post
title: CEABiGR lncRNA - RNAcentral BLASTn
date: '2023-08-31'
categories: bioinformatics CEABiGR
tags: CEABiGR
---

Started with a couple smaller datasets included in the umbrella of RNAcentral, but got almost no hits so scaled up. Blasted the entire RNAcentral database (over 3 million sequences) with the Crassostrea virginica lncRNA FASTA file.

- [Code](https://github.com/zbengt/oyster-lnc/blob/main/code/07-RNAcentral-database-comparison.Rmd)
- [RNAcentral link](https://rnacentral.org/)
- [Ensemble Metazoa link](https://metazoa.ensembl.org/index.html)

Key Points:
- This isn't really an annotation, just a database comparison to see if any of the lncRNAs are better described in other speacies
- Was able to get a hit for every lncRNA transcript, but this isn't surprising since the database includes uncharacterized RNAs from C. virginica in the Ensembl Metazoa database.
- Not sure trying to annotate this way does anything meaningful for us since the count matrix already uses the IDs associated with these catalogged but uncharacterized RNAs.

Summary stats of blast run:

- All 4750 lncRNA transcripts obtained 1 or greater hits
- [Results table](https://github.com/zbengt/oyster-lnc/blob/main/output/database_blasts/RNAcentral/rnacentral.tab) (hits are represented by RNAcentral IDs)

Distribution of Percentage Identity

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/RNAcentral-blast-percentage-identity.png?raw=true)

E-value Distribution

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/RNAcentral-blastn-e-value-distribution.png?raw=true)



