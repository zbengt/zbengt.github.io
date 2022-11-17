---
layout: post
title: Mussel Heat Stress Primer Design
date: '2022-08-12'
categories: mussels
tags: heat stress
---

### 08-02-2022
Designing heat stress specific primers for Mytilus trossulus/edulis…

Previous primers used in class lab were from this paper: [DNA Damage and Transcriptional Changes in the Gills of Mytilus galloprovincialis Exposed to Nanomolar Doses of Combined Metal Salts (Cd, Cu, Hg)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0054602#pone.0054602-Livak1)

Primers might not have been appropriate for heat stress since the study was directed more toward metal toxicity. The mussel species may have also been too dissimilar from those we used (primers were for galloprovincialis rather than trossulus or edulis).

<!---Currently looking at [this](https://cob.silverchair-cdn.com/cob/content_public/journal/jeb/213/20/10.1242_jeb.046094/3/3548.pdf?Expires=1662675798&Signature=yRo3Amp2FmeLdu6th0ysrMrBJO3ai3txBJOVHYyQOf41WC-082GbJnwM7gogrwiXX-S9WbqCzYV3NtWoVsghPmwnFQ0Crr4tzIHZKdlWp-4bqv5JxcMsl004fh9A2St3kfLDusO9ilkQV~1UjfCshuQMTwyptW0cLuWVywDa~LbBKfPZruEA~AXr0S9AZuT5-XCsKH-pabvUXERTaXoNcpnFNvtM3ewkOP3jY~uJhmcu1mWN8TlN-eC-Pzj00FPb-K9pzOB~OGSNKTEz1rqZM8pAD8Um5XID2PBhFrm3gjhAPELVvLBvHh6l8qIALNO-hHMLnT5MTDwEMaYfWVrrDg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) publication comparing differences in expression after acute heat stress between galloprovincialis and trossulus. Seems like some common players are the same between the two, like Hsp70s, Hsp90, Hsp60, DnaJ and chaperonin TCP1.-->

Hsp70s gene expression:

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel_HeatStress_Primers_Hsp70s.png)

Need to look a little closer at some of the figures to get species specific genes that saw up-regulation in trossulus and continue lit review.

### 08-05-2022 Update
qPCR runs from Dorothy's experiment showed successful amplification for all primers used in class except for ferritin. So I will likely re-rerun samples from class with the primers for MT10, MT20, small HSP24.1, SQSTM1, HSC70, HSP90, and GADD45 gamma.

Primers (will exclude ferritin):

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/FISH541_Mussel_HeatStress_Primers.png)

### 08-12-2022 Update
Melt curves from Dorothy's experiment (1 hour of heat stress) show good amplification of MT10, MT20, small HSP24.1, SQSTM1, HSC70, HSP90, and GADD45 gamma. Her findings show the most significant response in small HSP24.1. Dorothy's lab notebook can be found [here](https://genefish.wordpress.com/author/dorolar/).

**HSC70 & small HSP24**

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel_HeatStress_HSC70-sHSP24.png)

**HSP90 & SQSMT1**

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel_HeatStress_HSP90-SQSMT1.png)

**MT10 & MT20**

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel_HeatStress_MT10-MT20.png)

**GADD45 gamma**

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel_HeatStress_GADD45.png)

**Ferritin**

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel_HeatStress_Ferritin.png)

I will be running pPCR on existing cDNA from the class experiment using the small HSP24.1 primer.

### 08-18-2022 Update
We have 64 cDNA samples from experiment 2: replicates, A & B, for each individual (chronic heat of 18C, variable cycling heat 12C-18C, and control 12C). To run them with the HSP24 primer and for all subsequent analyses, I made 1:10 dilutions for each cDNA sample (40ul each, 4ul sample and 36ul DEPC H2O).

### 08-22-2022
Ran qPCR of A & B cDNA samples individually for HSP24. A & B are cDNA synthesis replicates (and A & B exist for each sample in Experiment 2).



### 08-31-2022 Update
Ran qPCR of pooled A & B cDNA samples for HSP24.

* [View data](https://docs.google.com/spreadsheets/d/1R94cfE804rsnnVULtMr3EdcrVhbc-2Zd813JzGFxB_8/edit?usp=sharing)
* [Full output folder](https://drive.google.com/drive/folders/1MGApIrCIWjuW92qNg4fLl7E1hM4wGON4?usp=sharing)

Melt Curve:

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel%20qPCR/08-31-2022/melt_curve.png)


Amplification:

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel%20qPCR/08-31-2022/amplification.png)

Melt curve shows potential primer-dimer and amplification in one of the NTCs. Standard deviations between replicates are suspiciously high. 

### 09-08-2022 Update
Ran qPCR of pooled A & B cDNA samples for HSP90.

* [View data](https://docs.google.com/spreadsheets/d/16AAFDXpL_6WE68svOdyC5H9ZU2gKQ5so9CqVwPqklgY/edit?usp=sharing)
* [Full output folder](https://drive.google.com/drive/folders/1PrBuZzbUfCa5KAIEYYW_ox13woH1kHVs?usp=sharing)

Melt Curve:

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel%20qPCR/09-08-2022/melt_curve.png)

Amplification:

![image](https://raw.githubusercontent.com/zbengt/zbengt.github.io/master/assets/img/Mussel%20qPCR/09-08-2022/amplification.png)

Melt curves look good, no amplification in NTCs, and standard deviations between replicates look alright.


**to be continued…**