---
layout: post
title:  Ditching full genomes, using hematodinium transcriptomev1.6
date: '2021-03-16'
categories: hematodinium
tags: hemat_transcriptomev1.6, kallisto, DESeq2, GO-MWU
---

## Using NCBI Genomes

Here were my planned steps last time:

1. BLASTx of transcriptome 2.0 against the protein sequences from the _C. opilio_ genome downloaded from NCBI

2. BLASTn of transcriptome 2.0 against the _Amoebophrya sp._ genome downloaded from NCBI

3. Use those two BLAST results to determine which genes are _Hematodinium_ in origin and which are _C. bairdi_ for all of Transcriptome 2.0

I completed all steps (plus a few tBLASTx of transcriptome 2.0 against _Amoebophrya_ genome), but ran into a major problem. Specifically, not enough of our sequences are finding matches - after combining all matches for the _C. opilio_ and _Amoebophrya sp._ genomes, only around half were matched, even after dropping the e-value bar all the way down to 10^-2 (and that's with double-counting all sequences that matched to both genomes!). Link to all scripts and the specific numbers of matches available [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/2_4_ncbi_genome_blasts.ipynb)

## So what now?

Well, we're still trying to diagnose what sequences are _Hematodinium_ and what sequences are _C. bairdi_. However, the best route to accomplish that may be redoing most of our previous steps. Specifically, we'll take our libraries, align them to a transcriptome using kallisto, analyze differential expression with DESeq2, examine our DEGs in detail, and then look at pathway enrichment with GO-MWU. However, unlike last time, where we aligned our libraries to a transcriptome containing a mix of _Hematodinium_ and _C. bairdi_, this time, we'll align to a transcriptome containing only sequences that are presumably _Hematodinium_. This transcriptome - hemat_transcriptomev1.6 - was created by Sam White, with details and a link to download is available [here](https://robertslab.github.io/sams-notebook/2021/03/08/Transcriptome-Assembly-Hematodinium-Transcriptomes-v1.6-and-v1.7-with-Trinity-on-Mox.html)

The scripts will be more or less identical to our past ones. For notation, we're deeming this Analysis 3. Rundown of past analyses:

**Note**: Scripts are annotated analysis/script number. Analysis 1, Script 5 = 15_script_title.

Analysis 1: Scripts 10-19. Putting mixed _C. bairdi_ / _Hematodinium_ transcriptome, specifically [cbai_transcriptomev2.0](https://robertslab.github.io/sams-notebook/2020/05/02/Transcriptome-Assembly-C.bairdi-All-RNAseq-Data-Without-Taxonomic-Filters-with-Trinity-on-Mox.html), available via [Genomic Resources](https://robertslab.github.io/resources/Genomic-Resources/), through the above pipeline (kallisto -> DESeq2 -> GO-MWU)

Analysis 2: Scripts 20-29: Misc analyses exploring output from Analysis 1, trying to determine what sequences are _Hematodinium_ and what sequences are _C. bairdi_.

Analysis 30-39: Putting a _Hematodinium_ only transcriptome, AKA hemat_transcriptomev1.6, through the same pipeline as Analysis 1.

## Analysis 3: Looking at _hemat_transcriptomev1.6_

As I mentioned, the plan is largely the same as that of Analysis 1. Put libraries through kallisto, then DESeq2, examine DEGs, and then finally analyze with GO-MWU. After that, we'll extend the analysis further by analyzing with [WGCNA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559). 

I've already progressed all the way through GO-MWU! Scripts available as follows:

[Creating BLASTx index for hemat_transcriptomev1.6](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/3_0_hemat1.6_indexcreation.ipynb)

[Running kallisto](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/3_1_download_libraries_run_kallisto.ipynb)

[Running DESeq2](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/3_2_kallisto_to_deseq_to_accessionIDs.Rmd)

[Obtaining GO terms](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/3_3_uniprot_to_GO.Rmd)

[Preparing one of the inputs for GO-MWU](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/3_4_eliminate_duplicates.ipynb)

[Preparing the other input for GO-MWU](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/3_5_GO-MWU_prep.Rmd)

[Running GO-MWU](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/3_6_running_GO-MWU/3_6_running_GO-MWU.R)

Two key comparisons were made - Elevated Day 0 vs. Elevated Day 2, individual libraries only (remember, Elevated Day 0 samples were taken prior to exposure to elevated-temperature water), and Ambient Day 2 vs. Elevated Day 2, individual libraries only. All 3 groups (Elev. Day 0, Elev. Day 2, Amb. Day 2) have 3 libraries.

## Results of Analysis 3

## [DESeq2 results](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/graphs/DESeq2_output/hemat_transcriptomev1.6)

The first thing to note is that our number of DEGs is quite low. For Elev. Day 0 vs. Elev. Day 2, we had [4 total differentially-expressed transcripts](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/hemat_transcriptomev1.6/elev0_vs_elev2_indiv/DEGlist_wcols.txt), of which [3 matched to an accession ID](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/accession_n_GOids/hemat_transcriptomev1.6/DEG_IDs/elev0_vs_elev2_indiv_DEG_IDs.txt). For Ambient Day 2 vs. Elev. Day 2, we had [only one differentially-expressed transcript](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/hemat_transcriptomev1.6/amb2_vs_elev2_indiv/DEGlist_wcols.txt), and [it matched to an accession ID](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/accession_n_GOids/hemat_transcriptomev1.6/DEG_IDs/amb2_vs_elev2_indiv_DEG_IDs.txt)

### DEG Rundown

#### Amb. Day 2 vs. Elev. Day 2

[**O99262**](https://www.uniprot.org/uniprot/O99252): This gene matches to the mitochondrial gene **cytochrome c oxidase subunit I, AKA COI**. This is...really interesting! I've got a bit of a history with COI - my paper on the European green crab found a linkage between thermal tolerance and COI genotype. This is obviously a bit different, as we're examining expression rather than genotype, and it's a lot more intuitive that mitochondrial activity could be linked to temperature, but still...neat! 

A caveat: I'm not completely sure this is real. Our p-value is 5.5x10^-7, which is nice, but our p-adj is just 0.0023. And looking at other ambient-temperature libraries, 5 out of 14 had a expression level equal or lower to that of the most-expressed elevated-temperature library. So just something to keep in mind. 

#### Elev. Day 0 vs. Elev. Day 2

[**Q06056**](https://www.uniprot.org/uniprot/Q06056): Another mitochondrial gene! This one matches to ATP synthase F(0) complex subunit C2. P-values look a bit better for this one (pval = 9.13x10^-8, padj = 4.0 x 10^-5). Interesting that we're seeing multiple mitochondrial genes here - definitely worth a deeper dive later on.

[**O23717**](https://www.uniprot.org/uniprot/O23717). The first of two proteosome subunits- this one is proteosome subunit beta type-5-a. Looks like it cleaves peptides, and some quick Google Scholar work found that in humans, it's involved in immunosuppression and may be oncogenic. Again, p-values are solid (pval = 6.1x10^-8, padj = 4.3x10^-5).

[**P85200**](https://www.uniprot.org/uniprot/P85200). The second of two proteosome subunits - this one is proteosome regulatory subunit 6b homolog. It's involved in the ATP-dependent degradation of ubiquinated proteins (whatever that means). Again, good pvals (pval = 7.0x10^-8, padj=4.3x10^-5)

## [GO-MWU Output](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/output/GO-MWU_output/hemat_transcriptomev1.6)

We weren't able to make graphs for elev. day 0 vs. elev. day 2, as there just weren't any enriched pathways. However, we did see two enriched pathways for Ambient Day 2 vs. Elevated Day 2 [NOTE: no longer the case as of 2021-11-08] specifically in organelle localization and cellular localization. This could certainly be due to different strains in the two different groups of crab, but is definitely worth looking into more.



## Main takeaways: 

I really need to do more investigation here, but it's quite interesting how much overlap there is in the function of our DEGs, with two being mitochondrial and two being proteosome subunits. Also interesting that there's no overlap in DEGs between our two comparisons! Anyway, I'll do some more investigation of this.



#### Next steps:

A few next steps to take:

1) Investigate function of DEGs listed above

2) Further examination of GO-MWU output

3) Perform same steps on a crab-specific transcriptome (currently under construction by Sam White)

4) Learn how to run WGCNA, then run on hemat_transcriptomev1.6 libraries (and potentially crab-specific libraries as well)

