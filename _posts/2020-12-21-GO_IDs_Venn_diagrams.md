---
layout: post
title:  GO terms and Venn Diagrams
date: '2020-12-21'
categories: hematodinium
tags: hematodinium, GO IDs
---

**Obtaining GO IDs**

My analysis prior to this used DESeq2 to obtain a series of tables of genes with significantly-different expressions for each of my four comparisons.

I wrote a [script in R](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/1_2_kallisto_to_deseq_to_accessionIDs.Rmd) that took in each DESeq2 output file, matched it against [this annotated transcriptome](https://gannet.fish.washington.edu/Atumefaciens/20200519_cbai_diamond_blastx_transcriptome_v3.0/20200518.C_bairdi.Trinity.blastx.outfmt6), and produced a  newline-separated file of UniProt accessions. I then took each newline-separated file and put it into [this Bash script from Sam](https://github.com/RobertsLab/code/blob/master/script-box/uniprot2go.sh), which retrieved the Gene Ontology terms for each Accession ID. 

This is progressing fairly well!
Next steps are as follows:
- Get the GOslim terms using the GSEAbase package
- Perform gene enrichment analysis with TopGO or GO-MWU
- When time allows, go back and annotate the transcriptome myself to understand it more fully


**Producing Venn Diagrams**

In order to determine the overlap between DEGs for my different comparisons, I produced a few Venn diagrams. Since a quick analysis found practically no overlap between the time comparison (Day 0 vs Day 17) and any of my temperature comparisons, the Venn diagram solely includes our three temperature comparisons (Day 0/2: Elevated vs Ambient, Ambient vs. Low, Elevated vs. Low).

[Overlap in transcript IDs](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/DEG_venn_diagrams/TranscriptID_DEGs.png)

Some stats for our DEGs:
- 2166 unique transcript IDs
- 2919 total transcript IDs
- 74.2% unique transcript IDs

[Overlap in accession IDs](https://github.com/afcoyle/hemat_bairdii_transcriptome/blob/main/graphs/DEG_venn_diagrams/AccessionIDs_DEGs.png)

Some stats for the accession IDs of our DEGs:
- 633 unique accession IDs
- 1061 total accession IDs
- 59.6% unique accession IDs


