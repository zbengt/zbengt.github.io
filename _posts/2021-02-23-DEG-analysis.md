---
layout: post
title:  Determining Taxa of DEGs
date: '2021-02-23'
categories: hematodinium
tags: overview, taxa
---

## Creating .fasta file of differentially-expressed transcripts

We previously created tables of differentially-expressed transcripts for all pairwise comparisons. We then decided to focus on two particularly relevant comparisons - [Ambient Day 2 vs. Elevated Day 2 (individual libraries only)](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/amb2_vs_elev2_indiv/DEGlist_wcols.txt) and [Elevated Day 0 vs. Elevated Day 2 (individual libraries only)](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/graphs/DESeq2_output/cbai_transcriptomev2.0/elev0_vs_elev2_indiv/DEGlist_wcols.txt). Both contrast samples taken at ambient temperatures vs. samples taken at elevated temperatures. 

I then used these tables to create a .fasta file containing all differentially-expressed transcripts by using [this script](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/2_2_DEG_blast.ipynb). I then used blastn to compare these .fasta files to a database containing all Alveolata nucleotide sequences. The output is [here for Ambient Day 2 vs. Elevated Day 2](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/BLASTs/alveolata_publicseqs/cbai_v2.0_amb2_vs_elev2_DEGs.tab) and [here for Elevated Day 0 vs. Elevated Day 2](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/output/BLASTs/alveolata_publicseqs/cbai_v2.0_elev0_vs_elev2_DEGs.tab). The database of all Alveolata nucleotide sequences was obtained from the NCBI Taxomony Browser at https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

I initially used megablast (the default) before changing to the more thorough blastn.

At this stage, I have a few next steps:
- Decide what our bar for evalue (expect value) should be to definitively say that a sequence is _Alveolata_
- Repeat this process, but with a .fasta file containing all Arthropoda sequences. Comparing the total number of matched sequences from both BLASTns should give us a good idea of how good this process is at determining taxa.
- Use BLAST with the taxonomy filter as another method of determining taxa

I also uncovered a potentially huge shortcut! Sam previously separated cbai_transcriptome_v2.0 (which is the one I'm currently using, and is unfiltered by taxa) into _Alveolata_ and non-_Alveolata_. Link [here](https://robertslab.github.io/sams-notebook/2020/06/05/Sequence-Extractions-C.bairdi-Transcriptomes-v2.0-and-v3.0-Excluding-Alveolata-with-MEGAN6-on-Swoose.html). I can do the following:
- Download the cbai_transcriptome_v2.1 and hemat_transcriptome_v2.1 .fasta files and see if they sum up to cbai_transcriptome_v2.0
- If yes, merge 2.1 with 2.0 and annotate non-_Alveolata_ (and potentially _Alveolata_?) by joining

Essentially, this should produce an annotated version of transcriptome 2.0, which will be super useful in future analyses!

