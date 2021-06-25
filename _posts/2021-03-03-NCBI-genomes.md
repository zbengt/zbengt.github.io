---
layout: post
title:  Using NCBI Genomes
date: '2021-03-03'
categories: hematodinium
tags: overview, taxa
---

## Transcriptome 2.1 update...

Last time I mentioned potentially using Sam's transcriptomes created by taking cbai_transcriptome_v2.0 and filtering _Alveolata_ and non-_Alveolata_ by taxa, thus producing [cbai_transcriptome_v2.1 and hemat_transcriptome_v2.1](https://robertslab.github.io/sams-notebook/2020/06/05/Sequence-Extractions-C.bairdi-Transcriptomes-v2.0-and-v3.0-Excluding-Alveolata-with-MEGAN6-on-Swoose.html). However, when I downloaded those two transcriptomes, they had ~230,000 and ~30,000 sequences respectively - far less than the ~1.4 million found in cbai_transcriptome_v2.0. Therefore, for the time being, we are neglecting this alternative route. Alright, on to my main update!

## Last time..

In my last post, I characterized the method of choice for obtaining taxa information. The plan was to do the following:

1. Download all Alveolata and Arthropoda nucleotide sequences from the NCBI database, available from the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)

2. Obtain a FASTA file of all DEGs for our two meaningful conditions - elev day 2 vs. amb. day 2, and elev day 0 vs. elev day 2 - by cross-referencing transcript IDs and transcriptome 2.0

3. BLASTn [those obtained FASTA files](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/output/BLASTs/input_seqs) twice - once against all Alveolata sequences and again against our Arthropoda sequences

4. Figure out a good e-value to set as our bar for taxa determination

5. Within each set of DEGs, see how many sequences appear as matches for both Alveolata and Arthropoda, and how many failed to match either

All five of the above steps were completed using [this script](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/2_2_DEG_blast.ipynb). 
 
[Here are my results for the alveolata BLAST](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/output/BLASTs/alveolata_publicseqs), and [for the arthropoda BLAST](https://github.com/afcoyle/hemat_bairdi_transcriptome/tree/main/output/BLASTs/arthropoda_publicseqs)

Here's what we got.


|                                                      | Amb. Day 2 vs. Elev. Day 2 | Elev. Day 0 vs. Elev. Day 2 |
|------------------------------------------------------|----------------------------|-----------------------------|
| Base Sequences                                       | 2069                       | 338                         |
| Alveolata (unfiltered)                               | 2391                       | 351                         |
| Arthropoda (unfiltered)                              | 2736                       | 349                         |
| Alveolata (eval <= 10^-4)                            | 569                        | 60                          |
| Arthropoda (eval <= 10^-4)                           | 1137                       | 181                         |
| Eval: Alveolata > Arthropoda (presumably Arthropoda) | 1175                       | 217                         |
| Eval: Alveolata < Arthropoda (presumably Alveolata)  | 741                        | 66                          |


Note: For the last two rows, duplicate transcript IDs were removed (some transcript IDs mapped to multiple genes within a single BLAST). The transcript ID with the higher e-value was retained

As a quick initial scan, this...is pretty good!! Definitely a good guide to which genes are likely _Hematodinium_ and which are _C. bairdi_. However, it has a downside - it exclusively examines differentially-expressed genes. After considering some whole-transcriptome alternatives, we made a discovery - as of January 2021, someone has uploaded a fairly complete _Chionoecetes opilio_ genome to the NCBI database - complete with a separate file of presumed protein sequences! We also located a full genome of a relatively-close species to _Hematodinium_ - _Amoebophrya sp._, a dinoflagellate that parasitizes other dinoflagellates. Therefore, our next steps are as follows:

1. BLASTx of transcriptome 2.0 against the protein sequences from the _C. opilio_ genome downloaded from NCBI

2. BLASTn of transcriptome 2.0 against the _Amoebophrya sp._ genome downloaded from NCBI

3. Use those two BLAST results to determine which genes are _Hematodinium_ in origin and which are _C. bairdi_ for all of Transcriptome 2.0

4. When time allows, BLAST all of Transcriptome 2.0 against the whole NCBI database with a taxonomy filter (run this on Mox).

5. Potential alternative step if needed: Take all _C. bairdi_ libraries of uninfected crab, and assemble a _C. bairdi_ transcriptome. Use that to determine which sequences are _Hematodinium_ and which are _C. bairdi_.

Steps 1-2 have already been completed - a script for those, including all data locations, is available [here](https://github.com/afcoyle/hemat_bairdi_transcriptome/blob/main/scripts/2_4_ncbi_genome_blasts.ipynb)! The slurm scripts are currently running on Mox. Once they're complete, time to move on to Step 3!

