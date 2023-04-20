---
layout: post
title: LncRNA Discovery Workflow
date: '2023-04-20'
categories: lncRNA
tags: analysis pipeline
---

This workflow completes lncRNA discovery for a single sample (C17) of data from [Danielle Becker's nutrient enrichment study](https://github.com/hputnam/Becker_E5).

* [_P. ver_ genome assembly (GFF3 and fasta)](http://pver.reefgenomics.org/download/)
* [_P. ver_ genome assembly GTF](https://gannet.fish.washington.edu/Atumefaciens/20230127-pver-gff_to_gtf/)
* [RNA-seq data for sample C17](https://gannet.fish.washington.edu/Atumefaciens/hputnam-Becker_E5/Becker_RNASeq/data/raw/)

Build an index for a reference genome using the HISAT2 software. The command takes the path to the HISAT2 executable, -f option specifies that the input is a fasta file, and the path to the fasta file to be indexed...

```{bash hisat}
/home/shared/hisat2-2.2.1/hisat2-build \
-f ../data/Pver_genome_assembly_v1.0.fasta \
../output/Pver_genome_assembly_v1.0-valid.index
```

Align paired-end RNA-seq reads to a reference genome using the HISAT2 software. The command takes the path to the HISAT2 executable, -x option specifies the path to the index file for the reference genome, -p option specifies the number of threads to use, -1 and -2 options specify the paths to the first and second read files respectively, and -S option specifies the output file in SAM format...

```
/home/shared/hisat2-2.2.1/hisat2 \
-x ../output/Pver_genome_assembly_v1.0-valid.index \
-p 8 \
-1 ../data/C17_R1_001.fastq.gz \
-2 ../data/C17_R2_001.fastq.gz \
-S ../output/C17-valid.sam
```

"View" converts a SAM file to a BAM file. The command takes the path to the SAM file as input using the -bS option, and outputs to stdout. The output is piped to the next command. "Sort" then sorts the input BAM file and outputs a sorted BAM file. The command takes the path to the sorted output file, specified with the -o option, and sorts the input BAM file which is received through the pipeline...

```
/home/shared/samtools-1.12/samtools \
view -bS ../output/C17-valid.sam | \
/home/shared/samtools-1.12/samtools sort \
-o ../output/C17-valid_sorted.bam
```

Assemble and quantify transcript expression using RNA-seq alignment data. The command takes the path to the StringTie executable, -p option specifies the number of threads to use, -G option specifies the path to the reference annotation file in GTF format, -o option specifies the path to the output GTF file, and the path to the input aligned BAM file...

```
/home/shared/stringtie-2.2.1.Linux_x86_64/stringtie \
-p 8 \
-G ../data/Pver_genome_assembly_v1.0-valid.gtf \
-o ../output/Pver-valid.gtf \
../output/C17-valid_sorted.bam
```

Compare and merge multiple GTF files, and generate a summary report. The command takes the path to the gffcompare executable, -r option specifies the path to the reference annotation file in GTF format, -G option enables merging of compatible transcripts, -o option specifies the output file prefix, and the paths to the input GTF files.

```
/home/shared/gffcompare-0.12.6.Linux_x86_64/gffcompare \
-r ../data/Pver_genome_assembly_v1.0-valid.gtf \
-G \
-o ../output/C17 \
../output/Pver-valid.gtf \
```

* ` awk '$3 == "transcript" && $1 !~ /^#/ {print}' ../output/C17.annotated.gtf: selects only lines that correspond to transcript features and do not start with a "#" (i.e., comments) in the input GTF file../output/GFF_C17.annotated.gtf.
* grep 'class_code "u"': searches for lines that contain the string class_code "u" in the output of the previous awk command. The u class code indicates transcripts that are not overlapping any reference transcript on the same strand.
* awk '$5 - $4 > 199 {print}': selects only lines where the difference between the fifth column (end position) and the fourth column (start position) is greater than 199 nucleotides. This filters out transcripts that are shorter than 200 nucleotides.
* > ../output/novel_lncRNA_candidates.gtf: redirects the filtered output to a new GTF file ../output/novel_lncRNA_candidates-valid.gtf.

```
awk '$3 == "transcript" && $1 !~ /^#/ {print}' ../output/C17.annotated.gtf | grep 'class_code "u"' | awk '$5 - $4 > 199 {print}' > ../output/novel_lncRNA_candidates-valid.gtf
```

Bedtools to get fasta. Use the -fi option to specify the input reference genome FASTA file, the -bed option to specify the input GTF file, the -fo option to specify the output FASTA file, the -name option to use the transcript IDs as the names of the output sequences, and the -split option to split exons into separate lines.

```
/home/shared/bedtools2/bin/bedtools \
getfasta -fi ../data/Pver_genome_assembly_v1.0.fasta -bed ../output/novel_lncRNA_candidates-valid.gtf -fo ../output/novel_lncRNA_candidates-valid.fasta -name -split
```

Used the web-based version of CPC2 tool rather than completing the Python installation. Copy and paste fasta into the web interface [here](http://cpc2.gao-lab.org/) to generate text file that scores coding potential. Download and use this file in the following steps. You may need to adjust the next code block depending on the placement of your columns. 

Use awk to extract the transcript IDs from the cpc2results-last.txt file that have a score less than 0 and saves them to noncoding_transcripts_ids.txt. It uses grep to search for the transcript IDs in noncoding_transcripts_ids.txt within novel_lncRNA_candidates-valid.gtf and saves the output to final_lncRNAs.gtf.

```
awk '$NF < 0 {print $1}' ../output/cpc2results-last.txt > ../output/noncoding_transcripts_ids.txt
grep -Fwf ../output/noncoding_transcripts_ids.txt ../output/novel_lncRNA_candidates-valid.fasta > ../output/final_lncRNAs.gtf
```

The restult is a GTF file listing the IDs of lncRNAs present in you RNA-seq data.