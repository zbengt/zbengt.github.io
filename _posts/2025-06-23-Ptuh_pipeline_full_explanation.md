---
title: "06.23.2025 Ptuh lncRNA Pipeline Updates: GFFcompare, duplication fixes, and GTF formatting"
date: 2025-06-23
layout: post
categories: lncRNA bioinformatics e5
---

TL;DR - Fixing the system environment settings for gffcompare and fixing formatting issues with out GTFs appears to have fixed issues with pipeline. The number of putative lncRNAs is now much lower (more similar to deep-dive descriptive numbers) and no duplicates are detected in bed, gtf, fasta, and count files.

# Ptuh lncRNA Pipeline: Step-by-Step Explanations

Below is a detailed breakdown of each major section of the Ptuh lncRNA pipeline, describing what each step does and why it’s important.

---

## 1. Variables and Environment Setup

```r
# Global R options
knitr::opts_chunk$set(echo = TRUE)

# Define key paths and tool directories
DATA_DIR     <- "../data/01.6-Ptuh-lncRNA-pipeline"
OUTPUT_DIR   <- "../output/01.6-Ptuh-lncRNA-pipeline"
THREADS      <- "24"
FASTQ_SOURCE <- "…/P_meandrina/trimmed/"
FASTQ_SUFFIX <- "fastq.gz"
GENOME_SOURCE<- "…/Pocillopora_meandrina_HIv1.assembly.fasta"
GTF_SOURCE   <- "…/Pocillopora_meandrina_HIv1.genes-validated.gtf"
GFF_SOURCE   <- "…/Pocillopora_meandrina_HIv1.genes-validated.gff3"
GFFPATTERN   <- 'class_code "u"|class_code "x"|class_code "o"|class_code "i"'

HISAT2_DIR     <- "/home/shared/hisat2-2.2.1/"
SAMTOOLS_DIR   <- "/home/shared/samtools-1.12/"
STRINGTIE_DIR  <- "/home/shared/stringtie-2.2.1.Linux_x86_64"
GFFCOMPARE_DIR <- "/home/shared/gffcompare-0.12.6.Linux_x86_64"
BEDTOOLS_DIR   <- "/home/shared/bedtools2/bin"
CPC2_DIR       <- "/home/shared/CPC2_standalone-1.0.1/bin"
CONDA_PATH     <- "/opt/anaconda/anaconda3/bin"

GENOME_FASTA <- file.path(DATA_DIR, "genome.fasta")
GENOME_GTF   <- file.path(DATA_DIR, "genome.gtf")
GENOME_GFF   <- file.path(DATA_DIR, "genome.gff")
FASTQ_DIR    <- file.path(DATA_DIR, "fastq")
GENOME_INDEX <- file.path(OUTPUT_DIR, "genome.index")

Sys.setenv(
  THREADS=THREADS, DATA_DIR=DATA_DIR, OUTPUT_DIR=OUTPUT_DIR,
  FASTQ_SOURCE=FASTQ_SOURCE, FASTQ_SUFFIX=FASTQ_SUFFIX,
  GENOME_SOURCE=GENOME_SOURCE, GTF_SOURCE=GTF_SOURCE,
  GFF_SOURCE=GFF_SOURCE, GFFPATTERN=GFFPATTERN,
  HISAT2_DIR=HISAT2_DIR, SAMTOOLS_DIR=SAMTOOLS_DIR,
  STRINGTIE_DIR=STRINGTIE_DIR, GFFCOMPARE_DIR=GFFCOMPARE_DIR,
  BEDTOOLS_DIR=BEDTOOLS_DIR, CPC2_DIR=CPC2_DIR,
  CONDA_PATH=CONDA_PATH, GENOME_FASTA=GENOME_FASTA,
  GENOME_GTF=GENOME_GTF, GENOME_GFF=GENOME_GFF,
  FASTQ_DIR=FASTQ_DIR, GENOME_INDEX=GENOME_INDEX
)
```

> **Why?**  
> Defining all variables and paths at the top centralizes configuration, making the pipeline reproducible and easy to update.

---

## 2. Directory Creation

```bash
mkdir -p "${DATA_DIR}" "${OUTPUT_DIR}"
```

- **Purpose:** Create the required data and output folders in one go.

---

## 3. Download Genome and Reads

```bash
wget -nv -r --no-directories --no-parent   -P "${FASTQ_DIR}" -A "*${FASTQ_SUFFIX}" "${FASTQ_SOURCE}"

curl -o "${GENOME_FASTA}" "${GENOME_SOURCE}"
curl -o "${GENOME_GTF}"   "${GTF_SOURCE}"
curl -o "${GENOME_GFF}"   "${GFF_SOURCE}"
```

- **`wget`**: Retrieves all trimmed FASTQ files.  
- **`curl`**: Downloads the reference genome, GTF, and GFF3.

---

## 4. Integrity Checks

```bash
head -1 "${GENOME_FASTA}"
head -2 "${GENOME_GFF}"
head -1 "${GENOME_GTF}"
```

- **Goal:** Quickly verify that the downloads are actual sequence/annotation files rather than HTML errors.

---

## 5. HISAT2 Indexing & Alignment

```bash
"${HISAT2_DIR}/hisat2_extract_exons.py" "${GENOME_GTF}" > exon.txt
"${HISAT2_DIR}/hisat2_extract_splice_sites.py" "${GENOME_GTF}" > splice_sites.txt

"${HISAT2_DIR}/hisat2-build" -p "${THREADS}"   "${GENOME_FASTA}" "${GENOME_INDEX}"   --exon exon.txt --ss splice_sites.txt
```

```bash
for r2 in "${FASTQ_DIR}"/*_R2_*."${FASTQ_SUFFIX}"; do
  sample="${r2##*/}"
  sample="${sample%%_R2_*}"
  r1="${r2/_R2_/_R1_}"
  "${HISAT2_DIR}/hisat2" -x "${GENOME_INDEX}" -p "${THREADS}"     -1 "$r1" -2 "$r2" -S "${OUTPUT_DIR}/${sample}.sam"
done
```

- **Extract hints** (exons, splice sites) for better alignment.  
- **Index** the genome with those hints.  
- **Align** each paired sample to generate SAM files.

---

## 6. SAM → Sorted & Indexed BAM

```bash
for sam in "${OUTPUT_DIR}"/*.sam; do
  bam="${sam%.sam}.bam"
  sorted="${sam%.sam}.sorted.bam"
  samtools view -bS -@ "${THREADS}" "$sam" > "$bam"
  samtools sort -@ "${THREADS}" "$bam" -o "$sorted"
  samtools index -@ "${THREADS}" "$sorted"
done
```

- **Convert** SAM → BAM,  
- **Sort** BAM for downstream compatibility,  
- **Index** for rapid access.

---

## 7. Transcript Assembly with StringTie

```bash
find "${OUTPUT_DIR}" -name "*sorted.bam"   | xargs -n1 -I{} "${STRINGTIE_DIR}/stringtie"       -p "${THREADS}" -G "${GENOME_GFF}"       -o "{}.gtf" "{}"

"${STRINGTIE_DIR}/stringtie" --merge -G "${GENOME_GFF}"   -o "${OUTPUT_DIR}/stringtie_merged.gtf"   "${OUTPUT_DIR}"/*.gtf
```

- **Per-sample assembly** of transcript models.  
- **Merge** across samples into a unified annotation.

At this point there are 552268 transcripts.

---

## 8. GFFCompare & Candidate Selection

```bash
"${GFFCOMPARE_DIR}/gffcompare" -r "${GENOME_GFF}"   -o "${OUTPUT_DIR}/gffcompare_merged"   "${OUTPUT_DIR}/stringtie_merged.gtf"

awk '$3=="transcript" && !/^#/ && /'"$GFFPATTERN"'/'   gffcompare_merged.annotated.gtf   | awk '($5 - $4) > 199'   > lncRNA_candidates.gtf
```

- **Class codes** filter novel/antisense classes.  
- **Length filter** ensures >200 bp.

After this filter step there are 18998 lncRNA candidates.

---

## 9. Fasta Extraction & CPC2

Bedtools

```bash
"${BEDTOOLS_DIR}"/bedtools getfasta \
-fi "${GENOME_FASTA}" \
-bed "${OUTPUT_DIR}/lncRNA_candidates.gtf" \
-fo "${OUTPUT_DIR}/lncRNA_candidates.fasta" \
-name -split
```

```bash
fgrep -c ">" ${OUTPUT_DIR}/lncRNA_candidates.fasta

head ${OUTPUT_DIR}/lncRNA_candidates.fasta
```

CPC2

```bash
eval "$(/opt/anaconda/anaconda3/bin/conda shell.bash hook)"
python /home/shared/CPC2_standalone-1.0.1/bin/CPC2.py \
-i "${OUTPUT_DIR}/lncRNA_candidates.fasta" \
-o "${OUTPUT_DIR}/CPC2"
```

Filter

```bash
awk '$8 == "noncoding" {print $1}' "${OUTPUT_DIR}/CPC2.txt" > "${OUTPUT_DIR}/noncoding_transcripts_ids.txt"
```

Subsetting new fasta

```bash
"${SAMTOOLS_DIR}samtools" faidx "${OUTPUT_DIR}/lncRNA_candidates.fasta" \
-r "${OUTPUT_DIR}/noncoding_transcripts_ids.txt" \
> "${OUTPUT_DIR}/lncRNA.fasta"
```

Generate new bed and new gtf

```bash
# Define input and output file paths using the OUTPUT_DIR variable
input="${OUTPUT_DIR}/noncoding_transcripts_ids.txt"
output="${OUTPUT_DIR}/lncRNA.bed"

# Process each line of the input file
while IFS= read -r line; do
    # Remove "transcript::" from the line
    line="${line//transcript::/}"
    
    # Split the line by ':' to get the chromosome and position string
    IFS=':' read -r chromosome pos <<< "$line"
    
    # Split the position string by '-' to separate start and end positions
    IFS='-' read -r start end <<< "$pos"
    
    # Convert the start position to 0-based by subtracting 1
    start=$((start - 1))
    
    # Write the chromosome, updated start, and end positions to the output file (tab-separated)
    printf "%s\t%s\t%s\n" "$chromosome" "$start" "$end"
done < "$input" > "$output"
```

```bash
awk 'BEGIN{OFS="\t"; count=1} {printf "%s\t.\tlncRNA\t%d\t%d\t.\t+\t.\tgene_id \"lncRNA_%03d\";\n", $1, $2, $3, count++;}' "${OUTPUT_DIR}/lncRNA.bed" \
> "${OUTPUT_DIR}/lncRNA.gtf"
```

- **Extract** nucleotide sequences of candidates.  
- **Predict** coding potential; retain only noncoding.
- **Generate** new gtf, bed, and fasta of lncRNA candidates

After this filter step there are 16153 lncRNA candidates in all three output file types.

---

## 10. GTF format editing & featureCounts for expression matrices

New reliable position change step. +1 to bed positions to fit GTF expectation for downstream analysis.

```bash
awk 'BEGIN{OFS="\t"} 
     { $4 = $4 + 1; print }' \
  ~/github/deep-dive-expression/M-multi-species/output/01.6-Ptuh-lncRNA-pipeline/lncRNA.gtf \
> ~/github/deep-dive-expression/M-multi-species/output/01.6-Ptuh-lncRNA-pipeline/Ptuh_lncRNA_fixed.gtf

```

The 16104th line has a start position error, with a true start of zero that cannot be fed into featureCounts. This chunk only fixes this one line to change start from 0 to 1.

```bash
awk 'BEGIN{OFS="\t"}
     # if the GTF start (col 4) is 0, set it to 1
     { if($4==0) $4=1
       print
     }' \
  ~/github/deep-dive-expression/M-multi-species/output/01.6-Ptuh-lncRNA-pipeline/Ptuh_lncRNA_fixed.gtf \
> ~/github/deep-dive-expression/M-multi-species/output/01.6-Ptuh-lncRNA-pipeline/Ptuh_lncRNA_zeroFix.gtf

```

Issue with column 9. Original gtf generation made lncRNA ID and gene_id two separate columns, need to be merge into one for proper formatting.

```bash

awk 'BEGIN{OFS="\t"}
     {
       # Gather everything from col 9 onward into one string
       attr = $9
       for(i=10; i<=NF; i++){
         attr = attr " " $i
       }
       # Now print exactly nine columns, with our recombined attr
       print $1, $2, $3, $4, $5, $6, $7, $8, attr
     }' \
  ~/github/deep-dive-expression/M-multi-species/output/01.6-Ptuh-lncRNA-pipeline/Ptuh_lncRNA_zeroFix.gtf \
> ~/github/deep-dive-expression/M-multi-species/output/01.6-Ptuh-lncRNA-pipeline/Ptuh_lncRNA_for_fc.gtf

```

GTF is now ready for input to featureCounts for expression matrices

```bash
/home/shared/subread-2.0.5-Linux-x86_64/bin/featureCounts \
-T 24 \
-a ../output/01.6-Ptuh-lncRNA-pipeline/Ptuh_lncRNA_for_fc.gtf \
-o ../output/01.6-Ptuh-lncRNA-pipeline/Ptuh_counts.txt \
-t lncRNA \
-g gene_id \
-p \
../output/01.6-Ptuh-lncRNA-pipeline/*sorted.bam

```

Edits will still need to be incorporated to the original final gtf generation now that formatting issues have been diagnosed for successful input into featureCounts.