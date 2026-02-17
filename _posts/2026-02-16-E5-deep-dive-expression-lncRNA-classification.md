---
layout: post
title: "E5 Deep Dive Expression lncRNA Classification"
author: "Zach Bengtsson"
date: "2026-02-16"
categories: [E5, lncRNA]
tags: [lncRNA, E5 Deep Dive Expression]
---

## Overview

This notebook documents the **Tier 1 lncRNA classification** workflow used for three coral species (**Apul, Peve, Ptuh**). The goal is to (1) classify lncRNAs by **genomic context** (intergenic/intronic/exonic; sense vs antisense) using **BEDTools**, and (2) classify lncRNAs as **cis** vs **non-cis/unknown** based on proximity to the nearest gene and expression correlation.

**Tier 1 outputs (per species):** - `lnc_tx.bed`, `gene_exons.bed`, `gene_body.bed` - overlap files: `lnc_vs_geneExons_sense.bed`, `lnc_vs_geneExons_antisense.bed`, `lnc_vs_geneBody_intronicCandidates.bed` - nearest-gene file: `lnc_nearestGene.bed` - annotation table: `<species>_lnc_Tier1_annotation.tsv`

**Combined output:** - `AllSpecies_lnc_Tier1_annotation.tsv`

------------------------------------------------------------------------

## 1) R environment

### 1.1 Global chunk options

These options keep the notebook verbose enough for debugging while still readable.

``` r
knitr::opts_chunk$set(
  echo    = TRUE,
  message = TRUE,
  warning = TRUE
)
```

### 1.2 Libraries

We use `tidyverse` for wrangling/plots and `purrr` for safe iteration over species.

``` r
library(tidyverse)
library(purrr)
```

------------------------------------------------------------------------

## 2) Paths + parameters

This block defines: - which species to run, - where GTFs + expression matrices live, - where outputs go, - and the Tier 1 cis-calling parameters.

``` r
# Species to process
species <- c("Apul", "Peve", "Ptuh")

# Base directory for GTFs
BASE_GTF_DIR <- "~/github/deep-dive-expression/M-multi-species/data/14-lncRNA-classification/GTFs"

# Base directory for classification outputs (BED + annotation tables)
OUT_BASE <- "~/github/deep-dive-expression/M-multi-species/output/14-lncRNA-classification/coral_lnc_classification"

# Base directory for expression matrices
# Expecting files:
#   <species>_lnc_expr.tsv   (tab-delimited, but may originate as .txt)
#   <species>_gene_expr.tsv  (comma-delimited CSV saved with .tsv extension here)
EXPR_BASE <- "~/github/deep-dive-expression/M-multi-species/data/14-lncRNA-classification/expression_matrices"

# Tier 1 parameters
cis_window_bp <- 10000   # ±10 kb around nearest gene
cor_threshold <- 0.6     # |r| >= 0.6 to call cis
```

------------------------------------------------------------------------

## 3) Download inputs (GTFs + expression matrices)

This section pins input files into the directory structure expected by the R code above.

### 3.1 Download GTFs

This shell chunk fetches both: - validated gene models (`*_genes.gtf`) - lncRNA annotation GTFs (`*-lncRNA.gtf`)

``` bash
# =========================
# GTF downloads
# =========================

# These must match the R paths above
BASE_GTF_DIR=$HOME/github/deep-dive-expression/M-multi-species/data/14-lncRNA-classification/GTFs

mkdir -p "${BASE_GTF_DIR}"

# ---- Apul ----
curl -L   https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/data/Apulchra-genome.gtf   -o "${BASE_GTF_DIR}/Apul_genes.gtf"

curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/D-Apul/output/31-Apul-lncRNA/Apul-lncRNA.gtf   -o "${BASE_GTF_DIR}/Apul-lncRNA.gtf"


# ---- Peve ----
curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/b0290e08af4eaeed30d74a758965debef6111801/E-Peve/data/Porites_evermanni_validated.gtf   -o "${BASE_GTF_DIR}/Peve_genes.gtf"

curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/E-Peve/output/17-Peve-lncRNA/Peve-lncRNA.gtf   -o "${BASE_GTF_DIR}/Peve-lncRNA.gtf"


# ---- Ptuh ----
curl -L   https://github.com/urol-e5/deep-dive-expression/raw/f62c6d01e04ef0007f2f53af84181481d64d29c1/F-Ptuh/data/Pocillopora_meandrina_HIv1.genes-validated.gtf   -o "${BASE_GTF_DIR}/Ptuh_genes.gtf"

curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/F-Ptuh/output/17-Ptuh-lncRNA/Ptuh-lncRNA.gtf   -o "${BASE_GTF_DIR}/Ptuh-lncRNA.gtf"

echo "GTFs downloaded to ${BASE_GTF_DIR}"
```

### 3.2 Download expression matrices

**Important detail:** lnc matrices are tab-delimited text; gene matrices are CSV.\
We keep the filenames ending in `.tsv` for consistency downstream, but read them using the correct import functions later (`read_tsv()` vs `read_csv()`).

``` bash
# =========================
# Expression matrices
# =========================

EXPR_BASE=$HOME/github/deep-dive-expression/M-multi-species/data/14-lncRNA-classification/expression_matrices
mkdir -p "${EXPR_BASE}"

# ---- Apul ----
# lncRNA counts: tab-delimited .txt
curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/main/M-multi-species/output/01.6-lncRNA-pipeline/Apul-lncRNA-counts-filtered.txt   -o "${EXPR_BASE}/Apul_lnc_expr.tsv"

# gene counts: CSV
curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/main/D-Apul/output/07-Apul-Hisat/Apul-gene_count_matrix.csv   -o "${EXPR_BASE}/Apul_gene_expr.tsv"


# ---- Peve ----
curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/main/M-multi-species/output/01.6-lncRNA-pipeline/Peve-lncRNA-counts-filtered.txt   -o "${EXPR_BASE}/Peve_lnc_expr.tsv"

curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/main/E-Peve/output/06.2-Peve-Hisat/Peve-gene_count_matrix.csv   -o "${EXPR_BASE}/Peve_gene_expr.tsv"


# ---- Ptuh ----
curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/main/M-multi-species/output/01.6-lncRNA-pipeline/Ptuh-lncRNA-counts-filtered.txt   -o "${EXPR_BASE}/Ptuh_lnc_expr.tsv"

curl -L   https://raw.githubusercontent.com/urol-e5/deep-dive-expression/main/F-Ptuh/output/06.2-Ptuh-Hisat/Ptuh-gene_count_matrix.csv   -o "${EXPR_BASE}/Ptuh_gene_expr.tsv"

echo "Expression matrices downloaded to ${EXPR_BASE}"
```

------------------------------------------------------------------------

## 4) Quick input sanity check (GTF feature counts)

This optional check confirms the GTFs contain the features we rely on: - gene GTFs: `exon`, `transcript`, `gene` - lncRNA GTFs: `lncRNA`

``` bash
set -e

BASE_GTF_DIR=$HOME/github/deep-dive-expression/M-multi-species/data/14-lncRNA-classification/GTFs

echo "BASE_GTF_DIR = ${BASE_GTF_DIR}"
echo
echo "Files in BASE_GTF_DIR:"
ls -l "${BASE_GTF_DIR}"

echo
echo "=== Apul_genes.gtf feature counts ==="
awk 'BEGIN{ex=0; tx=0; gene=0}
     !/^#/ {
       if ($3=="exon") ex++
       else if ($3=="transcript") tx++
       else if ($3=="gene") gene++
     }
     END{
       print "exon:", ex;
       print "transcript:", tx;
       print "gene:", gene;
     }' "${BASE_GTF_DIR}/Apul_genes.gtf"

echo
echo "=== Apul-lncRNA.gtf feature counts ==="
awk 'BEGIN{lnc=0}
     !/^#/ {
       if ($3=="lncRNA") lnc++
     }
     END{
       print "lncRNA:", lnc;
     }' "${BASE_GTF_DIR}/Apul-lncRNA.gtf"

echo
echo "=== Peve_genes.gtf feature counts ==="
awk 'BEGIN{ex=0; tx=0; gene=0}
     !/^#/ {
       if ($3=="exon") ex++
       else if ($3=="transcript") tx++
       else if ($3=="gene") gene++
     }
     END{
       print "exon:", ex;
       print "transcript:", tx;
       print "gene:", gene;
     }' "${BASE_GTF_DIR}/Peve_genes.gtf"

echo
echo "=== Peve-lncRNA.gtf feature counts ==="
awk 'BEGIN{lnc=0}
     !/^#/ {
       if ($3=="lncRNA") lnc++
     }
     END{
       print "lncRNA:", lnc;
     }' "${BASE_GTF_DIR}/Peve-lncRNA.gtf"

echo
echo "=== Ptuh_genes.gtf feature counts ==="
awk 'BEGIN{ex=0; tx=0; gene=0}
     !/^#/ {
       if ($3=="exon") ex++
       else if ($3=="transcript") tx++
       else if ($3=="gene") gene++
     }
     END{
       print "exon:", ex;
       print "transcript:", tx;
       print "gene:", gene;
     }' "${BASE_GTF_DIR}/Ptuh_genes.gtf"

echo
echo "=== Ptuh-lncRNA.gtf feature counts ==="
awk 'BEGIN{lnc=0}
     !/^#/ {
       if ($3=="lncRNA") lnc++
     }
     END{
       print "lncRNA:", lnc;
     }' "${BASE_GTF_DIR}/Ptuh-lncRNA.gtf"
```

------------------------------------------------------------------------

## 5) Build BED files + overlaps (BEDTools)

This is the “genomic context” step. For each species, we create:

1)  `gene_exons.bed`\
    All gene-model exons (using `gene_id`), strand-aware.

2)  `gene_body.bed`\
    “Gene bodies” assembled from `gene` and/or `transcript` features (depending on GTF content).

3)  `lnc_tx.bed`\
    One entry per lncRNA transcript interval (using `gene_id` in the lncRNA GTF), strand-aware.

Then we compute overlaps:

-   **sense exonic**: `bedtools intersect -s` between `lnc_tx` and `gene_exons`
-   **antisense exonic**: `bedtools intersect -S` between `lnc_tx` and `gene_exons`
-   **intronic candidates**: `bedtools intersect -s -f 1.0` where the lncRNA is fully contained in a gene body
-   **nearest gene**: `bedtools closest -d` to get nearest gene body and distance (bp)

``` bash
set -euo pipefail

species=("Apul" "Peve" "Ptuh")

BASE_GTF_DIR=$HOME/github/deep-dive-expression/M-multi-species/data/14-lncRNA-classification/GTFs
OUT_BASE=$HOME/github/deep-dive-expression/M-multi-species/output/14-lncRNA-classification/coral_lnc_classification

# Explicit bedtools path
BEDTOOLS="/home/shared/bedtools-v2.30.0/bin/bedtools"

if [[ ! -x "${BEDTOOLS}" ]]; then
  echo "ERROR: bedtools not executable at ${BEDTOOLS}" >&2
  exit 1
fi

echo "Using bedtools at: ${BEDTOOLS}"
"${BEDTOOLS}" --version
echo

make_beds_for_species () {
  local sp="$1"
  echo ">>> Processing species (bedtools step): ${sp}"

  local GENE_GTF="${BASE_GTF_DIR}/${sp}_genes.gtf"
  local LNC_GTF="${BASE_GTF_DIR}/${sp}-lncRNA.gtf"

  if [[ ! -f "${GENE_GTF}" ]]; then
    echo "ERROR: Missing gene GTF for ${sp}: ${GENE_GTF}" >&2
    exit 1
  fi
  if [[ ! -f "${LNC_GTF}" ]]; then
    echo "ERROR: Missing lncRNA GTF for ${sp}: ${LNC_GTF}" >&2
    exit 1
  fi

  local OUTDIR="${OUT_BASE}/${sp}"
  mkdir -p "${OUTDIR}"
  cd "${OUTDIR}"

  echo "  - Creating BED files for ${sp}"

  # 1) gene_exons.bed
  awk 'BEGIN{FS=OFS="   "}
       !/^#/ && $3=="exon" {
         attr = $9
         id   = attr
         sub(/.*gene_id "/, "", id)
         sub(/".*/, "", id)
         if (id != "") {
           print $1, $4-1, $5, id, ".", $7
         }
       }' "${GENE_GTF}"     | sort -k1,1 -k2,2n -k3,3n -k6,6     > gene_exons.bed

  echo "    gene_exons.bed lines: $(wc -l < gene_exons.bed)"

  # 2) gene_body.bed
  awk 'BEGIN{FS=OFS="   "}
       !/^#/ && ($3=="gene" || $3=="transcript") {
         attr = $9
         id   = ""
         if (index(attr, "gene_id ") > 0) {
           id = attr
           sub(/.*gene_id "/, "", id)
         } else if (index(attr, "transcript_id ") > 0) {
           id = attr
           sub(/.*transcript_id "/, "", id)
         }
         sub(/".*/, "", id)
         if (id != "") {
           print $1, $4-1, $5, id, ".", $7
         }
       }' "${GENE_GTF}"     | sort -k1,1 -k2,2n -k3,3n -k6,6     > gene_body.bed

  echo "    gene_body.bed lines: $(wc -l < gene_body.bed)"

  # 3) lnc_tx.bed
  awk 'BEGIN{FS=OFS="   "}
       !/^#/ && $3=="lncRNA" {
         attr = $9
         id   = attr
         sub(/.*gene_id "/, "", id)
         sub(/".*/, "", id)
         if (id != "") {
           print $1, $4-1, $5, id, ".", $7
         }
       }' "${LNC_GTF}"     | sort -k1,1 -k2,2n -k3,3n -k6,6     > lnc_tx.bed

  echo "    lnc_tx.bed lines: $(wc -l < lnc_tx.bed)"

  echo "  - Running bedtools for ${sp}"

  "${BEDTOOLS}" intersect     -s -wa -wb     -a lnc_tx.bed -b gene_exons.bed     > lnc_vs_geneExons_sense.bed

  "${BEDTOOLS}" intersect     -S -wa -wb     -a lnc_tx.bed -b gene_exons.bed     > lnc_vs_geneExons_antisense.bed

  "${BEDTOOLS}" intersect     -s -wa -wb     -f 1.0     -a lnc_tx.bed -b gene_body.bed     > lnc_vs_geneBody_intronicCandidates.bed

  "${BEDTOOLS}" closest     -d     -a lnc_tx.bed -b gene_body.bed     > lnc_nearestGene.bed

  echo "    sense exonic:      $(wc -l < lnc_vs_geneExons_sense.bed)"
  echo "    antisense exonic:  $(wc -l < lnc_vs_geneExons_antisense.bed)"
  echo "    intronic:          $(wc -l < lnc_vs_geneBody_intronicCandidates.bed)"
  echo "    nearestGene:       $(wc -l < lnc_nearestGene.bed)"

  echo ">>> Finished bedtools step for ${sp}"
  cd - >/dev/null
}

for sp in "${species[@]}"; do
  make_beds_for_species "${sp}"
done

echo "All species processed with bedtools. Outputs in ${OUT_BASE}/{Apul,Peve,Ptuh}"
```

### 5.1 Verify bedtools outputs exist and have lines

This R chunk confirms all expected files were generated for each species.

``` r
for (sp in species) {
  sp_dir <- file.path(OUT_BASE, sp)
  cat("\n====", sp, "====\n")
  files <- c("lnc_tx.bed",
             "gene_exons.bed",
             "gene_body.bed",
             "lnc_vs_geneExons_sense.bed",
             "lnc_vs_geneExons_antisense.bed",
             "lnc_vs_geneBody_intronicCandidates.bed",
             "lnc_nearestGene.bed")
  for (f in files) {
    path <- file.path(sp_dir, f)
    if (!file.exists(path)) {
      cat("  ", f, " -> MISSING\n")
    } else {
      n_lines <- length(readr::read_lines(path))
      cat("  ", f, " ->", n_lines, "lines\n")
    }
  }
}
```

------------------------------------------------------------------------

## 6) Tier 1 classification (genomic_class + cis_trans_class)

Tier 1 includes two labels:

### 6.1 `genomic_class` (from overlaps)

-   `sense_exonic`: overlaps an exon on the same strand (`intersect -s`)
-   `antisense_exonic`: overlaps an exon on the opposite strand (`intersect -S`)
-   `intronic`: fully contained in a gene body on the same strand (`intersect -s -f 1.0`), *and not already called exonic*
-   `intergenic`: anything else

### 6.2 `cis_trans_class` (from nearest gene + expression)

We call an lncRNA **cis** if:

-   its nearest gene is within `cis_window_bp` (default 10 kb), **and**
-   the lncRNA’s expression is correlated with that nearest gene with `|r| >= cor_threshold` (default 0.6)

Otherwise it becomes `non_cis_or_unknown`.

------------------------------------------------------------------------

## 7) Helper: safe correlation

This function avoids failures when either vector is all NA.

``` r
safe_cor <- function(x, y) {
  if (all(is.na(x)) | all(is.na(y))) return(NA_real_)
  suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
}
```

------------------------------------------------------------------------

## 8) Core function: `process_species()`

This function runs Tier 1 logic for a single species:

-   loads overlap/nearest BED outputs
-   assigns `genomic_class`
-   loads expression matrices and harmonizes sample names
-   computes nearest gene correlation (when possible)
-   assigns `cis_trans_class`
-   writes a per-species TSV

``` r
process_species <- function(sp,
                            bed_base = OUT_BASE,
                            expr_base = EXPR_BASE) {
  message(">>> Tier 1 classification for species: ", sp)

  sp_dir <- file.path(bed_base, sp)

  # BED / overlap files (from bash chunk)
  sense_file     <- file.path(sp_dir, "lnc_vs_geneExons_sense.bed")
  antisense_file <- file.path(sp_dir, "lnc_vs_geneExons_antisense.bed")
  intronic_file  <- file.path(sp_dir, "lnc_vs_geneBody_intronicCandidates.bed")
  nearest_file   <- file.path(sp_dir, "lnc_nearestGene.bed")
  lnc_tx_file    <- file.path(sp_dir, "lnc_tx.bed")

  # Expression matrices
  lnc_expr_file  <- file.path(expr_base, paste0(sp, "_lnc_expr.tsv"))
  gene_expr_file <- file.path(expr_base, paste0(sp, "_gene_expr.tsv"))

  # Sanity check
  files_to_check <- c(sense_file, antisense_file, intronic_file,
                      nearest_file, lnc_tx_file,
                      lnc_expr_file, gene_expr_file)
  missing <- files_to_check[!file.exists(files_to_check)]
  if (length(missing) > 0) {
    stop("Missing files for species ", sp, ":\n",
         paste(missing, collapse = "\n"))
  }

  #---------------------------
  # Load BED / overlap files
  #---------------------------
  sense <- readr::read_tsv(
    sense_file,
    col_names = c("lnc_chr","lnc_start","lnc_end","lnc_id","lnc_score","lnc_strand",
                  "gene_chr","gene_start","gene_end","gene_id","gene_score","gene_strand"),
    show_col_types = FALSE
  )

  antisense <- readr::read_tsv(
    antisense_file,
    col_names = c("lnc_chr","lnc_start","lnc_end","lnc_id","lnc_score","lnc_strand",
                  "gene_chr","gene_start","gene_end","gene_id","gene_score","gene_strand"),
    show_col_types = FALSE
  )

  intronic <- readr::read_tsv(
    intronic_file,
    col_names = c("lnc_chr","lnc_start","lnc_end","lnc_id","lnc_score","lnc_strand",
                  "gene_chr","gene_start","gene_end","gene_id","gene_score","gene_strand"),
    show_col_types = FALSE
  )

  nearest <- readr::read_tsv(
    nearest_file,
    col_names = c("lnc_chr","lnc_start","lnc_end","lnc_id","lnc_score","lnc_strand",
                  "gene_chr","gene_start","gene_end","gene_id","gene_score","gene_strand",
                  "distance_bp"),
    show_col_types = FALSE
  )

  # All lncRNAs
  lnc_all <- readr::read_tsv(
    lnc_tx_file,
    col_names = c("lnc_chr","lnc_start","lnc_end","lnc_id","lnc_score","lnc_strand"),
    show_col_types = FALSE
  ) %>%
    dplyr::distinct(lnc_id, .keep_all = TRUE)

  #---------------------------
  # Genomic class
  #---------------------------
  sense_ids     <- unique(sense$lnc_id)
  antisense_ids <- unique(antisense$lnc_id)
  intronic_ids  <- unique(intronic$lnc_id)

  genomic_class_df <- lnc_all %>%
    dplyr::mutate(
      genomic_class = dplyr::case_when(
        lnc_id %in% sense_ids ~ "sense_exonic",
        lnc_id %in% antisense_ids ~ "antisense_exonic",
        lnc_id %in% intronic_ids &
          !lnc_id %in% c(sense_ids, antisense_ids) ~ "intronic",
        TRUE ~ "intergenic"
      )
    )

  #---------------------------
  # Nearest gene info
  #---------------------------
  nearest_slim <- nearest %>%
    dplyr::select(lnc_id, nearest_gene_id = gene_id, distance_bp)

  lnc_annot <- genomic_class_df %>%
    dplyr::left_join(nearest_slim, by = "lnc_id")

  #---------------------------
  # Expression matrices
  #---------------------------
  lnc_expr <- readr::read_tsv(lnc_expr_file, show_col_types = FALSE)
  gene_expr <- readr::read_csv(gene_expr_file, show_col_types = FALSE)

  lnc_ids  <- lnc_expr[[1]]
  gene_ids <- gene_expr[[1]]

  # Drop ID col + lnc metadata cols if present
  lnc_data <- lnc_expr[, -1]
  meta_cols <- c("Chr", "Start", "End", "Strand", "Length")
  lnc_data <- dplyr::select(lnc_data, -dplyr::any_of(meta_cols))

  # Clean lnc sample names (pull RNA.ACR.140-like IDs, normalize separators)
  lnc_sample_names <- colnames(lnc_data)
  lnc_sample_names_clean <- vapply(
    lnc_sample_names,
    function(n) {
      extracted <- sub(".*(RNA[._-][A-Z]{3}[._-][0-9]+).*", "\\1", n)
      if (identical(extracted, n)) return(n)
      gsub("[._]", "-", extracted)
    },
    character(1)
  )
  colnames(lnc_data) <- lnc_sample_names_clean

  # Clean gene sample names (featureCounts prefixes)
  gene_data <- gene_expr[, -1]
  gene_sample_names_clean <- colnames(gene_data)
  gene_sample_names_clean <- sub("^transcript_counts\\.", "", gene_sample_names_clean)
  gene_sample_names_clean <- sub("^transcript.counts\\.", "", gene_sample_names_clean)
  gene_sample_names_clean <- sub("^counts\\.", "", gene_sample_names_clean)
  colnames(gene_data) <- gene_sample_names_clean

  # Align sample sets
  common_samples <- intersect(colnames(lnc_data), colnames(gene_data))
  if (length(common_samples) < 2) {
    stop(
      "Too few overlapping samples between lnc and gene matrices for ", sp, ".\n",
      "lnc samples (cleaned): ", paste(colnames(lnc_data), collapse = ", "), "\n",
      "gene samples (cleaned): ", paste(colnames(gene_data), collapse = ", ")
    )
  }

  lnc_mat <- as.matrix(lnc_data[, common_samples])
  rownames(lnc_mat) <- lnc_ids

  gene_mat <- as.matrix(gene_data[, common_samples])
  rownames(gene_mat) <- gene_ids

  stopifnot(identical(colnames(lnc_mat), colnames(gene_mat)))

  #---------------------------
  # cis vs non-cis
  #---------------------------
  lnc_cis_df <- lnc_annot %>%
    dplyr::mutate(
      in_cis_window = !is.na(distance_bp) & distance_bp <= cis_window_bp,
      nearest_gene_cor = purrr::pmap_dbl(
        list(lnc_id, nearest_gene_id, in_cis_window),
        function(lnc_id_i, gene_id_i, in_window) {
          if (is.na(gene_id_i) ||
              !lnc_id_i %in% rownames(lnc_mat) ||
              !gene_id_i %in% rownames(gene_mat)) {
            return(NA_real_)
          }
          safe_cor(lnc_mat[lnc_id_i, ], gene_mat[gene_id_i, ])
        }
      ),
      cis_flag = in_cis_window &
        !is.na(nearest_gene_cor) &
        abs(nearest_gene_cor) >= cor_threshold,
      cis_trans_class = dplyr::case_when(
        cis_flag ~ "cis",
        TRUE ~ "non_cis_or_unknown"
      )
    )

  #---------------------------
  # Final table for this species
  #---------------------------
  lnc_final <- lnc_cis_df %>%
    dplyr::mutate(species = sp) %>%
    dplyr::select(
      species,
      lnc_id,
      lnc_chr, lnc_start, lnc_end, lnc_strand,
      genomic_class,
      cis_trans_class,
      nearest_gene_id,
      distance_bp,
      nearest_gene_cor
    )

  out_file <- file.path(sp_dir, paste0(sp, "_lnc_Tier1_annotation.tsv"))
  readr::write_tsv(lnc_final, out_file)
  message("  - Written: ", out_file)

  lnc_final
}
```

------------------------------------------------------------------------

## 9) Dry run: test one species (Apul)

Before running everything, it’s useful to validate the pipeline on a single species.

``` r
tmp_Apul <- process_species("Apul")
head(tmp_Apul)
```

------------------------------------------------------------------------

## 10) Run all species + write combined table

This produces (1) three per-species TSVs and (2) the combined all-species TSV.

``` r
all_results <- purrr::map_df(species, process_species)

combined_out <- file.path(OUT_BASE, "AllSpecies_lnc_Tier1_annotation.tsv")
readr::write_tsv(all_results, combined_out)
message(">>> Combined Tier1 table written: ", combined_out)
```

### 10.1 Reload + confirm dimensions

This ensures the combined file can be read back cleanly (useful if you knit in a fresh session).

``` r
tier1_path <- file.path(
  OUT_BASE,
  "AllSpecies_lnc_Tier1_annotation.tsv"
)

all_results <- readr::read_tsv(tier1_path, show_col_types = FALSE)

print(dim(all_results))
print(head(all_results))
```

------------------------------------------------------------------------

## 11) Summary figures

### 11.1 Setup (read Tier 1 table + factor ordering)

``` r
library(tidyverse)

tier1_path <- "~/github/deep-dive-expression/M-multi-species/output/14-lncRNA-classification/coral_lnc_classification/AllSpecies_lnc_Tier1_annotation.tsv"

lnc_tier1 <- readr::read_tsv(tier1_path, show_col_types = FALSE) %>%
  mutate(
    species = factor(species, levels = c("Apul", "Peve", "Ptuh")),
    genomic_class = factor(
      genomic_class,
      levels = c("intergenic", "intronic", "antisense_exonic", "sense_exonic")
    ),
    cis_trans_class = factor(
      cis_trans_class,
      levels = c("cis", "non_cis_or_unknown")
    )
  )
```

------------------------------------------------------------------------

## Figure 1A — Counts by genomic class and species

**Question:** How many lncRNAs fall into each genomic class per species?

``` r
fig1a_counts <- lnc_tier1 %>%
  count(species, genomic_class) %>%
  ggplot(aes(x = species, y = n, fill = genomic_class)) +
  geom_col(color = "black", linewidth = 0.2) +
  labs(
    x = "Species",
    y = "Number of lncRNAs",
    fill = "Genomic class"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

fig1a_counts
# ggsave("fig1a_genomic_class_counts.png", fig1a_counts, width = 5, height = 4, dpi = 300)
```

![Figure 1A — Genomic class counts by species](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/fig1a_genomic_class_counts.png?raw=true)

*Figure 1A. Counts of Tier 1 lncRNAs per species, stratified by genomic class (intergenic, intronic, antisense_exonic, sense_exonic).*

------------------------------------------------------------------------

## Figure 1B — Proportions by genomic class and species

**Question:** What fraction of lncRNAs in each species are intergenic/intronic/exonic?

``` r
fig1b_props <- lnc_tier1 %>%
  count(species, genomic_class) %>%
  group_by(species) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = species, y = prop, fill = genomic_class)) +
  geom_col(color = "black", linewidth = 0.2) +
  labs(
    x = "Species",
    y = "Proportion of lncRNAs",
    fill = "Genomic class"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

fig1b_props
# ggsave("fig1b_genomic_class_proportions.png", fig1b_props, width = 5, height = 4, dpi = 300)
```

![Figure 1B — Genomic class proportions by species](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/fig1b_genomic_class_proportions.png?raw=true)

*Figure 1B. Proportions of Tier 1 lncRNAs per species, stratified by genomic class.*

------------------------------------------------------------------------

## Figure 2 — Cis vs non-cis by genomic class (faceted)

**Question:** Within each genomic class, how many lncRNAs are cis vs non-cis/unknown?

``` r
fig2_cis_by_class <- lnc_tier1 %>%
  count(species, genomic_class, cis_trans_class) %>%
  group_by(species, genomic_class) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = species, y = prop, fill = cis_trans_class)) +
  geom_col(color = "black", linewidth = 0.2, position = "fill") +
  facet_wrap(~ genomic_class, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Species",
    y = "Proportion within genomic class",
    fill = "Regulatory class"
  ) +
  scale_fill_brewer(palette = "Pastel1",
                    labels = c("cis", "non-cis or unknown")) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

fig2_cis_by_class
# ggsave("fig2_cis_by_genomic_class.png", fig2_cis_by_class, width = 9, height = 4, dpi = 300)
```

![Figure 2 — Cis vs non-cis within genomic class](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/fig2_cis_by_genomic_class.png?raw=true)

*Figure 2. Within each genomic class, stacked proportions of Tier 1 lncRNAs classified as cis vs non-cis/unknown across species.*

------------------------------------------------------------------------

## Figure 4 — Correlation strength for cis lncRNAs

**Question:** When we call something cis, how strong is its correlation with the nearest gene?

``` r
lnc_cis_only <- lnc_tier1 %>%
  filter(cis_trans_class == "cis") %>%
  filter(!is.na(nearest_gene_cor))
```

### Figure 4A — Cis correlation density (by species)

``` r
fig4_cis_cor <- lnc_cis_only %>%
  ggplot(aes(x = nearest_gene_cor, fill = species)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = c(-0.6, 0.6),
             linetype = "dashed", linewidth = 0.4) +
  labs(
    x = "Pearson correlation with nearest gene",
    y = "Density",
    fill = "Species"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic(base_size = 12)

fig4_cis_cor
# ggsave("fig4_cis_correlation_density.png", fig4_cis_cor, width = 6, height = 4, dpi = 300)
```

![Figure 4A — Cis correlation density](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/fig4_cis_correlation_by_class.png?raw=true)

*Figure 4A. Density of nearest-gene Pearson correlations among lncRNAs classified as cis (threshold lines at r = ±0.6), colored by species.*

### Figure 4B — Cis correlation by genomic class

``` r
fig4_cis_cor_facet <- lnc_cis_only %>%
  ggplot(aes(x = nearest_gene_cor, fill = species)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = c(-0.6, 0.6),
             linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ genomic_class, nrow = 2) +
  labs(
    x = "Pearson correlation with nearest gene",
    y = "Density",
    fill = "Species"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )

fig4_cis_cor_facet
# ggsave("fig4_cis_correlation_by_class.png", fig4_cis_cor_facet, width = 8, height = 6, dpi = 300)
```

![Figure 4B — Cis correlation by genomic class](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/fig4_cis_correlation_density.png?raw=true)

*Figure 4B. Density of nearest-gene Pearson correlations among cis lncRNAs, faceted by genomic class and colored by species.*

------------------------------------------------------------------------

## 12) Summary tables

### Table 1 — Genomic class counts + proportions

``` r
table_genomic_class <- lnc_tier1 %>%
  count(species, genomic_class, name = "n_lnc") %>%
  group_by(species) %>%
  mutate(
    total_lnc = sum(n_lnc),
    prop = n_lnc / total_lnc
  ) %>%
  ungroup()

table_genomic_class
```

| species | genomic_class    | n_lnc | total_lnc | prop       |
|---------|------------------|-------|-----------|------------|
| Apul    | intergenic       | 26249 | 33010     | 0.79518328 |
| Apul    | intronic         | 2404  | 33010     | 0.07282642 |
| Apul    | antisense_exonic | 2130  | 33010     | 0.06452590 |
| Apul    | sense_exonic     | 2227  | 33010     | 0.06746440 |
| Peve    | intergenic       | 15259 | 20095     | 0.75934312 |
| Peve    | intronic         | 1420  | 20095     | 0.07066434 |
| Peve    | antisense_exonic | 1936  | 20095     | 0.09634237 |
| Peve    | sense_exonic     | 1480  | 20095     | 0.07365016 |
| Ptuh    | intergenic       | 13552 | 17260     | 0.78516802 |
| Ptuh    | intronic         | 848   | 17260     | 0.04913094 |
| Ptuh    | antisense_exonic | 1355  | 17260     | 0.07850521 |
| Ptuh    | sense_exonic     | 1505  | 17260     | 0.08719583 |

### Table 2 — Cis vs non-cis within each genomic class

``` r
table_cis <- lnc_tier1 %>%
  count(species, genomic_class, cis_trans_class, name = "n_lnc") %>%
  group_by(species, genomic_class) %>%
  mutate(
    total_in_class = sum(n_lnc),
    prop_in_class = n_lnc / total_in_class
  ) %>%
  ungroup()

table_cis
```

| species | genomic_class     | cis_trans_class       | n_lnc | total_in_class | prop_in_class |
|----------|------------------|------------------------|-------|----------------|---------------|
| Apul     | intergenic       | cis                    | 6233  | 26249          | 0.23745666501581011 |
| Apul     | intergenic       | non_cis_or_unknown     | 20016 | 26249          | 0.7625433349841899  |
| Apul     | intronic         | cis                    | 803   | 2404           | 0.334026622296173   |
| Apul     | intronic         | non_cis_or_unknown     | 1601  | 2404           | 0.6659733777038269  |
| Apul     | antisense_exonic | cis                    | 1434  | 2130           | 0.6732394366197183  |
| Apul     | antisense_exonic | non_cis_or_unknown     | 696   | 2130           | 0.3267605633802817  |
| Apul     | sense_exonic     | cis                    | 1392  | 2227           | 0.6250561293219578  |
| Apul     | sense_exonic     | non_cis_or_unknown     | 835   | 2227           | 0.3749438706780422  |
| Peve     | intergenic       | cis                    | 5756  | 15259          | 0.3772200013107019  |
| Peve     | intergenic       | non_cis_or_unknown     | 9503  | 15259          | 0.6227799986892981  |
| Peve     | intronic         | cis                    | 682   | 1420           | 0.4802816901408451  |
| Peve     | intronic         | non_cis_or_unknown     | 738   | 1420           | 0.5197183098591549  |
| Peve     | antisense_exonic | cis                    | 1648  | 1936           | 0.8512396694214877  |
| Peve     | antisense_exonic | non_cis_or_unknown     | 288   | 1936           | 0.1487603305785124  |
| Peve     | sense_exonic     | cis                    | 1166  | 1480           | 0.7878378378378378  |
| Peve     | sense_exonic     | non_cis_or_unknown     | 314   | 1480           | 0.21216216216216216 |
| Ptuh     | intergenic       | cis                    | 3141  | 13552          | 0.23177390791027155 |
| Ptuh     | intergenic       | non_cis_or_unknown     | 10411 | 13552          | 0.7682260920897285  |
| Ptuh     | intronic         | cis                    | 277   | 848            | 0.3266509433962264  |
| Ptuh     | intronic         | non_cis_or_unknown     | 571   | 848            | 0.6733490566037735  |
| Ptuh     | antisense_exonic | cis                    | 887   | 1355           | 0.6546125461254613  |
| Ptuh     | antisense_exonic | non_cis_or_unknown     | 468   | 1355           | 0.34538745387453873 |
| Ptuh     | sense_exonic     | cis                    | 761   | 1505           | 0.5056478405315614  |
| Ptuh     | sense_exonic     | non_cis_or_unknown     | 744   | 1505           | 0.4943521594684385  |


------------------------------------------------------------------------

## Notes / gotchas

-   **Sample name harmonization** is doing a lot of work here (especially for lncRNA matrices with embedded metadata in column names). If any species fails with “Too few overlapping samples…”, print `colnames(lnc_data)` and `colnames(gene_data)` immediately after cleaning and compare.
-   **Intronic calling** is conservative: it requires the lncRNA interval to be *fully contained* in a gene body (`-f 1.0`) and on the same strand (`-s`).
-   **cis calling** is also conservative: proximity **and** correlation are required.

------------------------------------------------------------------------
