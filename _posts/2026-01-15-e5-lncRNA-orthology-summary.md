---
layout: post
title: "E5 Cross-species lncRNA conservation analysis & summary"
author: "Zachary Bengtsson"
date: 2026-01-15
categories: [lncRNA, comparative-genomics, expression]
tags: [lncRNA, BLAST, coral, E5, RNA-seq]
editor_options: 
  markdown: 
    wrap: 72
---

[Repo
link](https://github.com/urol-e5/deep-dive-expression/blob/main/M-multi-species/code/13-lncRNA-cross-species.Rmd)

## Overview

This notebook documents a cross-species analysis of long non-coding RNAs
(lncRNAs) across *Acropora pulchra*, *Porites evermanni*, and
*Pocillopora tuahiniensis*. The goals were to:

1.  Quantify basic properties of lncRNAs within each species (length and
    expression)
2.  Identify putatively conserved lncRNAs using sequence similarity
3.  Visualize whether conserved lncRNAs show distinct length or
    expression patterns

------------------------------------------------------------------------

## 1. Libraries and directory structure

Loading core R packages and defining all output directories up front.
Centralizing paths avoids hard-coding and makes the workflow portable
across machines and HPC environments.

``` r
library(tidyverse)
library(Biostrings)
library(purrr)

blast_dir      <- "~/github/deep-dive-expression/M-multi-species/output/13-lncRNA-cross-species/blast"
out_table_dir  <- "~/github/deep-dive-expression/M-multi-species/output/13-lncRNA-cross-species/length_expr_tables"
plot_dir       <- "~/github/deep-dive-expression/M-multi-species/output/13-lncRNA-cross-species/plots"

dir.create(out_table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir,      showWarnings = FALSE, recursive = TRUE)
```

------------------------------------------------------------------------

## 2. Data acquisition

Didn't want to go back and confirm files so used a quick and dirty cd
download DON'T JUDGE ME. lncRNA FASTA files and filtered count matrices
were downloaded directly from GitHub.

``` bash
DEST="$HOME/github/deep-dive-expression/M-multi-species/data/13-lncRNA-cross-species"
cd "$DEST"

URLS=(
  "https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/D-Apul/output/31-Apul-lncRNA/Apul-lncRNA.fasta"
  "https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/E-Peve/output/17-Peve-lncRNA/Peve-lncRNA.fasta"
  "https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/F-Ptuh/output/17-Ptuh-lncRNA/Ptuh-lncRNA.fasta"
  "https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/M-multi-species/output/01.6-lncRNA-pipeline/Apul-lncRNA-counts-filtered.txt"
  "https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/M-multi-species/output/01.6-lncRNA-pipeline/Peve-lncRNA-counts-filtered.txt"
  "https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/M-multi-species/output/01.6-lncRNA-pipeline/Ptuh-lncRNA-counts-filtered.txt"
)

for u in "${URLS[@]}"; do
  curl -L --fail -O "$u"
done
```

------------------------------------------------------------------------

## 3. Species metadata

Short species codes are used internally to simplify joins and iteration.
Publication-ready names are introduced *only at the plotting stage*.

``` r
species_meta <- tibble::tribble(
  ~species, ~fasta, ~counts,
  "Apul", "~/github/deep-dive-expression/M-multi-species/data/13-lncRNA-cross-species/Apul-lncRNA.fasta",
          "~/github/deep-dive-expression/M-multi-species/data/13-lncRNA-cross-species/Apul-lncRNA-counts-filtered.txt",
  "Peve", "~/github/deep-dive-expression/M-multi-species/data/13-lncRNA-cross-species/Peve-lncRNA.fasta",
          "~/github/deep-dive-expression/M-multi-species/data/13-lncRNA-cross-species/Peve-lncRNA-counts-filtered.txt",
  "Ptuh", "~/github/deep-dive-expression/M-multi-species/data/13-lncRNA-cross-species/Ptuh-lncRNA.fasta",
          "~/github/deep-dive-expression/M-multi-species/data/13-lncRNA-cross-species/Ptuh-lncRNA-counts-filtered.txt"
)
```

------------------------------------------------------------------------

## 4. Helper functions

### Extracting lncRNA lengths

FASTA headers include species-specific prefixes (e.g.
`Apul_lncRNA_00001`). These prefixes are stripped to ensure consistent
IDs across FASTA, counts, and BLAST.

``` r
get_lnc_lengths <- function(fasta_path) {
  dna <- Biostrings::readDNAStringSet(fasta_path)
  ids_clean <- sub("^[^_]+_", "", names(dna))
  tibble(lnc_id = ids_clean, length_nt = as.integer(width(dna)))
}
```

### Summarizing expression

Expression is summarized as the **mean across samples**, providing a
single value per lncRNA suitable for cross-species visualization.

``` r
get_lnc_expression <- function(counts_path) {
  counts <- readr::read_tsv(counts_path, show_col_types = FALSE)

  if ("Geneid" %in% names(counts)) {
    counts <- dplyr::rename(counts, lnc_id = Geneid)
  } else if (!"lnc_id" %in% names(counts)) {
    names(counts)[1] <- "lnc_id"
  }

  meta_cols <- c("lnc_id", "Chr", "Start", "End", "Strand", "Length")

  expr_mat <- counts %>%
    dplyr::select(-dplyr::any_of(meta_cols)) %>%
    dplyr::mutate(across(everything(), as.numeric))

  tibble(
    lnc_id          = counts$lnc_id,
    mean_expr       = rowMeans(expr_mat, na.rm = TRUE),
    log10_mean_expr = log10(rowMeans(expr_mat, na.rm = TRUE) + 1)
  )
}
```

------------------------------------------------------------------------

## 5. Conservation inference using BLAST

Pairwise BLASTn comparisons were performed between all species. For each
query lncRNA, **only the best hit** was retained after filtering.

**Key filtering criteria** - Percent identity ≥ 70% - Query coverage ≥
50% - E-value ≤ 1e-5

This avoids inflating conservation counts due to paralogs or weak
matches.

``` r
read_blast_best_hits <- function(path,
                                 pident_min = 70,
                                 qcov_min   = 0.5,
                                 evalue_max = 1e-5) {
  ...
}
```

------------------------------------------------------------------------

## 6. Building a unified lncRNA table

Lengths, expression summaries, and conservation status are merged into a
single tidy table (`lnc_df`), with one row per lncRNA per species. This
table underlies all downstream analyses.

------------------------------------------------------------------------

## 7. Visualization

### Length distributions

lncRNA lengths span multiple orders of magnitude, so a log10 scale is
used. Histograms are faceted by species and conservation status.

![](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/lncRNA-length-distribution.png?raw=true)

### Length vs expression

Scatterplots highlight conserved lncRNAs in bright red, while
non-conserved lncRNAs are visually de-emphasized using transparency.

![](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/lncRNA-length-v-expression.png?raw=true)

------------------------------------------------------------------------

## 8. Summary statistics

Per-species summary statistics are computed to support results text and
tables.

``` r
lncRNA_summary <- lnc_df %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_lncRNAs = n(),
    mean_length_bp = round(mean(length_nt, na.rm = TRUE), 0),
    median_length_bp = round(median(length_nt, na.rm = TRUE), 0),
    mean_expression = round(mean(mean_expr, na.rm = TRUE), 1),
    median_expression = round(median(mean_expr, na.rm = TRUE), 1),
    .groups = "drop"
  )
```

| Species | n_lncRNAs | Mean length (bp) | Median length (bp) | Mean expression | Median expression |
|--------:|----------:|-----------------:|-------------------:|----------------:|------------------:|
| Apul    | 31,491    | 2,397            | 761                | 305.8           | 34.3              |
| Peve    | 10,090    | 2,591            | 977                | 2,907.7         | 46.5              |
| Ptuh    | 16,153    | 3,124            | 700                | 663.9           | 28.3              |

------------------------------------------------------------------------

## E5 Deep Dive Expression Manuscript Text Summarizing Results

A total of 31491 putative long non-coding RNAs (lncRNAs) were identified in A. pulchra, with an average transcript length of 2397 base pairs (bp), a mean expression level of 306 counts, and a median expression of 34 counts. In P. evermanni, 10090 lncRNAs were detected, with an average length of 2591 bp, a mean expression level of 2908 counts, and a median of 47 counts. For P. tuahiniensis, 16153 lncRNAs were identified, with an average length of 3124 bp, a mean expression level of 664 TPM, and a median of 28 TPM. Shorter lncRNAs, those below 2000 bp, are more highly expressed across all three species.  Multi-species comparison identified 46 unique lncRNA transcripts shared among the three species, while 205, 190, and 65 were shared among A. pulchra and P. evermanni, A. pulchra and P. tuahiniensis, and P. evermanni and P. tuahiniensis, respectively. The majority of conserved lncRNAs are below 2000 bp  and vary widely in mean expression.

------------------------------------------------------------------------

## Common pitfalls encountered

-   Mismatched lncRNA IDs due to species prefixes
-   Treating all BLAST hits as conserved (resolved via best-hit
    filtering)
-   Overplotting in scatterplots (resolved via alpha and color choices)

------------------------------------------------------------------------

## Take-home message

This workflow integrates sequence similarity, expression summarization,
and visualization to provide a reproducible and biologically
interpretable view of lncRNA conservation across coral species. The
resulting figures and tables directly support comparative and
evolutionary analyses of non-coding RNAs.
