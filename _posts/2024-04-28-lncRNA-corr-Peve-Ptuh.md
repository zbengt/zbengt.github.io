---
title: 04.28.2025 lncRNA–mRNA Correlation Analysis (Apul, Peve, Ptuh)
date: '2025-04-28'
categories: E5
tags: pipeline
---

Full code in E5 deep-dive-expression repo [here](https://github.com/urol-e5/deep-dive-expression/blob/main/M-multi-species/code/03-expression-matrix-correlations.qmd)

## Overview

This post documents a complete R-based pipeline for analyzing the correlation between long non-coding RNAs (lncRNAs) and mRNAs across three marine invertebrate species: **Apul**, **Peve**, and **Ptuh**. The analysis includes:

-   Reading in raw count matrices for lncRNAs and mRNAs
-   Cleaning and harmonizing sample names
-   Filtering low-expression genes
-   Applying DESeq2 normalization and variance-stabilizing transformation
-   Computing Pearson correlations in a block-wise fashion using RAM-efficient matrices
-   Exporting full and significant correlation results
-   Visualizing the top 20 significant lncRNA–mRNA correlations via heatmaps

------------------------------------------------------------------------

## 1. Libraries and Setup

All three analyses begin by loading the necessary R packages:

``` r
library(data.table)
library(DESeq2)
library(WGCNA)      # for correlation p-values
if (!requireNamespace("bigmemory", quietly = TRUE)) install.packages("bigmemory")
library(bigmemory)
library(dplyr)
```

------------------------------------------------------------------------

## 2. Data Input and Cleaning

### Apul

-   Loads lncRNA and mRNA count matrices from the local filesystem
-   Removes blank Gene IDs and enforces uniqueness
-   Harmonizes sample names using regex
-   Rounds Kallisto-generated mRNA counts

### Peve and Ptuh

-   Uses `curl` to download mRNA matrices from GitHub
-   Performs the same harmonization and cleaning steps
-   Each sample name is standardized via regex from raw filenames

------------------------------------------------------------------------

## 3. Matching and Filtering

For all three species: - Retain only samples common to both lncRNA and mRNA matrices - Apply filters to exclude genes expressed in fewer than 3 samples with fewer than 10 counts

``` r
min_samples <- 3
min_counts  <- 10
keep_lnc  <- rowSums(lnc_counts  >= min_counts) >= min_samples
keep_mrna <- rowSums(mrna_counts >= min_counts) >= min_samples
```

------------------------------------------------------------------------

## 4. DESeq2 Normalization and VST

DESeq2 is used to normalize combined count matrices and extract variance-stabilized data:

``` r
dds <- DESeqDataSetFromMatrix(countData = combined,
                              colData = data.frame(cond = factor(rep("all", ncol(combined))),
                              row.names = colnames(combined)),
                              design = ~1)
dds <- estimateSizeFactors(dds)
vst_mat <- assay(vst(dds, blind = TRUE))
```

------------------------------------------------------------------------

## 5. Pearson Correlation (Block-wise)

To manage memory, correlations are calculated in blocks:

``` r
for (start in seq(1, nr_lnc, by = block_size)) {
  ...
  cor_mat [start:end, ] <- rblk
  cor_padj[start:end, ] <- p.adjust(pblk, method = "BH")
}
```

-   Correlations are calculated using `cor()`
-   P-values adjusted with Benjamini–Hochberg (BH) procedure

------------------------------------------------------------------------

## 6. Export Results

-   Full correlation and adjusted p-value matrices are exported to `.tsv`
-   Significant pairs (\|r\| ≥ 0.7 & padj ≤ 0.05) are extracted and saved separately

``` r
fwrite(sig, "significant_pairs.tsv")
```

------------------------------------------------------------------------

## 7. Visualization

Top 20 most significant lncRNA–mRNA pairs are visualized using a heatmap:

``` r
ggplot(heat_long, aes(x = mRNA, y = lncRNA, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(...) +
  labs(title = "Top lncRNA–mRNA Correlations")
```

------------------------------------------------------------------------

## 8. Session Info

Each species’ script ends with:

``` r
sessionInfo()
```

This ensures reproducibility and version tracking.

------------------------------------------------------------------------

## Summary

This workflow efficiently computes transcript correlations across three species while accounting for memory constraints. By combining DESeq2 and WGCNA in a block-wise fashion and harmonizing sample names from disparate pipelines, we ensure robust comparisons of lncRNA-mRNA co-expression patterns across marine invertebrate datasets.
