---
title: "07.01.2025 Adaptation of Pearson's correlation code for E5 lncRNA"
date: 2025-07-01
layout: post
categories: E5 lncRNA bioinformatics
---

TL;DR - Adapted the PCC code Jill and Kathleen have been using. Should work on raven for all mRNA-lncRNA correlations. Ptuh used in example.

## 1. Setup

Load all the packages we’ll need and turn on code echoing.

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

library(tidyverse)    # for data wrangling & piping
library(reshape2)     # for wide⇄long reshaping (optional)
library(igraph)       # graph object
library(tidygraph)    # tidy wrapper around igraph
library(ggraph)       # ggplot2-style network plotting
```

------------------------------------------------------------------------

## 2. Read & Filter Count Matrices

### 2.1 lncRNA

```{r read-lncrna}
# Adjust path to wherever your lncRNA counts live
lncRNA_counts_raw <- read.delim("../output/18-Ptuh-lncRNA-matrix/Ptuh-lncRNA-counts.txt", skip=1)

# Rename & drop coordinate columns, then filter out any lncRNA with zero total counts
lncRNA_counts_df <- lncRNA_counts_raw %>%
  rename(lncRNA_id = Geneid) %>%
  select(lncRNA_id, starts_with("RNA-POC")) %>% 
  column_to_rownames("lncRNA_id") %>%
  filter(rowSums(.) != 0)
```

### 2.2 mRNA

```{r read-mrna}
# Replace with your own path & column names
mrna_counts_raw <- read.delim("../output/03.1-Ptuh-sRNA-summary/Ptuh_mRNA_counts.txt")

mrna_counts_df <- mrna_counts_raw %>%
  rename(mRNA_id = GeneID) %>%     # adjust GeneID to match your sheet
  column_to_rownames("mRNA_id") %>%
  select(everything()) %>%
  filter(rowSums(.) != 0)
```

------------------------------------------------------------------------

## 3. Normalize to RPM

Define a simple RPM function and apply it to both matrices.

```{r normalize}
normalize_counts <- function(counts) {
  t(t(counts) / colSums(counts)) * 1e6
}

lncRNA_norm <- normalize_counts(lncRNA_counts_df)
mrna_norm   <- normalize_counts(mrna_counts_df)
```

------------------------------------------------------------------------

## 4. Generate All lncRNA–mRNA Pairs

We build a two‐column data frame where each row is one possible lncRNA–mRNA combination.

```{r pair-generation}
pairs_lnc_mrna <- expand.grid(
  lncRNA = rownames(lncRNA_norm),
  mRNA   = rownames(mrna_norm),
  stringsAsFactors = FALSE
)
```

------------------------------------------------------------------------

## 5. Compute Pearson’s *r* & *p*-value

1.  Define a helper that runs `cor.test()` on two numeric vectors.\
2.  Rowwise, feed each pair into that function.\
3.  Unpack the results, adjust *p*-values for FDR, and then filter significant hits.

```{r}
# 5.1 Helper function
calc_pcc <- function(x, y) {
  res <- cor.test(x, y, method = "pearson")
  c(PCC       = unname(res$estimate),
    p_value   = res$p.value)
}

# 5.2 Run one test per pair
pcc_lnc_mrna <- pairs_lnc_mrna %>%
  rowwise() %>%
  mutate(
    stats = list(calc_pcc(
      lncRNA_norm[lncRNA, ],
      mrna_norm[mRNA, ]
    ))
  ) %>%
  unnest_wider(stats) %>%
  ungroup()

# 5.3 FDR adjustment & filtering
pcc_lnc_mrna <- pcc_lnc_mrna %>%
  mutate(adj_p = p.adjust(p_value, method = "fdr"))

sig_lnc_mrna <- pcc_lnc_mrna %>%
  filter(p_value < 0.05)
```

------------------------------------------------------------------------

## 6. Save Results

Write out both the full table and the significant‐only subset to CSV.

```{r save-results}
write.csv(pcc_lnc_mrna, "../output/15-Ptuh-lncRNA-mRNA-PCC/all_pairs.csv",         row.names = FALSE)
write.csv(sig_lnc_mrna, "../output/15-Ptuh-lncRNA-mRNA-PCC/significant_pairs.csv", row.names = FALSE)
```

------------------------------------------------------------------------

## 7. Quick Inspection

Confirm how many unique genes made it through significance filtering.

```{r inspect}
cat("Unique lncRNAs with significant partners: ", n_distinct(sig_lnc_mrna$lncRNA), "
")
cat("Unique mRNAs  with significant partners: ", n_distinct(sig_lnc_mrna$mRNA),   "
")
```
