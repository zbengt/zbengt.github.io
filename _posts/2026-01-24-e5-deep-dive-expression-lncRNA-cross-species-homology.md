---
layout: post
title: "E5 Deep-Dive Expression Cross-species lncRNA homology in corals (strict vs relaxed)"
date: 2026-01-24
tags: [lncRNA, deep-dive]
categories: [coral, comparative-genomics]
---

Code repo link: [13.1-lncRNA-cross-species.Rmd](https://github.com/urol-e5/deep-dive-expression/blob/main/M-multi-species/code/13.1-lncRNA-cross-species.Rmd)

Compared lncRNAs across *Acropora pulchra* (Apul), *Porites evermanni* (Peve), and *Porites tuahiniensis* (Ptuh) — to identify:

-   **1:1 homologous lncRNAs between species** (pairwise BRHs), and\
-   **1:1:1 homologous lncRNAs across all three species** (triplets).

First used **strict, protein-like BLAST filters**, then repeated the analysis with **relaxed, lncRNA-appropriate thresholds**, all using the same BLAST results.

## Summary

-   **Strict criteria** (≥70% identity, ≥50% coverage, ≤ 1e-5) identify a very small set of strongly conserved lncRNAs.

![](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/13.1-STRICT-histo.png?raw=true)

![](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/13.1-STRICT-scatter.png?raw=true)

-   **Relaxed criteria** (≥50% identity, ≥30% coverage, ≤ 1e-5) expand the set of shared lncRNAs while still enforcing **1:1 and 1:1:1 relationships**.

![](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/13.1-relaxed-histo.png?raw=true)

![](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/E5-deep-dive-expression/13.1-relaxed.scatter.png?raw=true)

| Species | Category           | Strict n | Relaxed n | Species total |
|---------|--------------------|----------|-----------|---------------|
| Apul    | Species-only       | 31,465   | 31,426    | 31,491        |
| Apul    | Apul–Peve          | 9        | 28        | 31,491        |
| Apul    | Apul–Ptuh          | 16       | 36        | 31,491        |
| Apul    | Apul–Peve–Ptuh     | 1        | 1         | 31,491        |
| Peve    | Species-only       | 10,074   | 10,048    | 10,090        |
| Peve    | Apul–Peve          | 9        | 28        | 10,090        |
| Peve    | Peve–Ptuh          | 6        | 13        | 10,090        |
| Peve    | Apul–Peve–Ptuh     | 1        | 1         | 10,090        |
| Ptuh    | Species-only       | 16,130   | 16,103    | 16,153        |
| Ptuh    | Apul–Ptuh          | 16       | 36        | 16,153        |
| Ptuh    | Peve–Ptuh          | 6        | 13        | 16,153        |
| Ptuh    | Apul–Peve–Ptuh     | 1        | 1         | 16,153        |

Overall, this workflow makes me feel better about the actual process of getting 1:1 and 1:1:1 hits, but I think this might be too stringent to really be meaningful. Since we don't expect much sequence conservation, it seems silly to apply such strict sequence similarity standards. As you can see the number of conserved sequences is extremely small to the point of what feels like not comparing sequences in a biologically relevant way.

## Next steps

-   I will do some re-runs with settings more similar to E5 deep-dive comparisons
-   Orthogroup analysis (like used in my previous post) needs further consideration as it likely allows for the expected biological reality of low conservation while still meaningfully comparing lncRNAs
-   Run blasts for the handful of conserved lncRNAs identified. These might actually be a really cool test of whether or not a small number of lncRNAs are conserved across species. Maybe we can get a hit with a previously described lncRNA of known function!

# Breakdown of code for this task

------------------------------------------------------------------------

## Data download and setup

I started by downloading species-specific lncRNA FASTAs and filtered count matrices from the `deep-dive-expression` repository, then set up directories and loaded core packages:

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

The FASTAs and count files are downloaded with a bash chunk in the Rmd, into:

``` bash
~/github/deep-dive-expression/M-multi-species/data/13.1-lncRNA-cross-species
```

------------------------------------------------------------------------

## Helper functions: length, expression, and BLAST parsing

I defined helpers to extract **lncRNA length**, **mean expression**, and to **filter BLAST hits**:

``` r
get_lnc_lengths <- function(fasta_path) {
  dna <- Biostrings::readDNAStringSet(fasta_path)
  ids_raw   <- names(dna)
  ids_clean <- sub("^[^_]+_", "", ids_raw)  # Apul_lncRNA_001 -> lncRNA_001

  tibble(
    lnc_id    = ids_clean,
    length_nt = as.integer(width(dna))
  )
}

get_lnc_expression <- function(counts_path) {
  counts <- readr::read_tsv(counts_path, show_col_types = FALSE)

  if ("Geneid" %in% names(counts)) {
    counts <- counts %>% dplyr::rename(lnc_id = Geneid)
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

For BLAST parsing, I used a function that keeps **only the single best hit per query**, after filtering on percent identity, query coverage, and e-value:

``` r
read_blast_best_hits <- function(path,
                                 pident_min = 70,
                                 qcov_min   = 0.5,
                                 evalue_max = 1e-5) {
  empty <- tibble(
    qseqid   = character(),
    sseqid   = character(),
    pident   = numeric(),
    aln_len  = integer(),
    qlen     = integer(),
    slen     = integer(),
    evalue   = numeric(),
    bitscore = numeric(),
    qcov     = numeric()
  )

  if (!file.exists(path)) {
    warning("BLAST file not found: ", path)
    return(empty)
  }

  x <- readr::read_tsv(
    path,
    col_names = c("qseqid", "sseqid", "pident",
                  "aln_len", "qlen", "slen", "evalue", "bitscore"),
    show_col_types = FALSE
  )

  if (nrow(x) == 0) {
    warning("BLAST file is empty: ", path)
    return(empty)
  }

  x %>%
    dplyr::mutate(qcov = aln_len / qlen) %>%
    dplyr::filter(
      pident >= pident_min,
      qcov   >= qcov_min,
      evalue <= evalue_max
    ) %>%
    dplyr::group_by(qseqid) %>%
    dplyr::arrange(
      dplyr::desc(bitscore),
      evalue,
      dplyr::desc(qcov),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
}
```

------------------------------------------------------------------------

## Building the core lncRNA table

Using the FASTAs and count matrices, I built a combined table with **length and expression** for each lncRNA:

``` r
species_meta <- tibble::tribble(
  ~species, ~fasta, ~counts,
  "Apul", "~/github/deep-dive-expression/M-multi-species/data/13.1-lncRNA-cross-species/Apul-lncRNA.fasta",
  "Peve", "~/github/deep-dive-expression/M-multi-species/data/13.1-lncRNA-cross-species/Peve-lncRNA.fasta",
  "Ptuh", "~/github/deep-dive-expression/M-multi-species/data/13.1-lncRNA-cross-species/Ptuh-lncRNA.fasta"
)

all_species_df <- species_meta %>%
  dplyr::mutate(
    lengths = purrr::map(fasta,  get_lnc_lengths),
    expr    = purrr::map(counts, get_lnc_expression)
  ) %>%
  dplyr::mutate(
    merged = purrr::map2(
      lengths, expr,
      ~ dplyr::left_join(.x, .y, by = "lnc_id")
    )
  )

lnc_df <- all_species_df %>%
  dplyr::select(species, merged) %>%
  tidyr::unnest(merged)
```

This `lnc_df` is the base object I keep joining onto throughout.

------------------------------------------------------------------------

## Strict homology: 1:1 and 1:1:1 using BRHs

To define **1:1 homology between species**, I used **best reciprocal hits (BRHs)**. For each species pair, I read in filtered BLAST hits and constructed BRH pairs:

``` r
apul_vs_peve <- read_blast_best_hits(file.path(blast_dir, "Apul_vs_Peve.blastn.tsv"))
peve_vs_apul <- read_blast_best_hits(file.path(blast_dir, "Peve_vs_Apul.blastn.tsv"))

apul_vs_ptuh <- read_blast_best_hits(file.path(blast_dir, "Apul_vs_Ptuh.blastn.tsv"))
ptuh_vs_apul <- read_blast_best_hits(file.path(blast_dir, "Ptuh_vs_Apul.blastn.tsv"))

peve_vs_ptuh <- read_blast_best_hits(file.path(blast_dir, "Peve_vs_Ptuh.blastn.tsv"))
ptuh_vs_peve <- read_blast_best_hits(file.path(blast_dir, "Ptuh_vs_Peve.blastn.tsv"))

get_brh_pairs <- function(q_to_s, s_to_q, q_species, s_species) {
  dplyr::inner_join(
    q_to_s %>% dplyr::select(qseqid, sseqid),
    s_to_q %>% dplyr::select(qseqid, sseqid),
    by = c("qseqid" = "sseqid", "sseqid" = "qseqid")
  ) %>%
    dplyr::transmute(
      !!q_species := qseqid,
      !!s_species := sseqid
    ) %>%
    dplyr::mutate(
      !!q_species := sub("^[^_]+_", "", .data[[q_species]]),
      !!s_species := sub("^[^_]+_", "", .data[[s_species]])
    )
}

brh_Apul_Peve <- get_brh_pairs(apul_vs_peve, peve_vs_apul, "Apul", "Peve")
brh_Apul_Ptuh <- get_brh_pairs(apul_vs_ptuh, ptuh_vs_apul, "Apul", "Ptuh")
brh_Peve_Ptuh <- get_brh_pairs(peve_vs_ptuh, ptuh_vs_peve, "Peve", "Ptuh")
```

**Three-way 1:1:1 triplets** (Apul–Peve–Ptuh) were then defined as consistent BRHs across all three pairs:

``` r
triplets_1to1to1 <- brh_Apul_Peve %>%
  dplyr::inner_join(brh_Apul_Ptuh, by = "Apul") %>%
  dplyr::inner_join(brh_Peve_Ptuh, by = c("Peve", "Ptuh"))
```

------------------------------------------------------------------------

## Per-species categories under strict criteria

For each species, every lncRNA was assigned to **exactly one category**:

-   species-only\
-   pairwise shared (Apul–Peve, Apul–Ptuh, or Peve–Ptuh)\
-   shared across all three species (Apul–Peve–Ptuh)

Example for Apul:

``` r
apul_all <- lnc_df %>%
  dplyr::filter(species == "Apul") %>%
  dplyr::pull(lnc_id) %>% unique()

apul_membership <- tibble::tibble(
  species = "Apul",
  lnc_id  = apul_all
) %>%
  dplyr::mutate(
    in_triplet = lnc_id %in% triplets_1to1to1$Apul,
    in_APEVE  = lnc_id %in% brh_Apul_Peve$Apul,
    in_APPT   = lnc_id %in% brh_Apul_Ptuh$Apul,
    category  = dplyr::case_when(
      in_triplet             ~ "Apul_Peve_Ptuh",
      in_APEVE & !in_triplet ~ "Apul_Peve",
      in_APPT  & !in_triplet ~ "Apul_Ptuh",
      TRUE                   ~ "Apul_only"
    )
  )
```

I did the same for Peve and Ptuh, then combined:

``` r
membership_all <- dplyr::bind_rows(
  apul_membership,
  peve_membership,
  ptuh_membership
) %>%
  dplyr::select(species, lnc_id, category)
```

A quick sanity check confirms that categories **sum to the FASTA totals per species**:

``` r
membership_all %>%
  dplyr::count(species, category, name = "n_lncRNAs") %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(species_total = sum(n_lncRNAs)) %>%
  dplyr::ungroup()
```

These strict categories were then joined back into `lnc_df` as `sharing_simple`.

------------------------------------------------------------------------

## Relaxed lncRNA-appropriate thresholds (post-hoc)

The strict criteria (≥70% identity, ≥50% coverage) produced very few shared lncRNAs, which is plausible for lncRNAs but too restrictive for exploring broader conservation. To loosen this **without rerunning BLAST**, I defined relaxed thresholds:

``` r
PIDENT_MIN_RELAXED <- 50
QCOV_MIN_RELAXED   <- 0.3
EVALUE_MAX_RELAXED <- 1e-5
```

Then I wrapped the BLAST reader:

``` r
read_blast_best_hits_relaxed <- function(path) {
  read_blast_best_hits(
    path,
    pident_min = PIDENT_MIN_RELAXED,
    qcov_min   = QCOV_MIN_RELAXED,
    evalue_max = EVALUE_MAX_RELAXED
  )
}
```

Using these new thresholds, I rebuilt pairwise BRHs and 1:1:1 triplets:

``` r
apul_vs_peve_r <- read_blast_best_hits_relaxed(file.path(blast_dir, "Apul_vs_Peve.blastn.tsv"))
peve_vs_apul_r <- read_blast_best_hits_relaxed(file.path(blast_dir, "Peve_vs_Apul.blastn.tsv"))

apul_vs_ptuh_r <- read_blast_best_hits_relaxed(file.path(blast_dir, "Apul_vs_Ptuh.blastn.tsv"))
ptuh_vs_apul_r <- read_blast_best_hits_relaxed(file.path(blast_dir, "Ptuh_vs_Apul.blastn.tsv"))

peve_vs_ptuh_r <- read_blast_best_hits_relaxed(file.path(blast_dir, "Peve_vs_Ptuh.blastn.tsv"))
ptuh_vs_peve_r <- read_blast_best_hits_relaxed(file.path(blast_dir, "Ptuh_vs_Peve.blastn.tsv"))

brh_Apul_Peve_r <- get_brh_pairs(apul_vs_peve_r, peve_vs_apul_r, "Apul", "Peve")
brh_Apul_Ptuh_r <- get_brh_pairs(apul_vs_ptuh_r, ptuh_vs_apul_r, "Apul", "Ptuh")
brh_Peve_Ptuh_r <- get_brh_pairs(peve_vs_ptuh_r, ptuh_vs_peve_r, "Peve", "Ptuh")

triplets_1to1to1_relaxed <- brh_Apul_Peve_r %>%
  dplyr::inner_join(brh_Apul_Ptuh_r, by = "Apul") %>%
  dplyr::inner_join(brh_Peve_Ptuh_r, by = c("Peve", "Ptuh"))
```

------------------------------------------------------------------------

## Relaxed category assignment and inspection

I repeated the per-species category assignment using the relaxed BRHs and triplets, storing labels as `category_relaxed`:

``` r
apul_membership_relaxed <- tibble::tibble(
  species = "Apul",
  lnc_id  = apul_all
) %>%
  dplyr::mutate(
    in_triplet = lnc_id %in% triplets_1to1to1_relaxed$Apul,
    in_APEVE  = lnc_id %in% brh_Apul_Peve_r$Apul,
    in_APPT   = lnc_id %in% brh_Apul_Ptuh_r$Apul,
    category_relaxed = dplyr::case_when(
      in_triplet             ~ "Apul_Peve_Ptuh",
      in_APEVE & !in_triplet ~ "Apul_Peve",
      in_APPT  & !in_triplet ~ "Apul_Ptuh",
      TRUE                   ~ "Apul_only"
    )
  )
```

Peve and Ptuh were handled analogously; then I combined:

``` r
membership_relaxed <- dplyr::bind_rows(
  apul_membership_relaxed,
  peve_membership_relaxed,
  ptuh_membership_relaxed
) %>%
  dplyr::select(species, lnc_id, category_relaxed)
```

Sanity check (again, per-species totals):

``` r
membership_relaxed %>%
  dplyr::count(species, category_relaxed, name = "n_lncRNAs") %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(species_total = sum(n_lncRNAs)) %>%
  dplyr::ungroup()
```

Finally, I joined relaxed labels back into `lnc_df`:

``` r
lnc_df_relaxed <- lnc_df %>%
  dplyr::left_join(
    membership_relaxed %>% dplyr::select(species, lnc_id, category_relaxed),
    by = c("species", "lnc_id")
  )
```

------------------------------------------------------------------------

## Visualizing relaxed patterns

To explore how relaxed homology relates to length and expression, I derived a simplified sharing label:

``` r
lnc_df_plot_relaxed <- lnc_df_relaxed %>%
  dplyr::mutate(
    species_facet = dplyr::case_when(
      species == "Apul" ~ "italic('A. pulchra')",
      species == "Peve" ~ "italic('P. evermanni')",
      species == "Ptuh" ~ "italic('P. tuahiniensis')",
      TRUE ~ species
    ),
    sharing_simple_relaxed = dplyr::case_when(
      grepl("only$", category_relaxed)     ~ "Species-only",
      category_relaxed == "Apul_Peve_Ptuh" ~ "All three",
      TRUE                                 ~ "Pairwise"
    )
  )
```

Then I plotted **lncRNA length distributions by relaxed sharing pattern**, and **length vs expression**, colored by sharing pattern, in the Rmd using `ggplot2`.
