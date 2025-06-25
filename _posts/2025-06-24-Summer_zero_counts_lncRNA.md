---
title: "06.24.2025 0 lncRNA counts in E5 deep-dive-expression matrices due to FeatureCounts settings"
date: 2025-06-24
layout: post
categories: lncRNA bioinformatics e5
---

TL;DR - I was able to resolve the issue of all 0 count rows in the Ptuh count matrix by adding options for including multi-overlap and multi-mapping. These deviate from the default settings of FeatureCounts and I wonder if we just want to eliminate the 25 lncRNAs unable to pass these default settings. We can, however, proceed with all counts.

### FeatureCounts

FeatureCounts is fundamentally a read (or fragment) counter, not an “expression estimator,” so it doesn’t apply any hidden rounding or expression‐level filtering at count time. Here’s what actually happens under the hood:

### Whole‐read (fragment) counting

-   By default featureCounts asks “does this fragment overlap a feature by at least one base?” (the --minOverlap 1 rule). If any base of either mate in a paired‐end fragment overlaps your GTF feature, you get 1 count for that fragment subread.sourceforge.net biostars.org .
-   There is no notion of “fractional read = 0.2 counts” unless you explicitly turn on the --fraction option (see below).

### No built-in expression threshold

-   If zero fragments overlap a feature, featureCounts simply reports 0.
-   There is no further “filter” at the counting step to drop low‐coverage features or to require a minimum count threshold. Everything in your GTF shows up, zero or more.

### Handling multi-overlaps and multi-mapping

-   By default, fragments that overlap more than one feature are considered ambiguous and are not counted at all (you get no partial credit) bioconductor.org .
-   If you really want each overlapping feature to get a hit, you can add -O (allowMultiOverlap). If you also add --fraction, each feature receives 1/x where x = number of overlapping features.
-   Similarly, multi-mapping reads are ignored by default; -M (countMultiMappingReads) and --fraction can be combined to assign fractional weights to each alignment.

#### Re-ran Ptuh with FeatureCounts settings for multi-overlap and fractions

```{bash}
/home/shared/subread-2.0.5-Linux-x86_64/bin/featureCounts \
  -T 24 \
  -O \
  --fraction \
  -a ../output/01.6-Ptuh-lncRNA-pipeline/Ptuh_lncRNA_for_fc.gtf \
  -o ../output/01.6-Ptuh-lncRNA-pipeline/Ptuh_counts_fractional.txt \
  -t lncRNA \
  -g gene_id \
  -p \
  ../output/01.6-Ptuh-lncRNA-pipeline/*sorted.bam

```

`-O` allows fragments to count toward every overlapping feature. `--fraction` assign 1/x count if a fragment overlaps x features.

Ran these options together and separately with the **result of bringing the original number of full 0 entries in the count matrix down from 25 to 21.**

#### Re-run FeatureCounts with even more settings for multi-mapping and pair-ended reads

```{bash}
/home/shared/subread-2.0.5-Linux-x86_64/bin/featureCounts \
  -T 24 \
  -M \
  -O \
  --fraction \
  -a ../output/01.6-Ptuh-lncRNA-pipeline/Ptuh_lncRNA_for_fc.gtf \
  -o ../output/01.6-Ptuh-lncRNA-pipeline/Ptuh_counts_fractional_overlap-multi-map_pair-end.txt \
  -t lncRNA \
  -g gene_id \
  -p \
  ../output/01.6-Ptuh-lncRNA-pipeline/*sorted.bam
  
```

Ran with `-O` and `--fraction` as well as `-M` for multi-mapped features with the **result of bringing the original number of full 0 entries in the count matrix down from 21 to 0.**

#### Summary

The addition of multi-overlap and multi-mapped features to obtain counts at the FeatureCounts step resulted in counts being attributed to those rows with 0s across all samples. This is the desired outcome, but I wonder if we even want to include these counts for overly repetitive lncRNAs and those that overlap more than one feature. We could leave the FeatureCount settings at the default and eliminate all 0 rows as a matter of filtering.

Summary of what each additional option does:

| Option       | Long form                  | Description                                                                                       | Notes / Combinations                                                       |
|-------------|-------------|---------------------------|--------------------|
| `-O`         | `--allowMultiOverlap`      | **Allow multi-feature overlap.**<br>Counts a fragment for every feature it overlaps.              | Use when features overlap; each overlapping feature receives a full count. |
| `-M`         | `--countMultiMappingReads` | **Include multi-mapping reads.**<br>Counts reads that map to multiple genomic locations.          | Without `--fraction`, each alignment receives a full count.                |
| `--fraction` | —                          | **Fractional assignment.**<br>Divides a fragment’s count by its number of overlaps or alignments. | Must be combined with `-O` and/or `-M` to avoid inflating total counts.    |
