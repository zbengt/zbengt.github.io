---
title: "07.07.2025 Converting genind to genpop file for updated NeEstimator run"
date: 2025-07-07
layout: post
categories: oly-genetics bioinformatics
---

TL;DR - Had to re-run Oly RAD-seq data for effective population size after outliers were filtered out. Had to convert from genind to GENEPOP, which was much more annoying than I expected.

## Converting a `genind` Object to GENEPOP Format

This code shows how to export a **genind** object (`oly.genind`) to a GENEPOP‐format file for use in NeEstimator2.

### 1. Install and load `graph4lg`

``` r
if (!requireNamespace("graph4lg", quietly = TRUE)) {
  install.packages("graph4lg")
}
library(graph4lg)
```

### 2. Load the `.Rdata` File and Identify the `genind` Object

``` r
load("/Users/zacharybengtsson/Downloads/oly.no23.Rdata")
ls()
# [1] "oly.genind"  # your genind object
```

### 3. Export to GENEPOP Format (`.txt`)

``` r
genind_to_genepop(
  oly.genind,
  output = "/Users/zacharybengtsson/Desktop/NeEstimator2.X/oly_no23.txt"
)
```

### 4. (Optional) Rename to `.gen` Extension

``` r
file.rename(
  "/Users/zacharybengtsson/Desktop/NeEstimator2.X/oly_no23.txt",
  "/Users/zacharybengtsson/Desktop/NeEstimator2.X/oly_no23.gen"
)
```

### 5. (Optional) Fix Missing‐Data Padding with `sed`

In case missing genotypes are six-digit codes (`000000`), collapse them to four digits (`0000`):

``` bash
sed 's/000000/0000/g' \
  /Users/zacharybengtsson/Desktop/NeEstimator2.X/oly_no23.gen \
  > /Users/zacharybengtsson/Desktop/NeEstimator2.X/oly_no23_fixed.gen
```

Now `/Users/.../oly_no23_fixed.gen` is ready for NeEstimator2.

### 6. Run in NeEstimator2

The file is now ready for NeEstimator2 which can be run with the [linked instructions in my lab notebook](https://zbengt.github.io/2025-04-05-OlyPopGen-NeEstimatorv2/).

The result is update effective population sizes (Ne) for wild and hatchery individuals in the North, Central, and South Sound.

Here’s the table in Markdown format:

| Population | N   | LD Ne (Pcrit 0.05–0.01; parametric CI; jackknife CI) | HE Nb (range; upper CI) | MC Nb (95% CI)   |
|-----------|-----------|----------------------------|-----------|-----------|
| CSMB17W    | 95  | 2 994–3 466 (2 476–4 383; 1 417–∞)                   | 219–447 (upper ∞)       | 12.4 (10.1–15.0) |
| CSMB18H    | 95  | 485–531 (469–548; 278–1 590)                         | ∞                       | 7.1 (6.0–8.4)    |
| NSFB16H    | 89  | 331–374 (323–384; 239–586)                           | ∞                       | 22.5 (16.4–29.7) |
| NSFB18W2   | 95  | 1 716–1 880 (1 540–2 115; 1 301–2 725)               | ∞                       | 32.4 (23.0–43.3) |
| SSNB18H    | 94  | 298–333 (291–341; 200–623)                           | 1 051–2 797 (upper ∞)   | 2.6 (2.3–2.9)    |
| SSNB18W    | 93  | 3 338–4 054 (2 714–5 361; 1 698–∞)                   | ∞                       | 7.6 (6.5–8.6)    |
