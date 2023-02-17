---
layout: post
title: DESeq2 Differential Expression Analysis Pipeline
date: '2023-02-16'
categories: genomics expression
tags: analysis
---

# Basic code to build off of...


### Load your expression data into R as a matrix or data frame

```r
### txt
expression_data <- read.table("expression_data.txt", header=TRUE, sep="\t")

### csv
expression_data <- read.csv("expression_data.csv")
```

### Normalize the expression data to account for technical variability and scale differences

```r
library(DESeq2)
expression_dds <- DESeqDataSetFromMatrix(countData = expression_data,
                                          colData = sample_info,
                                          design = ~ condition)
```

Note that the DESeqDataSetFromMatrix function also requires a data frame sample_info with information about the samples, such as the experimental condition, and a formula specifying the design of the experiment. The DESeq2 package provides several normalization methods, including the default "SizeFactors" method.

### Perform statistical tests to identify differentially expressed genes: DESeq() function

```r
dds <- DESeq(expression_dds)
```

The DESeq2 package provides a likelihood-ratio test for the negative binomial distribution, which is the default test for identifying differentially expressed genes.

After the statistical tests have been performed, the resulting p-values can be adjusted for multiple testing to control the false discovery rate (FDR) using methods such as the Benjamini-Hochberg procedure. The adjusted p-values can then be used to identify differentially expressed genes. The limma, DESeq2, and EdgeR packages all offer functions for multiple testing correction.

It is important to understand the assumptions and limitations of each statistical test and to choose the most appropriate test based on the characteristics of your data and the question you are trying to answer.

### Correct for multiple testing to control for false positive results.

Multiple testing correction is an important step in differential expression analysis to control for false positive results due to the large number of hypothesis tests performed. Multiple testing correction adjusts the p-values obtained from the statistical tests to account for the number of tests performed and the probability of observing false positive results.

```r
res <- results(dds)
res <- res[which(res$padj<0.05),]
```

Here, res$padj is a vector of adjusted p-values, and the resulting data frame res only includes rows for genes with adjusted p-values less than 0.05.

### Visualize and interpret the results
Volcano plots: Volcano plots display the log2 fold-change of each gene on the x-axis and its corresponding -log10 adjusted p-value on the y-axis. Genes that are significantly differentially expressed are plotted as dots, and the dots are colored based on the direction of the log2 fold-change.

Heatmaps: Heatmaps display the expression levels of a subset of genes as a matrix, where each row represents a gene and each column represents a sample. The expression levels are represented as colors, with red indicating high expression and green indicating low expression.

Principal component analysis (PCA) plots: PCA plots display the first two principal components of the expression data, which summarize the variability in the data. Samples with similar expression profiles will cluster together in the plot, and the first two principal components explain the most variability in the data.

Clustered heatmaps: Clustered heatmaps display the expression levels of a subset of genes as a matrix, where each row represents a gene and each column represents a sample. The expression levels are represented as colors, with red indicating high expression and green indicating low expression. The rows and columns are also clustered based on the similarity of the expression profiles.

```r
### volcano plot

library(ggplot2)
ggplot(results, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color=ifelse(log2FoldChange>0 & padj<0.05, "red", "black")), alpha=0.5) +
  scale_color_manual(values=c("red", "black")) +
  theme_classic() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 p-value") +
  ggtitle("Volcano Plot")
```

```r
### heatmap

### counts matrix
counts_matrix <- counts(dds, normalized=TRUE)
de_genes <- rownames(results[which(results$padj < 0.05 & abs(results$log2FoldChange) > 1),])
de_counts <- counts_matrix[de_genes,]

### transform data
rld <- rlog(dds, blind=FALSE)
rld_de_genes <- rld[de_genes,]

### create heatmap
library(pheatmap)
pheatmap(assay(rld_de_genes), scale="row", show_rownames=FALSE, fontsize=8, main="Heatmap of Differentially Expressed Genes")
```

```r
### PCA

### normalize
rld <- rlog(dds, blind=FALSE)

### plot
plotPCA(rld, intgroup="condition", returnData=TRUE)
```













