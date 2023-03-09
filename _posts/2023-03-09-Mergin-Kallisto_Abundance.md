---
layout: post
title: Merging Kallisto Abundance Data into Count Matrix
date: '2023-03-09'
categories: genomics expression
tags: analysis
---

Ended up choosing to go with a more step by step method for merging count matrices...

Get list of all folders containing abundance data:
```r
setwd("~/github/oyster-lnc/output/01-lncRNA-kallisto")
getwd()
folders <- list.files(pattern = "S\\d+[FM]$")
print(folders)
```
Load GraphQL:
```r
install.packages("ghql")
```
Use a for loop to read in the abundance data from each folder using read.table(), and combine the data into a single count matrix using cbind():
```r
count_matrix <- NULL
for (folder in folders) {
  # set the working directory to the current folder
  setwd(paste0("/home/shared/8TB_HDD_02/zbengt/github/oyster-lnc/output/01-lncRNA-kallisto/", folder))
  
  # read in the abundance data
  data <- read.table("abundance.tsv", header = TRUE, sep = "\t")
  
  # add the abundance data to the count matrix
  count_matrix <- cbind(count_matrix, data[, "est_counts"])
}
print(count_matrix)
```
Add row names and column names to the count matrix:
```r
rownames(count_matrix) <- data[, "target_id"]
colnames(count_matrix) <- gsub(".*/", "", folders)
print(count_matrix)
```

```r
setwd("/home/shared/8TB_HDD_02/zbengt/github/oyster-lnc/output/01-lncRNA-kallisto")
write.table(count_matrix, "merged_counts.txt", sep = "\t", quote = FALSE)
```
Save the count matrix to a file using write.table():
```r
write.csv(count_matrix, file = "merged_counts.csv", row.names = TRUE)
```
