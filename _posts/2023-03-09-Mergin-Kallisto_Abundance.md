---
layout: post
title: Merging Kallisto Abundance Data into Count Matrix
date: '2023-03-09'
categories: genomics expression
tags: analysis
---

Ended up choosing to go with a more step by step method for merging count matrices...

```{r}
### set folder path to location of your kallisto outputs
folder_path <- "~/github/oyster-lnc/output/01-lncRNA-kallisto"
folders <- list.files(path = folder_path, pattern = "S\\d+[FM]$")
print(folders)
```
```{r}
install.packages("ghql")
library("ghql")
```


```{r}
count_matrix <- NULL
for (folder in folders) {
  # specify the full file path to the abundance data
  file_path <- paste0("/home/shared/8TB_HDD_02/zbengt/github/oyster-lnc/output/01-lncRNA-kallisto/", folder, "/abundance.tsv")
  
  # read in the abundance data
  data <- read.table(file_path, header = TRUE, sep = "\t")
  
  # add the abundance data to the count matrix
  count_matrix <- cbind(count_matrix, data[, "est_counts"])
}
print(count_matrix)
```

```{r}
rownames(count_matrix) <- data[, "target_id"]
colnames(count_matrix) <- gsub(".*/", "", folders)
print(count_matrix)
```

```{r}
write.csv(count_matrix, file = "merged_kallisto.csv", row.names = TRUE)
```