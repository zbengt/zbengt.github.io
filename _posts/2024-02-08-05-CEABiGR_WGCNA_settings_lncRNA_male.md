---
layout: post
title: 05 - WGCNA - CEABiGR WGCNA Settings with male lncRNAs
date: '2024-02-14'
categories: lncRNA results status
tags: CEABiGR WGCNA
---

TL;DR - Updated WGCNA settings to limit the number of modules. Code [here](https://github.com/zbengt/oyster-lnc/blob/main/code/12-WGCNA-lncRNA.Rmd) 

This work was completed for male lncRNAs. Ariana and I decided on WGCNA settings we think make sense with experimental design, overall sensitivity of analysis, and reasonable filtering to limit the scope of analysis. The settings under consideration are proportion of samples in which a lncRNA transcript must occur, minimum number of counts, minimum module size, and sensitivity of splitting for module generation. These decisions were made to accurately reflect treatments within the experiment, reduce noise from lowly and infrequently expressed lncRNAs, and ensure module generation was sensitive but did not generate an unmanageable number of modules. Settings are as follows:

* pOverA(0.5, 20) - Setting pOverA to this means a lncRNA must be expressed in half of the samples to be included in the analysis with a count of at least 20. Half of samples are control and half are exposed to OA conditions, so we might expect similar expression patterns based on treatment. Limiting to at least 20 removes lowly expressed lncRNAs.

* Minimum module size: 40 - This means a module must include at least 40 lncRNAs. We tested multiple min module sizes and noted an increasing trend of more modules as we increased the number to 30. There was a drop of and stabilization of modules generated between 30 and 40. So we went with 40.

* deepSplit: 3 - deepSplit influences how aggressively the clustering dendrogram is cut into separate modules. It can take on values typically in the range of 0 to 4. 1 to 3 Indicates an increasingly aggressive approach to splitting, with higher values resulting in more and smaller modules. These intermediate values allow for a balance between sensitivity and specificity in module detection. We went with 3 to ensure sensitivity in module separation, but not so much as to inflate the number of modules created.

We started with 4750 lncRNAs, removed rows in the count matrix with a sum of 0 (indicating no expression). This brought the number of lncRNAs down to 4515. After application of pOverA with the above settings, 3862 lncRNAs remained. After setting min module size and deepSplit, we were left with 31 modules of lncRNAs. This is heavily reduced from less strict parameters which generated 185 modules.




