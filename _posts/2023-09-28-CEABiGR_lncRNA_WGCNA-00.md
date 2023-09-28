---
layout: post
title: CEABiGR lncRNA - WGCNA-00
date: '2023-09-27'
categories: bioinformatics CEABiGR
tags: CEABiGR WGCNA
---

Running individual WGCNA analyses on male and female sets of oysters. Initial code can be found [here](https://github.com/zbengt/oyster-lnc/blob/main/code/05-WGCNA-edgeR.Rmd).
General notes about first run:
* Using signed networks (those that take into account the direction of expression as well as the value of the count - up or down)
* Blockwise clustering for module formation, usually good for large datasets, probably not necessary since there isn't too much data
* Current status: Switched from this analysis at module merging stage since none of the modules merged with a stringent similarity setting

Male modules and clustering:

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/WGCNA-00-male-modules.png?raw=true)

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/WGCNA-00-male-clustering.png?raw=true)

Female modules and clustering:

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/WGCNA-00-female-modules.png?raw=true)

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/WGCNA-00-female-clustering.png?raw=true)

Attempt at merging, showing males only in this case. Module similarity did not meet threshold to merge any of the many modules:

![image](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/oyster-lnc/WGCNA-00-male-merge.png?raw=true)

Paused this to attempt Ariana's methodology which has more stringent filtering steps.