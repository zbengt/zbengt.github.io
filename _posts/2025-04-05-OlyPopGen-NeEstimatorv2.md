---
title: "04.05.2025 Oly Pop Gen - Estimating Effective Population Size with NeEstimator v2"
author: "Zachary Bengtsson"
date:   2025-04-05
layout: post
tags:   [oly, mapping, pop gen]
---

## 🧬 NeEstimator v2 Workflow

This guide describes how to install and use **NeEstimator v2** to estimate effective population size (Ne) from a Genepop file on macOS (Apple Silicon). Note that the online host of this tool is no longer available, so you have to clone [this github repo](https://github.com/bunop/NeEstimator2.X) and build the tool on your own machine.

------------------------------------------------------------------------

### 🔧 Setup Instructions

#### 1. Install Prerequisites

**Install SDKMAN and Java 8 (Zulu JDK)**:

``` bash
curl -s "https://get.sdkman.io" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"
sdk install java 8.0.442-zulu
sdk use java 8.0.442-zulu
```

**Install Apache Ant**:

``` bash
brew install ant
```

------------------------------------------------------------------------

#### 2. Download and Build NeEstimator

Clone or download the NeEstimator v2 source code, then run:

``` bash
cd NeEstimator2.X
make
cd NeEstimator2x
ant
cd ..
mv Ne2x Ne2M
cp NeEstimator2x/dist/NeEstimator2x.jar .
```

------------------------------------------------------------------------

#### 3. Run the GUI

``` bash
java -jar NeEstimator2x.jar
```

![](https://github.com/zbengt/zbengt.github.io/blob/master/assets/img/NeEstimatorv2_GUI.png?raw=true)

Make sure `Ne2M` is in the same directory as the `.jar` file. Select your Genepop file and configure settings in the GUI. Then click **“Run Ne”**.

------------------------------------------------------------------------

## 📊 Results Summary

| Population  | LD Method Ne  | Het. Excess Ne | Molecular Coancestry Ne |
|-------------|---------------|----------------|-------------------------|
| CSMB17W.01  | \~3422 – 3926 | 119 – 153      | 13.1                    |
| CSMB18H.01  | \~490 – 539   | Infinite       | 6.8                     |
| NSFB16H.01  | \~333 – 376   | Infinite       | 19.8                    |
| NSFB18W2.01 | \~1719 – 1903 | Infinite       | 23.6                    |
| SSNB18H.01  | \~300 – 336   | 464 – 620      | 2.4                     |
| SSNB18W.01  | \~3378 – 4046 | Infinite       | 9.3                     |

------------------------------------------------------------------------

### 📝 Notes

-   LD method estimates are generally more stable and preferred for comparisons.
-   Heterozygote Excess often returned "Infinite", suggesting little to no detectable excess.
-   Molecular Coancestry consistently yielded very small Ne values, indicating possible recent bottlenecks or high relatedness.

------------------------------------------------------------------------
