---
layout: post
title: Linear Modeling, for real
date: '2021-09-09'
categories: hematodinium
tags: linear modeling
---

Longtime readers will remember my post from June on [using linear modeling to model Hematodinium infection](https://afcoyle.github.io/2021-06-09-QERM_514_project/). In case you don't, I'll give you a quick recap

As part of the QERM 514 class I took, I created a linear model from the survey data that Pam Jensen sent over. These data were Southeast AK Tanner crab surveys from 2007 to 2012, and I focused on modeling _Hematodinium_ infection status. My model found the following linkages to _Hematodinium_ infection.

- Shell condition: strong linkage. Low for fresh-shell crabs, high for new-shell, then decline as shell age increases

- Depth: Infection rates higher for crabs in shallower waters

- Carapace width: Larger crabs were more likely to be infected

- Julian day: Included in final model, but likely unimportant (p-val = 0.85)

- Black Mat model: Included in final model, not judged to be significant (p-val = 0.54). However, despite sampling ~800 crabs with Black Mat, zero had both Black Mat and _Hematodinium_ infections. Indicates that we just didn't have enough data for this to be significant. Also definitely interesting - may want to see if fungal infections exclude dinoflagellate infections and vice-versa.

However, there's a big ol' caveat to all this. This model failed the Hosmer-Lemeshow test, indicating that there's variance in the data that isn't accounted for in the model. This menas that we're missing at least one important predictor variable, and the conclusions of the model should be taken with caution.

Basically, the results from my model itself were really neat! But they were also sadly incomplete. To draw better conclusions, I needed more data with more variables. So I got ahold of two possible data sources.

I'll start with the less interesting one - the NOAA EBS trawl data from 2014-2019. This is decent, as it does contain temperature data (a likely possibility for the source of our missing variation). However, there are some issues with it. 
- Latitude, longitude, depth, and Julian day are all correlated. This is because the survey follows the same basic pattern from year to year, the overall survey proceeds in a northwesterly direction, and depth increases as you move west. Since these are all tightly correlated, proceeding means a) ditching lat/long in favor of a categorical location, and b) tossing out Julian day.

- Crab species is also correlated with lat/long and location. This makes sense - snow crab live further north than Tanner crab. However, it also means that we need to either toss out location as a variable (impractical) or make separate models for each species. I'm currently proceeding with that now.

Now, the more interesting data. I got in touch with the ADF&G crab team in Juneau, and they sent over all survey data from their red king crab and Tanner crab surveys going all the way back to 1978!! This is an immense amount of data. However, it has both upsides and downsides. 

- Upside: Tons and tons of data. Like, over four decades worth of data. Also contains temperature data!

- Downside: That temperature data is all contained within separate CSV files, rather than being part of the sheet. That is...not optimal, to say the least. But hey, at least we have it!

Next step: Try to figure out how to extract the temperature data in a reasonably-organized manner and merge it with the main sheet. It's gonna be tricky - we have dozens and dozens of separate CSV files - but it seems doable at least!