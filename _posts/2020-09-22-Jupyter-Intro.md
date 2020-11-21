---
layout: post
title: Probing Jupyter
date: '2020-09-22'
categories: onboarding hematodinium Jupyter
tags: hematodinium Ubuntu Jupyter
---

**Final Results (I Think) For Hematodinium BLAST**

Again, this all refers to my project to locally blast hematodinium transcriptome data against the Swiss-prot database 

I got some help in the morning lab meeting, and learned how to get Jupyter Notebooks to accept Unix commands! Just gotta start each cell with ! or %%bash. Once I learned that, things got a lot simpler. I spent a while screwing things up in Jupyter (kept accidentally running cells by pressing shift+enter), but eventually got all the code for the hematodinium blast written out. Then I ran the blast, linked my new github repo, and boom - [final results produced!](https://github.com/afcoyle/jupyter-notebook/tree/master/hematodinium-blast-data)

As before, I ran two separate blasts - one with a value of 1 for max_hsps, the other with no specified max_hsps value. 

Now that I know how to use it, wow, Jupyter is pretty fantastic. Love how easy it is to go back and edit the code! Very impressed, excited to learn more.


