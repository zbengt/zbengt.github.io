---
layout: post
title: Getting A Grep On Things
date: '2020-09-23'
categories: onboarding hematodinium Jupyter
tags: hematodinium Ubuntu Jupyter
---

**Analyzing Hematodinium BLAST results**

Again, this all refers to my project to locally blast hematodinium transcriptome data against the Swiss-prot database 

Spent the morning at the SAFS Coffee Hour! Met most of the new grad students, learned about all the best bagel places, and got some info on how the semester would go. Overall, a great morning. 

After that, I started work on analyzing my blast results. Here are the goals:
    - Compare no_max_hsps to max_hsps blasts
        - Count lines of output files
        - Determine which queries didn't produce alignments
    - Re-run using DIAMOND BLASTx
    - Compare results of DIAMOND BLASTx and NCBI BLASTx
        - Speed
        - Number of matches
        - SPID matches
        - Anything else
It was fairly straightforward to count the output file lines, as I could just use wc -l. However, determining which queries didn't produce alignments has proved more difficult. After quite a bit of trial and error and more error, I figured out how to use grep to select just the query name. All queries begin with "TRINITY_", so I used the following command:
grep -o 'TRINITY[^ ]*' filename > query_list

If I understand it correctly, the [^ ] indicates that it searched through the line until locating the space following TRINITY, and the * meant that all characters were accepted until a space was found. 

So at this point, I have two files - a list of queries and a list of subjects. Now I just need to scrub them both of whitespace (apparently - the results are different when I run grep with -w) and cross-reference the two to isolate differences! Much easier said than done, but feeling fairly good. However, unlikely to make too much progress tomorrow, as I have orientation for much of the day
