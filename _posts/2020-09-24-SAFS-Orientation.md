---
layout: post
title: Orientation!
date: '2020-09-24'
categories: onboarding hematodinium Jupyter
tags: hematodinium Ubuntu Jupyter
---

**Continuing to Analyze Hematodinium BLAST results**

Again, this all refers to my project to locally blast hematodinium transcriptome data against the Swiss-prot database 

Spent most of the morning at SAFS orientation! Chatted with the new grad students, learned some info about the program and schedule - typical orientation stuff

Given the limited time, made a surprising amount of progress today! Hit a substantial roadblock late last night - couldn't figure out a way to scrub the two files (queries and matches) of whitespace and properly run a comparison function. Realized that it'd be optimal for both the query file and match file to be the same filetype, so made new .txt files. That completely fixed the problems and let me use diff to compare the two!

Diff eliminated all lines that were identical in both files, but still left our data a bit messy. Cleaned up the existing data, and got our results - 3486 of our 6348 queries had no match!

Next task: rerun the transcriptome through DIAMOND BLASTx! Got started on it, but having some problems getting it to run. Ah well, good project for tomorrow
