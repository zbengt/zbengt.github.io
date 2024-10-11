---
layout: post
title: 10.10.24 rclone for Google Drive
date: '2024-10-10'
categories: oly
tags: popgen download
---

Mac provided Olympia oyster RADSeq data via Google Drive. BAMs are corrupted when compressed by drive's go to zipping tool, so had to use rclone to get the data on my computer.

### Download rclone

[rclone installation page](https://rclone.org/downloads/)

My Mac has M2 so installed the ARM - 64 Bit version

### Configure rclone

``` bash
rclone config
```

Follow the prompts: \* Select n for a new remote. \* Name the remote (e.g., mygdrive). \* Choose drive as the storage type. \* Follow the prompts to authenticate with your Google account. It will provide a link for you to authorize rclone access to your Google Drive. Paste the authorization code back into the terminal.

[Step by step instructions](https://rclone.org/drive/#making-your-own-client-id) for making your own Client ID

### Copy Shared with Me Folder

```bash
rclone copy mygdrive:"SharedFolderName" /local/destination --drive-shared-with-me
```

"SharedFolderName" can be replaced with file path starting with the shared folder. This is good for only downloading target subfolders. I used "OlyRAD_6plates/SEDNA_Olympia_oyster_data/bowtie2".

### Hosting data on Gannet

Working with data a little on my own computer but will eventually host BAMs on Gannet for easy linking.
