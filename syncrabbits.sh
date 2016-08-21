#!/usr/bin/env bash

# rsync all the rabbit data
rsync -avz --partial \
           --exclude '.git' \
           --exclude 'nohup.out' \
           --exclude 'fastq' \
           --exclude 'sam' \
           --exclude 'bam' \
           --exclude '*.py*' \
           --exclude '*.old*' \
           evan@palaeoprime:~/rabbits/ ~/Dropbox/Code/rabbits/
