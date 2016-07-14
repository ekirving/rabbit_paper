#!/usr/bin/env bash

# read in all the accession codes
readarray accessions < "$1"

for var in "${accessions[@]}"
do
  printf "Downloading ${var}... "

  # dump as fastq files, with paired-end reads in separate files
  fastq-dump --gzip --split-files "${var//[[:space:]]/}"

  printf "done!\n"
done