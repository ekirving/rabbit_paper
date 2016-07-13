#!/usr/bin/env bash

# read in all the accession codes
readarray accessions < "$1"

for var in "${accessions[@]}"
do
#  echo "${var}"
  fastq-dump "${var//[[:space:]]/}"
done