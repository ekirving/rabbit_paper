#!/usr/bin/env bash

# list of all th accession codes
accessions=( SRR997303 SRR997304 SRR997305 SRR997316 SRR997317 SRR997318 SRR997319 SRR997320 SRR997321 SRR997322 SRR997323 SRR997324 SRR997325 SRR997326 SRR997327 )

for var in "${accessions[@]}"
do
  echo "${var}"
  wget -m ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/"${var}"/
done