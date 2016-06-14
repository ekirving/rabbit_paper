#!/usr/bin/env bash

# extract the reference genome
gunzip -k ORYCTOLAGUS/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz

# build the alignment index, using the bwtsw algorithm
bwa index -a bwtsw ORYCTOLAGUS/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa

# find all the fastq.gz files
FASTQS=$(find -name \*.fastq.gz)

# align each of the fastq files
for file in $FASTQS
do
    # get the sample code from the filename
    code=$(basename "$file" .fastq.gz)

    echo "Processing $code..."

    # TODO trim ends

    # unzip the file and pipe the output to bwa to align it (use 20 threads)
    gunzip -c "$file" | bwa mem -t 20 Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa - 1> "$code".sam 2> "$code".sam.log

    # TODO convert the sam file into a bam file

    # TODO samtools to merge/index/sort

done


 find -name \*SRR997319_1.fastq.gz | xargs gunzip -c | bwa mem -t 20 Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa - 1> PRJNA221358.sam 2> PRJNA221358.sam.log

 samtools view -hb PRJNA221358.sam 1> PRJNA221358.bam 2>> PRJNA221358.bam.log


 # --------------------------------------------


 bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' <ref.fa> <read1.fa> <read1.fa> > lane.sam
          -r '@RG\tID:OXF20\tLB:OXF20\tPL:ILLUMINA\tSM:AA133'


bwa mem -t 20 -R '@RG\tID:RABBITS\tLB:RABBITS\tPL:ILLUMINA\tSM:SRR997327' Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa SRR997327_1.fastq SRR997327_2.fastq > SRR997327.sam