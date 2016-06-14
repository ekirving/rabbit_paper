#!/usr/bin/env bash

# extract the reference genome
gunzip -k ORYCTOLAGUS/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz

# build the alignment index, using the bwtsw algorithm
bwa index -a bwtsw ORYCTOLAGUS/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa

# make all the temp directories
mkdir fastq sam bam

# unzip all the fastq files
find -name \*.fastq.gz | xargs gunzip -k

# move all the fastq files (gunzip doesn't support extracting to another directory)
find -name \*.fastq -exec mv {} ~/rabbits/fastq \;

# get all the unique codes
codes=$(ls fastq/*_*.fastq | xargs -n1 basename | awk 'BEGIN{FS="_"}{ print $1 }' | uniq)

# mapping of codes to IDs
declare -A mapping
mapping[SRR997319]=Avey36
mapping[SRR997325]=BH23
mapping[SRR997320]=FA801
mapping[SRR997322]=REX12
mapping[SRR997305]=Vau73
mapping[SRR997321]=AC100
mapping[SRR997317]=Fos6
mapping[SRR997326]=FG3
mapping[SRR997318]=Lan7
mapping[SRR997327]=A93015
mapping[SRR997304]=Fos2
mapping[SRR997323]=FL920
mapping[SRR997316]=Lan8
mapping[SRR997324]=FG4
mapping[SRR997303]=Her65

# align each of the fastq files
for code in $codes
do
    echo "Processing $code..."

    # TODO read trimming
    # Heng Li says read trimming isn't necessary with bwa mem
    # https://sourceforge.net/p/bio-bwa/mailman/message/32030167/

    # align all the fastq files (paired end reads, use 20 threads)
    bwa mem -t 20 -R "@RG\tID:$code\tSM:${mapping[$code]}\tPL:ILLUMINA" ./ORYCTOLAGUS/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa \
        ./fastq/"$code"_1.fastq ./fastq/"$code"_2.fastq 1> ./sam/"$code".sam 2> ./sam/"$code".sam.log

    # convert the sam file into a bam file
    samtools view -Shb ./sam/"$code".sam 1> ./bam/"$code".bam 2> ./bam/"$code".bam.log

done

# merge all the bam files
samtools merge merged.bam ./bam/*.bam

# sort the merged file
samtools sort merged.bam merged_sorted

# now index it
samtools index merged_sorted.bam

# -------------------------------

# remove duplicates
samtools rmdup merged_sorted.bam merged_sorted_unique.bam

# index the new bam
samtools index merged_sorted_unique.bam

# index reference sequence, so we can make the vcf
samtools faidx ORYCTOLAGUS/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa

# create a binary vcf
samtools mpileup -g -f ORYCTOLAGUS/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa \
    merged_sorted_unique.bam > merged_sorted_unique.bcf

# variant call with bcftools
bcftools view -bvcg merged_sorted_unique.bcf > variant_called.bcf

# TODO make SNP data file for dadi

# TODO tidy up old fastq, sam and bam files



samtools view -T ORYCTOLAGUS/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa -C -o yeast.cram yeast.bam


for i in `seq 1 10`;
do
    echo "Iteration $i";
    sleep 5
done

echo "Program finished" > finished.log
