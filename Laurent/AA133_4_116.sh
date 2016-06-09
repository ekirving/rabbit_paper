#!/usr/bin/env bash
/usr/local/bin/AdapterRemoval --file1 ../../screen/run_CPN_11052016/ftp.dna.ku.dk/20160503_hiseq4b/Project_PalaeoBARN/Sample_TOG_4UNA_OXF20_AA133/TOG_4UNA_OXF20_AA133_GCAGAG_L008_R1_006.fastq.gz --trimns --output1 AA133/AA133_OXF20_L008_006-adRm.fq.gz --gzip --basename AA133/AA133_OXF20_L008_006-adRm.fq.gz
/usr/local/bin/AdapterRemoval --file1 ../../screen/run_CPN_11052016/ftp.dna.ku.dk/20160503_hiseq4b/Project_PalaeoBARN/Sample_TOG_4UNA_OXF20_AA133/TOG_4UNA_OXF20_AA133_GCAGAG_L008_R1_004.fastq.gz --trimns --output1 AA133/AA133_OXF20_L008_004-adRm.fq.gz --gzip --basename AA133/AA133_OXF20_L008_004-adRm.fq.gz
/usr/local/bin/AdapterRemoval --file1 ../../screen/run_CPN_11052016/ftp.dna.ku.dk/20160503_hiseq4b/Project_PalaeoBARN/Sample_TOG_4UNA_OXF20_AA133/TOG_4UNA_OXF20_AA133_GCAGAG_L008_R1_007.fastq.gz --trimns --output1 AA133/AA133_OXF20_L008_007-adRm.fq.gz --gzip --basename AA133/AA133_OXF20_L008_007-adRm.fq.gz
/usr/local/bin/AdapterRemoval --file1 ../../screen/run_CPN_11052016/ftp.dna.ku.dk/20160503_hiseq4b/Project_PalaeoBARN/Sample_TOG_4UNA_OXF20_AA133/TOG_4UNA_OXF20_AA133_GCAGAG_L008_R1_002.fastq.gz --trimns --output1 AA133/AA133_OXF20_L008_002-adRm.fq.gz --gzip --basename AA133/AA133_OXF20_L008_002-adRm.fq.gz
/usr/local/bin/AdapterRemoval --file1 ../../screen/run_CPN_11052016/ftp.dna.ku.dk/20160503_hiseq4b/Project_PalaeoBARN/Sample_TOG_4UNA_OXF20_AA133/TOG_4UNA_OXF20_AA133_GCAGAG_L008_R1_003.fastq.gz --trimns --output1 AA133/AA133_OXF20_L008_003-adRm.fq.gz --gzip --basename AA133/AA133_OXF20_L008_003-adRm.fq.gz
/usr/local/bin/AdapterRemoval --file1 ../../screen/run_CPN_11052016/ftp.dna.ku.dk/20160503_hiseq4b/Project_PalaeoBARN/Sample_TOG_4UNA_OXF20_AA133/TOG_4UNA_OXF20_AA133_GCAGAG_L008_R1_005.fastq.gz --trimns --output1 AA133/AA133_OXF20_L008_005-adRm.fq.gz --gzip --basename AA133/AA133_OXF20_L008_005-adRm.fq.gz
/usr/local/bin/AdapterRemoval --file1 ../../screen/run_CPN_11052016/ftp.dna.ku.dk/20160503_hiseq4b/Project_PalaeoBARN/Sample_TOG_4UNA_OXF20_AA133/TOG_4UNA_OXF20_AA133_GCAGAG_L008_R1_001.fastq.gz --trimns --output1 AA133/AA133_OXF20_L008_001-adRm.fq.gz --gzip --basename AA133/AA133_OXF20_L008_001-adRm.fq.gz
gunzip -c  AA133/AA133_OXF20_L008_006-adRm.fq.gz AA133/AA133_OXF20_L008_004-adRm.fq.gz AA133/AA133_OXF20_L008_007-adRm.fq.gz AA133/AA133_OXF20_L008_002-adRm.fq.gz AA133/AA133_OXF20_L008_003-adRm.fq.gz AA133/AA133_OXF20_L008_005-adRm.fq.gz AA133/AA133_OXF20_L008_001-adRm.fq.gz
    | bwa aln -l 1024 -n 0.03 -t 20 /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa - 1> AA133/AA133_OXF20.sai 2> AA133/AA133_OXF20.log
gunzip -c  AA133/AA133_OXF20_L008_006-adRm.fq.gz AA133/AA133_OXF20_L008_004-adRm.fq.gz AA133/AA133_OXF20_L008_007-adRm.fq.gz AA133/AA133_OXF20_L008_002-adRm.fq.gz AA133/AA133_OXF20_L008_003-adRm.fq.gz AA133/AA133_OXF20_L008_005-adRm.fq.gz AA133/AA133_OXF20_L008_001-adRm.fq.gz
    | bwa samse -r '@RG\tID:OXF20\tLB:OXF20\tPL:ILLUMINA\tSM:AA133' /home/laurent/REF/SUS/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa AA133/AA133_OXF20.sai -
    | samtools view -Shu - 1> AA133/AA133_OXF20.bam 2>> AA133/AA133_OXF20.log

samtools sort  AA133/AA133_OXF20.bam AA133/AA133_sorted
samtools index AA133/AA133_sorted.bam
samtools rmdup -s AA133/AA133_sorted.bam AA133/AA133_rmdup.bam 2>&1 >/dev/null | cut -d' ' -f 6 >  AA133/AA133_dup.txt
samtools index AA133/AA133_rmdup.bam
samtools view AA133/AA133_rmdup.bam | perl /media/raid5/laurent/full_run_results/bin/compute_l.pl > AA133/AA133_RC.txt
samtools view -q 20 -b AA133/AA133_rmdup.bam | bedtools genomecov -ibam stdin | grep -v genome > AA133/AA133.cov
