from subprocess import check_output
import subprocess
import luigi
import os
import random
import multiprocessing

# Based on example pipleine from http://coderscrowd.com/app/public/codes/view/229

# the reference genome
GENOME = "OryCun2.0"
GENOME_URL = "ftp://ftp.ensembl.org/pub/release-84/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz"

# the list of sample codes
# SAMPLES = ['SRR997303','SRR997304','SRR997305','SRR997316','SRR997317','SRR997318','SRR997319','SRR997320',
#            'SRR997321','SRR997322','SRR997323','SRR997324','SRR997325','SRR997326','SRR997327']
SAMPLES = ['SRR997303','SRR997304']

# the URLs for the paried end fastq files
pair1_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_1.fastq.gz"
pair2_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_2.fastq.gz"

# the samtools flag for using no comression
NO_COMPRESSION = 0

# no single worker should use more than 50% of the available cores
MAX_CPU_CORES = multiprocessing.cpu_count() * 0.5

def run_cmd(cmd):
    """
    Executes the given command in a system subprocess
    """
    # TODO remove when done debugging
    print cmd

    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # run the command
    proc = subprocess.Popen(cmd,
                            shell=False,
                            # universal_newlines=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # fetch the output and error
    (stdout, stderr) = proc.communicate()

    # TODO remove when done debugging
    print vars(proc)

    # TODO do we need to delete the parent task's output if this is triggered?
    if proc.returncode:
        print stderr
        exit()

    return stdout

class Curl_Download(luigi.Task):
    """
    Downloads a remote url to a local file path using cURL
    """
    url = luigi.Parameter()
    file = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.file)

    def run(self):
        # TODO this should be piped too
        tmp = run_cmd(["curl", "-so", self.file, self.url])

        print "====== cURL Downloading file ======"

class PairedEnd_Fastq(luigi.Task):
    """
    Fetches the two paired-end FASTQ files for a given sample ID
    """
    sample = luigi.Parameter()

    def requires(self):
        global pair1_url, pair2_url
        # download the gzipped fastq files
        return [Curl_Download(pair1_url.format(sample=self.sample), "fastq/"+self.sample+"_1.fastq.gz"),
                Curl_Download(pair2_url.format(sample=self.sample), "fastq/"+self.sample+"_2.fastq.gz")]

    def output(self):
        return [luigi.LocalTarget("fastq/"+self.sample+"_1.fastq"),
                luigi.LocalTarget("fastq/"+self.sample+"_2.fastq")]

    def run(self):
        # unzip the files
        (file1, file2) = self.output()
        run_cmd(["gunzip", "-k", file1.path, file2.path]) # TODO use multithreading

        print "====== Downloaded Paired-end FASTQ ======"

class Genome_Fasta(luigi.Task):
    """
    Fetches the reference genome
    """
    genome = luigi.Parameter()
    # url = luigi.Parameter()

    def requires(self):
        global GENOME_URL
        # download the gzipped fasta file of the reference genome
        return Curl_Download(GENOME_URL, "fasta/"+self.genome+".fa.gz")

    def output(self):
        return luigi.LocalTarget("fasta/"+self.genome+".fa")

    def run(self):
        # unzip the file
        run_cmd(["gunzip", "-k", self.output().path])

        print "====== Downloaded genome FASTA ======"

class Samtools_Fasta_Index(luigi.Task):
    """
    Builds the samtools fasta index for the reference genome, needed for downstream CRAM files
    """
    genome = luigi.Parameter()

    def requires(self):
        return Genome_Fasta(self.genome)

    def output(self):
        return luigi.LocalTarget("fasta/"+self.genome+".fa")

    def run(self):
        run_cmd(["samtools",
                 "faidx",                     # index needed for CRAM files
                 "fasta/"+self.genome+".fa"]) # input file

        print "====== Built the samtools FASTA index ======"

class Bwa_Index(luigi.Task):
    """
    Builds the BWA index for the reference genome, needed for performing alignments
    """
    genome = luigi.Parameter()

    def requires(self):
        return Genome_Fasta(self.genome)

    def output(self):
        extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
        return [luigi.LocalTarget("fasta/"+self.genome+".fa."+ext) for ext in extensions]

    def run(self):
        run_cmd(["bwa",
                 "index",                     # index needed for bwa alignment
                 "-a", "bwtsw",               # algorithm suitable for mammals
                 "fasta/"+self.genome+".fa"]) # input file

        print "====== Built the BWA index ======"

class Bwa_Mem(luigi.Task):
    """
    Align the fastq files to the reference genome
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return [PairedEnd_Fastq(self.sample),
                Bwa_Index(self.genome)]

    def output(self):
        return luigi.LocalTarget("sam/"+self.sample+".sam")

    def run(self):
        # TODO add the sample name to the readgroup header
        readgroup= "@RG\\tID:{sample}\\tSM:{sample}".format(sample=self.sample)

        # perform the alignment
        sam = run_cmd(["bwa",
                       "mem",                            # align using the mem algorithm
                       "-t", cores,                      # number of cores
                       "-R", readgroup,                  # readgroup metadata
                       "fasta/"+self.genome+".fa",       # reference genome
                       "fastq/"+self.sample+"_1.fastq",  # pair 1
                       "fastq/"+self.sample+"_2.fastq"]) # pair 2

        # save the SAM file
        with self.output().open('w') as fout:
            fout.write(sam)

        print "====== Aligned the FASTQ using BWA ======"

class Convert_Sam_Sorted_Bam(luigi.Task):
    """
    Convert an uncompressed SAM file into BAM format, and sort it
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()
    compression = luigi.Parameter(default=6)

    def requires(self):
        return Bwa_Mem(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bam/"+self.sample+".bam")

    def run(self):
        # get the count of CPUs
        cores = multiprocessing.cpu_count()  # TODO should this be limited?

        # perform the SAM -> BAM conversion and sorting
        bam = run_cmd(["samtools",
                 "sort",                          # sort the reads
                 "-l", self.compression,          # level of compression
                 "-@", cores,                     # number of cores
                 "-O", "bam",                     # output a BAM file
                 "sam/"+self.sample+".sam"])      # input a SAM file

        # save the BAM file
        with self.output().open('w') as fout:
            fout.write(bam)

        print "===== Converted SAM file to BAM ======"

# class Convert_Sam_Bam(luigi.Task):
#     """
#     Convert an uncompressed SAM file into BAM format
#     """
#     sample = luigi.Parameter()
#     genome = luigi.Parameter()
#
#     def requires(self):
#         return Bwa_Mem(self.sample, self.genome)
#
#     def output(self):
#         return luigi.LocalTarget("bam/"+self.sample+".bam")
#
#     def run(self):
#         # perform the SAM -> BAM conversion
#         bam = run_cmd(["samtools",
#                        "view", "-bhS",
#                        "sam/"+self.sample+".sam"])
#
#         # save the BAM file
#         with self.output().open('w') as fout:
#             fout.write(bam)
#
#         print "===== Converting Sam file to Bam file ======"
#
# class Sort_Bam(luigi.Task):
#
#     sample = luigi.Parameter()
#
#     def requires(self):
#         return Convert_Sam_Bam(self.sample)
#
#     def run(self):
#         run_cmd(["samtools",
#                  "sort",
#                  "bam/"+self.sample+".bam",
#                  "bam/"+self.sample+".sorted"])
#
#         print "===== Sorting Bam ======="

class Convert_Bam_Cram(luigi.Task):
    """
    Convert a BAM file into the smaller CRAM format
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return [Convert_Sam_Sorted_Bam(self.sample, self.genome, NO_COMPRESSION),
                Samtools_Fasta_Index(self.genome)]

    def output(self):
        return luigi.LocalTarget("cram/"+self.sample+".cram")

    def run(self):
        # get the count of CPUs
        cores = multiprocessing.cpu_count()  # TODO should this be limited?

        # perform the SAM -> BAM conversion
        cram = run_cmd(["samtools",
                        "view",
                        "-@", cores,                       # number of cores
                        "-C",                              # output a CRAM file
                        "-o", "cram/"+self.sample+".cram", # output location
                        "bam/"+self.sample+".bam"])        # input BAM file

        # save the CRAM file
        with self.output().open('w') as fout:
            fout.write(cram)

        print "===== Converted BAM file to CRAM ======"

class Index_Cram(luigi.Task):
    """
    Create an index for the CRAM file for fast random access
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("cram/" + self.sample + ".cram.crai")

    def requires(self):
        return Convert_Bam_Cram(self.sample, self.genome)

    def run(self):
        tmp = run_cmd(["samtools",
                       "index",
                       "cram/" + self.sample + ".cram"])

        print "==== Indexing CRAM ===="

class Samtools_MPileup(luigi.Task):
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Index_Cram(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("pileup/"+self.sample+".pileup")

    def run(self):
        pileup = run_cmd(["samtools",
                          "mpileup",                             # output a pileup file
                          "-o", "pileup/"+self.sample+".pileup", # output location
                          "-f", "fasta/"+self.genome+".fa",      # reference genome
                          "cram/"+self.sample+".cram"])          # input CRAM file

        # TODO can't use output pipe because it fails with error "IOError: [Errno 22] Invalid argument"

        # save the pileup file
        # with self.output().open('w') as fout:
        #     fout.write(pileup)

        print "===== Converted CRAM file to pileup ======="

# class Call_Variant(luigi.Task):
#
#     sample = luigi.Parameter()
#
#     def requires(self):
#         return Index_Bam(self.sample)
#
#     def output(self):
#         return luigi.LocalTarget("bcf/"+self.sample+".bcf")
#
#     def run(self):
#         tmp = run_cmd(["samtools",
#                        "mpileup",
#                        "-u",
#                        "-f",
#                        genome,
#                        "bam/"+self.sample+".sorted.bam"])
#
#         with self.output().open("w") as bcf_file:
#             bcf_file.write(tmp)
#
#         print "=====  Calling Variants ======="
#
#
# class Convert_Bcf_Vcf(luigi.Task):
#
#     sample = luigi.Parameter()
#
#     def requires(self):
#         return Call_Variant(self.sample)
#
#     def output(self):
#         return luigi.LocalTarget("vcf/"+self.sample+".vcf")
#
#     def run(self):
#         tmp = run_cmd(["bcftools",
#                        "view",
#                        "-v",
#                        "-c",
#                        "-g",
#                        "bcf/"+self.sample+".bcf"])
#
#         with self.output().open("w") as vcf_file:
#             vcf_file.write(tmp)
#
#         print "===== Creating VCF ======="

class Custom_Genome_Pipeline(luigi.Task):
    """
    Run all the samples through the pipeline
    """
    def requires(self):
        for sample in SAMPLES:
            yield Samtools_MPileup(sample, GENOME)

if __name__=='__main__':
    luigi.run()