import luigi
import multiprocessing
import numpy
import os
import json
import subprocess
from collections import OrderedDict

# import the custom vcf parser
from vcfparser import *

# Based on example pipleine from http://coderscrowd.com/app/public/codes/view/229

# the reference genome
GENOME = "OryCun2.0"
GENOME_URL = "ftp://ftp.ensembl.org/pub/release-84/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz"

SAMPLES = OrderedDict()

# wild population, w/ accession codes and sample ID
SAMPLES['SRR997319'] = 'Avey36' # Aveyron
SAMPLES['SRR997317'] = 'Fos6'   # Fos-su-Mer
SAMPLES['SRR997304'] = 'Fos2'   # Fos-su-Mer
SAMPLES['SRR997303'] = 'Her65'  # Herauld
SAMPLES['SRR997318'] = 'Lan7'   # Lancon
SAMPLES['SRR997316'] = 'Lan8'   # Lancon
SAMPLES['SRR997305'] = 'Vau73'  # Vaucluse

# record the index of the last wild sample
WILD_THRESHOLD = len(SAMPLES)

# domestic population, w/ accession codes and sample ID
SAMPLES['SRR997325'] = 'BH23'    # Belgian hare
SAMPLES['SRR997320'] = 'FA801'   # Champagne d'argent
SAMPLES['SRR997321'] = 'AC100'   # Angora
SAMPLES['SRR997327'] = 'A93015'  # Angora
SAMPLES['SRR997323'] = 'FL920'   # French lop
SAMPLES['SRR997326'] = 'FG3'     # Flemish giant
SAMPLES['SRR997324'] = 'FG4'     # Flemish giant
SAMPLES['SRR997322'] = 'REX12'   # Rex

# the URLs for the paried end fastq files
PAIR1_URL = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_1.fastq.gz"
PAIR2_URL = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_2.fastq.gz"

# the samtools flag for BAM file comression
DEFAULT_COMPRESSION = 6

# TODO find out how many workers there are
# no single worker should use more than 50% of the available cores
MAX_CPU_CORES = int(multiprocessing.cpu_count() * 0.5)


def run_cmd(cmd):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :return: The stdout stream
    """

    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # TODO write all bash commands to a script and log it
    print cmd

    # run the command
    proc = subprocess.Popen(cmd,
                            shell=False,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # fetch the output and error
    (stdout, stderr) = proc.communicate()

    # TODO remove when done debugging
    print vars(proc)

    # bail if something went wrong
    if proc.returncode:
        print stderr
        exit()

    return stdout

def unzip_file(gzip):
    """
    Unzip a gzipped file, using multi-threading when available

    :param gzip: The path to the gzip file
    :return: The uncompressed data stream
    """
    try:
        # use unpigz for multithreaded unzipping (if installed)
        return run_cmd(["unpigz",
                        "-c",                # output to stdout
                        "-p", MAX_CPU_CORES, # use the maximum cores
                        gzip])               # gzip file

    except OSError as e:
        # if it's not installed
        if e.errno == os.errno.ENOENT:

            # use single-threaded gunzip
            return run_cmd(["gunzip",
                            "-c",            # output to stdout
                            gzip])           # gzip file
        else:
            # escalate the exception
            raise e


class Curl_Download(luigi.Task):
    """
    Downloads a remote url to a local file path using cURL
    """
    url = luigi.Parameter()
    file = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.file)

    def run(self):
        # download the file
        stdout = run_cmd(["curl",
                          "-s",       # download silently
                          self.url])  # from this url

        # save the file
        with self.output().open('w') as fout:
            fout.write(stdout)

        print "====== cURL Downloading file ======"


class Sample_PairedEnd_Fastq(luigi.Task):
    """
    Fetches the two paired-end FASTQ files for a given sample ID
    """
    sample = luigi.Parameter()

    def requires(self):
        global PAIR1_URL, PAIR2_URL
        # download the gzipped fastq files
        return [Curl_Download(PAIR1_URL.format(sample=self.sample), "fastq/" + self.sample + "_1.fastq.gz"),
                Curl_Download(PAIR2_URL.format(sample=self.sample), "fastq/" + self.sample + "_2.fastq.gz")]

    def output(self):
        return [luigi.LocalTarget("fastq/" + self.sample + "_1.fastq"),
                luigi.LocalTarget("fastq/" + self.sample + "_2.fastq")]

    def run(self):
        # fetch the downloaded files
        (pair1, pair2) = self.requires()

        # unzip them
        fastq1 = unzip_file(pair1.file)
        fastq2 = unzip_file(pair2.file)

        with self.output()[0].open('w') as fout1, self.output()[1].open('w') as fout2:
            fout1.write(fastq1)
            fout2.write(fastq2)

        print "====== Downloaded Paired-end FASTQ ======"

class Reference_Genome_Fasta(luigi.Task):
    """
    Fetches the reference genome
    """
    genome = luigi.Parameter()
    # url = luigi.Parameter()

    def requires(self):
        global GENOME_URL
        # download the gzipped fasta file of the reference genome
        return Curl_Download(GENOME_URL, "fasta/" + self.genome + ".fa.gz")

    def output(self):
        return luigi.LocalTarget("fasta/" + self.genome + ".fa")

    def run(self):
        # unzip the file
        fasta = unzip_file(self.requires().file)

        with self.output().open('w') as fout:
            fout.write(fasta)

        print "====== Downloaded genome FASTA ======"

class Samtools_Faidx(luigi.Task):
    """
    Create the samtools fasta index for the reference genome, needed for GATK
    """
    genome = luigi.Parameter()

    def requires(self):
        return Reference_Genome_Fasta(self.genome)

    def output(self):
        return luigi.LocalTarget("fasta/" + self.genome + ".fa.fai")

    def run(self):
        run_cmd(["samtools1.3",
                 "faidx",                     # index needed for CRAM files
                 "fasta/" + self.genome + ".fa"]) # input file

        print "====== Built the samtools FASTA index ======"

class Picard_CreateSequenceDictionary(luigi.Task):
    """
    Create the sequence dictionary for the reference genome, needed for GATK
    """
    genome = luigi.Parameter()

    def requires(self):
        return Reference_Genome_Fasta(self.genome)

    def output(self):
        return luigi.LocalTarget("fasta/" + self.genome + ".fa")

    def run(self):
        run_cmd(["java", "-jar",
                 "/usr/local/picard-tools-2.5.0/picard.jar",
                 "CreateSequenceDictionary",
                 "R=fasta/" + self.genome + ".fa",        # reference fasta file
                 "O=fasta/" + self.genome + ".dict"])     # dictionary file

        print "====== Built the picard sequence dictionary ======"

class Bwa_Index_Bwtsw(luigi.Task):
    """
    Builds the BWA index for the reference genome, needed for performing alignments
    """
    genome = luigi.Parameter()

    def requires(self):
        return Reference_Genome_Fasta(self.genome)

    def output(self):
        extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
        return [luigi.LocalTarget("fasta/" + self.genome + ".fa." + ext) for ext in extensions]

    def run(self):
        run_cmd(["bwa",
                 "index",                         # index needed for bwa alignment
                 "-a", "bwtsw",                   # algorithm suitable for mammals
                 "fasta/" + self.genome + ".fa"]) # input file

        print "====== Built the BWA index ======"

class Bwa_Mem(luigi.Task):
    """
    Align the fastq files to the reference genome
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return [Sample_PairedEnd_Fastq(self.sample),
                Bwa_Index_Bwtsw(self.genome)]

    def output(self):
        return luigi.LocalTarget("sam/" + self.sample + ".sam")

    def run(self):
        global SAMPLES
        readgroup= "@RG\\tID:{id}\\tSM:{sample}".format(id=SAMPLES[self.sample],
                                                        sample=self.sample)

        # perform the alignment
        sam = run_cmd(["bwa",
                       "mem",                                # align using the mem algorithm
                       "-t", MAX_CPU_CORES,                  # number of cores
                       "-R", readgroup,                      # readgroup metadata
                       "fasta/" + self.genome + ".fa",       # reference genome
                       "fastq/" + self.sample + "_1.fastq",  # pair 1
                       "fastq/" + self.sample + "_2.fastq"]) # pair 2

        # save the SAM file
        with self.output().open('w') as fout:
            fout.write(sam)

        print "====== Aligned the FASTQ using BWA ======"

class Samtools_Sort_Bam(luigi.Task):
    """
    Convert an uncompressed SAM file into BAM format, and sort it
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Bwa_Mem(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bam/" + self.sample + ".bam")

    def run(self):
        # perform the SAM -> BAM conversion and sorting
        bam = run_cmd(["samtools1.3",
                       "sort",                         # sort the reads
                       "-l", DEFAULT_COMPRESSION,      # level of compression
                       "-@", MAX_CPU_CORES,            # number of cores
                       "-O", "bam",                    # output a BAM file
                       "sam/" + self.sample + ".sam"]) # input a SAM file

        # save the BAM file
        with self.output().open('w') as fout:
            fout.write(bam)

        print "===== Converted SAM file to BAM ======"

class Picard_MarkDuplicates(luigi.Task):
    """
    Remove PCR duplicates, so we don't overestimate coverage
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Samtools_Sort_Bam(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bam/" + self.sample + ".rmdup.bam")

    def run(self):
        run_cmd(["java", "-jar",
                 "/usr/local/picard-tools-2.5.0/picard.jar",
                 "MarkDuplicates",
                 "INPUT=bam/" + self.sample + ".bam",
                 "OUTPUT=bam/" + self.sample + ".rmdup.bam",
                 "METRICS_FILE=bam/" + self.sample + ".rmdup.txt",
                 "REMOVE_DUPLICATES=true",
                 "QUIET=true"])

        print "==== Removed Duplicates from BAM ===="

class Samtools_Index_Bam(luigi.Task):
    """
    Create an index for the BAM file, for fast random access
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("bam/" + self.sample + ".rmdup.bam.bai")

    def requires(self):
        return Picard_MarkDuplicates(self.sample, self.genome)

    def run(self):
        run_cmd(["samtools1.3",
                 "index",
                 "-b",  # create a BAI index
                 "bam/" + self.sample + ".bam"])  # file to index

        print "==== Indexing CRAM ===="

# class GATK_Variant_Call(luigi.Task):
#     """
#     Convert the BAM file into a BCF
#     """
#     sample = luigi.Parameter()
#     genome = luigi.Parameter()
#
#     def requires(self):
#         return Picard_MarkDuplicates(self.sample, self.genome)
#
#     def output(self):
#         return luigi.LocalTarget("bcf/" + self.sample + ".bcf")
#
#     def run(self):
#         bcf = run_cmd(["samtools1.3",
#                        "mpileup",
#                        "-g",                                 # call genotypes
#                        "-o", self.output().path,             # output location
#                        "-f", "fasta/" + self.genome + ".fa", # reference genome
#                        "bam/" + self.sample + ".rmdup.bam"]) # input BAM file
#
#         # TODO can't use output pipe because multithreading casues server to crash due to buffering hundreds of GB in RAM
#
#         # # save the BCF file
#         # with self.output().open('w') as fout:
#         #     fout.write(bcf)
#
#         print "===== Converted BAM file to BCF ======="

class Custom_Genome_Pipeline(luigi.Task):
    """
    Run all the samples through the pipeline
    """

    def requires(self):
        for sample in SAMPLES:
            yield Samtools_Index_Bam(sample, GENOME)

        yield Samtools_Faidx(GENOME)
        yield Picard_CreateSequenceDictionary(GENOME)

if __name__=='__main__':
    luigi.run()