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


class PairedEnd_Fastq(luigi.Task):
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

class Genome_Fasta(luigi.Task):
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

class Samtools_Fasta_Index(luigi.Task):
    """
    Builds the samtools fasta index for the reference genome, needed for downstream CRAM files
    """
    genome = luigi.Parameter()

    def requires(self):
        return Genome_Fasta(self.genome)

    def output(self):
        return luigi.LocalTarget("fasta/" + self.genome + ".fa")

    def run(self):
        run_cmd(["samtools1.3",
                 "faidx",                     # index needed for CRAM files
                 "fasta/" + self.genome + ".fa"]) # input file

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
        return [luigi.LocalTarget("fasta/" + self.genome + ".fa." + ext) for ext in extensions]

    def run(self):
        run_cmd(["bwa",
                 "index",                     # index needed for bwa alignment
                 "-a", "bwtsw",               # algorithm suitable for mammals
                 "fasta/" + self.genome + ".fa"]) # input file

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

class Convert_Sam_Sorted_Bam(luigi.Task):
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

class Index_Bam(luigi.Task):
    """
    Create an index for the BAM file, for fast random access
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("bam/" + self.sample + ".bam.bai")

    def requires(self):
        return Convert_Sam_Sorted_Bam(self.sample, self.genome)

    def run(self):
        run_cmd(["samtools1.3",
                 "index",
                 "-b",                           # create a BAI index
                 "bam/" + self.sample + ".bam"]) # file to index

        print "==== Indexing CRAM ===="

class Remove_Duplicates_BAM(luigi.Task):
    """
    Remove PCR duplicates, so we don't overestimate coverage
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Index_Bam(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bam/" + self.sample + ".rmdup.bam")

    def run(self):
        run_cmd(["java", "-jar", "/usr/local/picard-tools-2.5.0/picard.jar",
                 "MarkDuplicates"
                 "INPUT=" + self.sample + ".bam",
                 "OUTPUT=" + self.sample + ".rmdup.bam",
                 "METRICS_FILE=" + self.sample + ".rmdup.txt",
                 "REMOVE_DUPLICATES=true",
                 "QUIET=true"])

        print "==== Removed Duplicates from BAM ===="

# class Convert_Bam_Cram(luigi.Task):
#     """
#     Convert a BAM file into the smaller CRAM format
#     """
#     sample = luigi.Parameter()
#     genome = luigi.Parameter()
#
#     def requires(self):
#         return [Index_Bam(self.sample, self.genome),
#                 Samtools_Fasta_Index(self.genome)]
#
#     def output(self):
#         return luigi.LocalTarget("cram/" + self.sample + ".cram")
#
#     def run(self):
#         # perform the BAM -> CRAM conversion
#         cram = run_cmd(["samtools1.3",
#                         "view",
#                         "-@", MAX_CPU_CORES,                  # number of cores
#                         "-C",                                 # output a CRAM file
#                         "-T", "fasta/" + self.genome + ".fa", # reference genome
#                         "bam/" + self.sample + ".bam"])       # input BAM file
#
#         # save the CRAM file
#         with self.output().open('w') as fout:
#             fout.write(cram)
#
#         print "===== Converted BAM file to CRAM ======"
#
# class Index_Cram(luigi.Task):
#     """
#     Create an index for the CRAM file, for fast random access
#     """
#     sample = luigi.Parameter()
#     genome = luigi.Parameter()
#
#     def output(self):
#         return luigi.LocalTarget("cram/" + self.sample + ".cram.crai")
#
#     def requires(self):
#         return Convert_Bam_Cram(self.sample, self.genome)
#
#     def run(self):
#         run_cmd(["samtools1.3",
#                  "index",
#                  "cram/" + self.sample + ".cram"])
#
#         print "==== Indexing CRAM ===="
#
# class Samtools_MPileup(luigi.Task):
#     """
#     Convert the BAM file into a pileup
#     """
#     sample = luigi.Parameter()
#     genome = luigi.Parameter()
#
#     def requires(self):
#         return Remove_Duplicates_BAM(self.sample, self.genome)
#
#     def output(self):
#         return luigi.LocalTarget("pileup/" + self.sample + ".pileup")
#
#     def run(self):
#         pileup = run_cmd(["samtools1.3",
#                           "mpileup",
#                           "-o",  self.output().path,            # output location
#                           "-f", "fasta/" + self.genome + ".fa", # reference genome
#                           "bam/" + self.sample + ".rmdup.bam"]) # input BAM file
#
#         # TODO can't use output pipe because multithreading casues server to crash due to buffering hundreds of GB in RAM
#
#         # save the pileup file
#         # with self.output().open('w') as fout:
#         #     fout.write(pileup)
#
#         print "===== Converted CRAM file to pileup ======="
#
# class Coverage_Stats(luigi.Task):
#     """
#     Parse the pileup and compute depth of coverage statistics
#
#     NB. Method is very memory hungry
#     """
#     sample = luigi.Parameter()
#     genome = luigi.Parameter()
#
#     def requires(self):
#         return Samtools_MPileup(self.sample, self.genome)
#
#     def output(self):
#         return luigi.LocalTarget("pileup/" + self.sample + ".pileup.stats")
#
#     def run(self):
#
#         # parse the pileup and load the coverage data into memory (expensive!!)
#         coverage = numpy.loadtxt("pileup/" + self.sample + ".pileup", dtype=int, delimiter=' ', usecols=[3])
#
#         # make a dictionary
#         stats = {}
#
#         # calcluate the 1st and 3rd quartiles
#         (stats['q1'], stats['q2']) = numpy.percentile(coverage, [25, 75])
#
#         # seralize the dictionary for later use
#         with self.output().open('w') as fout:
#             fout.write(json.dumps(stats))
#
#         print "===== Calcualted depth of coverage ======="


class Convert_Bam_Bcf(luigi.Task):
    """
    Convert the BAM file into a BCF
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Remove_Duplicates_BAM(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bcf/" + self.sample + ".bcf")

    def run(self):
        bcf = run_cmd(["samtools1.3",
                       "mpileup",
                       "-g",                                 # call genotypes
                       "-o", self.output().path,             # output location
                       "-f", "fasta/" + self.genome + ".fa", # reference genome
                       "bam/" + self.sample + ".rmdup.bam"]) # input BAM file

        # TODO can't use output pipe because multithreading casues server to crash due to buffering hundreds of GB in RAM

        # # save the BCF file
        # with self.output().open('w') as fout:
        #     fout.write(bcf)

        print "===== Converted BAM file to BCF ======="

class Convert_Bcf_Vcf(luigi.Task):
    """
    Convert the BCF file into a VCF
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Convert_Bam_Bcf(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("vcf/" + self.sample + ".vcf")

    def run(self):
        vcf = run_cmd(["bcftools",
                        "call",
                        "-c",                           # consensus caller
                        "-o", self.output().path,       # output location
                        "-O", "v",                      # output uncompressed VCF
                        "bcf/" + self.sample + ".bcf"]) # input BCF file

        # TODO can't use output pipe because multithreading casues server to crash due to buffering hundreds of GB in RAM

        # # save the VCF file
        # with self.output().open('w') as fout:
        #     fout.write(vcf)

        print "===== Called variant sites ======="

class Site_Frequency_Spectrum(luigi.Task):
    """
    Produce the site frequency spectrum, based on genotype calls from bcftools
    """
    samples = luigi.ListParameter()
    genome = luigi.Parameter()

    def requires(self):
        for sample in self.samples:
            yield Bcftools_Call(sample, self.genome)

    def output(self):
        return luigi.LocalTarget("fsdata/" + self.genome + ".data")

    def run(self):

        # log everything to file
        logging.basicConfig(filename="fsdata/" + self.genome + ".log", level=logging.DEBUG)

        # generate the frequency spectrum
        fsdata = generate_frequency_spectrum(self.samples, WILD_THRESHOLD)

        # save the fsdata file
        with self.output().open('w') as fout:
            fout.write(fsdata)

        print "===== Generated the site frequency spectrum  ======="

class Custom_Genome_Pipeline(luigi.Task):
    """
    Run all the samples through the pipeline
    """

    def requires(self):
        return Site_Frequency_Spectrum(SAMPLES, GENOME)

if __name__=='__main__':
    luigi.run()