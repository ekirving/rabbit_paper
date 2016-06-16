from subprocess import check_output
import subprocess
import luigi
import os
import random
import multiprocessing
import numpy
import json

# Based on example pipleine from http://coderscrowd.com/app/public/codes/view/229

# the reference genome
GENOME = "OryCun2.0"
GENOME_URL = "ftp://ftp.ensembl.org/pub/release-84/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz"

# the list of sample codes
# SAMPLES = ['SRR997303','SRR997304','SRR997305','SRR997316','SRR997317','SRR997318','SRR997319','SRR997320',
#            'SRR997321','SRR997322','SRR997323','SRR997324','SRR997325','SRR997326','SRR997327']

SAMPLES = ['SRR997317','SRR997321', 'SRR997322']

# the URLs for the paried end fastq files
pair1_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_1.fastq.gz"
pair2_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_2.fastq.gz"

# the samtools flag for BAM file comression
DEFAULT_COMPRESSION = 6

# TODO find out how many workers there are
# no single worker should use more than 50% of the available cores
MAX_CPU_CORES = int(multiprocessing.cpu_count() * 0.5)

# consume all the available cores
# MAX_CPU_CORES = multiprocessing.cpu_count()

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

def Unzip_File(gzip):
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
        global pair1_url, pair2_url
        # download the gzipped fastq files
        return [Curl_Download(pair1_url.format(sample=self.sample), "fastq/"+self.sample+"_1.fastq.gz"),
                Curl_Download(pair2_url.format(sample=self.sample), "fastq/"+self.sample+"_2.fastq.gz")]

    def output(self):
        return [luigi.LocalTarget("fastq/"+self.sample+"_1.fastq"),
                luigi.LocalTarget("fastq/"+self.sample+"_2.fastq")]

    def run(self):
        # fetch the downloaded files
        (pair1, pair2) = self.requires()

        # unzip them
        fastq1 = Unzip_File(pair1.file)
        fastq2 = Unzip_File(pair2.file)

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
        return Curl_Download(GENOME_URL, "fasta/"+self.genome+".fa.gz")

    def output(self):
        return luigi.LocalTarget("fasta/"+self.genome+".fa")

    def run(self):
        # unzip the file
        fasta = Unzip_File(self.requires().file)

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
        return luigi.LocalTarget("fasta/"+self.genome+".fa")

    def run(self):
        run_cmd(["samtools1.3",
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
                       "-t", MAX_CPU_CORES,              # number of cores
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

    def requires(self):
        return Bwa_Mem(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bam/"+self.sample+".bam")

    def run(self):
        # perform the SAM -> BAM conversion and sorting
        bam = run_cmd(["samtools1.3",
                       "sort",                     # sort the reads
                       "-l", DEFAULT_COMPRESSION,  # level of compression
                       "-@", MAX_CPU_CORES,        # number of cores
                       "-O", "bam",                # output a BAM file
                       "sam/"+self.sample+".sam"]) # input a SAM file

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
                 "-b",                       # create a BAI index
                 "bam/"+self.sample+".bam"]) # file to index

        print "==== Indexing CRAM ===="

class Convert_Bam_Cram(luigi.Task):
    """
    Convert a BAM file into the smaller CRAM format
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return [Index_Bam(self.sample, self.genome),
                Samtools_Fasta_Index(self.genome)]

    def output(self):
        return luigi.LocalTarget("cram/"+self.sample+".cram")

    def run(self):
        # perform the BAM -> CRAM conversion
        cram = run_cmd(["samtools1.3",
                        "view",
                        "-@", MAX_CPU_CORES,              # number of cores
                        "-C",                             # output a CRAM file
                        "-T", "fasta/"+self.genome+".fa", # reference genome
                        "bam/"+self.sample+".bam"])       # input BAM file

        # save the CRAM file
        with self.output().open('w') as fout:
            fout.write(cram)

        print "===== Converted BAM file to CRAM ======"

class Index_Cram(luigi.Task):
    """
    Create an index for the CRAM file, for fast random access
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("cram/" + self.sample + ".cram.crai")

    def requires(self):
        return Convert_Bam_Cram(self.sample, self.genome)

    def run(self):
        run_cmd(["samtools1.3",
                 "index",
                 "cram/" + self.sample + ".cram"])

        print "==== Indexing CRAM ===="

class Samtools_MPileup(luigi.Task):
    """
    Convert the CRAM file into a pileup
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Index_Cram(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("pileup/"+self.sample+".pileup")

    def run(self):
        pileup = run_cmd(["samtools1.3",
                          "mpileup",                             # output a pileup file
                          "-o", "pileup/"+self.sample+".pileup", # output location
                          "-f", "fasta/"+self.genome+".fa",      # reference genome
                          "cram/"+self.sample+".cram"])          # input CRAM file

        # TODO can't use output pipe because it fails with error "IOError: [Errno 22] Invalid argument"

        # save the pileup file
        # with self.output().open('w') as fout:
        #     fout.write(pileup)

        print "===== Converted CRAM file to pileup ======="

class Coverage_Stats(luigi.Task):
    """
    Parse the pileup and compute depth of coverage statistics

    NB. Method is very memory hungry
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Samtools_MPileup(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("pileup/"+self.sample+".pileup.stats")

    def run(self):

        # parse the pileup and load the coverage data into memory (expensive!!)
        coverage = numpy.loadtxt("pileup/" + self.sample + ".pileup", dtype=int, delimiter=' ', usecols=[3])

        # make a dictionary
        stats = {}

        # calcluate the 1st and 3rd quartiles
        (stats['q1'], stats['q2']) = numpy.percentile(coverage, [25, 75])

        # seralize the dictionary for later use
        with self.output().open('w') as fout:
            fout.write(json.dumps(stats))

        print "===== Calcualted depth of coverage ======="


class Convert_Cram_Bcf(luigi.Task):
    """
    Convert the BAM file into a BCF
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Index_Cram(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bcf/"+self.sample+".bcf")

    def run(self):
        pileup = run_cmd(["samtools1.3",
                          "mpileup",
                          "-g",                              # call genotypes
                          "-o", "bcf/"+self.sample+".bcf",   # output location
                          "-f", "fasta/"+self.genome+".fa",  # reference genome
                          "cram/"+self.sample+".cram"])      # input CRAM file

        # TODO can't use output pipe because it fails with error "IOError: [Errno 22] Invalid argument"

        # save the pileup file
        # with self.output().open('w') as fout:
        #     fout.write(pileup)

        print "===== Converted CRAM file to pileup ======="

class Bcftools_Call(luigi.Task):
    """
    Convert the BAM file into a BCF
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Convert_Cram_Bcf(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("vcf/"+self.sample+".vcf")

    def run(self):
        vcf = run_cmd(["bcftools",
                        "call",
                        "-v",                              # output variant sites only
                        "-m",                              # multiallelic-caller
                        "-O", "v",                         # output uncompressed VCF
                        # "-o", "vcf/"+self.sample+".vcf", # output location
                        "bcf/"+self.sample+".bcf"])        # input BCF file

        # save the VCF file
        with self.output().open('w') as fout:
            fout.write(vcf)

        print "===== Called variant sites ======="

class Custom_Genome_Pipeline(luigi.Task):
    """
    Run all the samples through the pipeline
    """
    def requires(self):
        for sample in SAMPLES:
            yield Bcftools_Call(sample, GENOME)

if __name__=='__main__':
    luigi.run()