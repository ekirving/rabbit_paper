#!/usr/bin/python
import luigi
import multiprocessing
import os
import os.path
import subprocess
from shutil import copyfile
import glob
import random

# import the custom vcf parser
from vcfparser import *

# Based on example pipeline from http://coderscrowd.com/app/public/codes/view/229

# the reference genome
GENOME = "OryCun2.0"
GENOME_URL = "ftp://ftp.ensembl.org/pub/release-84/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz"

# file containing the list of sequence capture regions, format <chr>:<start>-<stop>
TARGETS = './targets.interval_list'

# population, accession codes
POPULATIONS = dict() # (n=28)

# Wild mountain hare / Lepus timidus (n=1)
POPULATIONS['OUT'] = ['SRR824842']

# Domestic breeds (n=8)
POPULATIONS['DOM'] = ['SRR997325','SRR997320','SRR997321','SRR997327','SRR997323','SRR997326','SRR997324','SRR997322']

# Wild French (n=7)
POPULATIONS['WLD-FRE'] = ['SRR997319','SRR997317','SRR997304','SRR997303','SRR997318','SRR997316','SRR997305']

# Wild Iberian / Oryctolagus cuniculus algirus (n=6)
POPULATIONS['WLD-IB1'] = ['SRR827758', 'SRR827761', 'SRR827762', 'SRR827763', 'SRR827764', 'SRR827765']

# Wild Iberian / Oryctolagus cuniculus cuniculus (n=6)
POPULATIONS['WLD-IB2'] = ['SRR827759', 'SRR827760', 'SRR827766', 'SRR827767', 'SRR827768', 'SRR827769']

# the samtools flag for BAM file comression
DEFAULT_COMPRESSION = 6

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_GENOTYPE_QUAL = 30

# the maximum number of ancestral populatons to run admiture for
MAX_ANCESTRAL_K = 10

# no single worker should use more than 50% of the available cores
MAX_CPU_CORES = int(multiprocessing.cpu_count() * 0.5)

def run_cmd(cmd, shell=False):
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
                            shell=shell,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # fetch the output and error
    (stdout, stderr) = proc.communicate()

    # TODO remove when done debugging
    print vars(proc)

    # bail if something went wrong
    if proc.returncode:
        raise Exception(stderr)

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

    def output(self):
        return [luigi.LocalTarget("fastq/" + self.sample + "_1.fastq.gz"),
                luigi.LocalTarget("fastq/" + self.sample + "_2.fastq.gz")]

    def run(self):

        # NCBI SRA toolkit
        run_cmd(["fastq-dump",
                 "--gzip",              # outpute gzipped files
                 "--split-files",       # split paired end files
                 "--outdir", "./fastq", # output directory
                 self.sample])

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
        readgroup= "@RG\\tID:{sample}\\tSM:{sample}".format(sample=self.sample)

        # perform the alignment
        sam = run_cmd(["bwa",
                       "mem",                                   # align using the mem algorithm
                       "-t", MAX_CPU_CORES,                     # number of cores
                       "-R", readgroup,                         # readgroup metadata
                       "fasta/" + self.genome + ".fa",          # reference genome
                       "fastq/" + self.sample + "_1.fastq.gz",  # pair 1
                       "fastq/" + self.sample + "_2.fastq.gz"]) # pair 2

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

        # https://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation
        run_cmd(["java", "-jar",
                 "/usr/local/picard-tools-2.5.0/picard.jar",
                 "FixMateInformation",
                 "INPUT=bam/" + self.sample + ".bam"])

        # https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
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
                 "bam/" + self.sample + ".rmdup.bam"])  # file to index

        print "==== Indexing CRAM ===="

class GATK_HaplotypeCaller(luigi.Task):
    """
    Convert the BAM file into a gVCF
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        return [Picard_CreateSequenceDictionary(self.genome),
                Samtools_Faidx(self.genome),
                Samtools_Index_Bam(self.sample, self.genome)]

    def output(self):
        return luigi.LocalTarget("vcf/" + self.sample + ".vcf")

    def run(self):

        run_cmd(["java", "-Xmx8G", "-Xms8G", "-jar",
                 "../GenomeAnalysisTK.jar",
                 "-T", "HaplotypeCaller",                # use the HaplotypeCaller to call variants
                 "-R", "fasta/" + self.genome + ".fa",   # the indexed reference genome
                 "--genotyping_mode", "DISCOVERY",       # variant discovery
                 "--emitRefConfidence", "GVCF",          # reference model emitted with condensed non-variant blocks
                 "-variant_index_type", "LINEAR",
                 "-variant_index_parameter", "128000",
                 "--output_mode", "EMIT_ALL_SITES",      # produces calls at any callable site regardless of confidence
                 "-L", self.targets,                     # limit to the list of regions defined in the targets file
                 "-stand_emit_conf", "10",               # min confidence threshold
                 "-stand_call_conf", MIN_GENOTYPE_QUAL,  # min call threshold
                 "-I", "bam/" + self.sample + ".rmdup.bam",
                 "-o", "vcf/" + self.sample + ".vcf"])

        print "===== Created gVCF for sample ======="

class GATK_GenotypeGVCFs(luigi.Task):
    """
    Joint genotype all the samples in a population
    """
    population = luigi.Parameter()
    samples = luigi.ListParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        for sample in self.samples:
            yield GATK_HaplotypeCaller(sample, self.genome, self.targets)

    def output(self):
        return luigi.LocalTarget("vcf/" + self.population + ".vcf")

    def run(self):
        # make a list of input files
        vcf_files = sum([["-V", "vcf/" + sample + ".vcf"] for sample in self.samples], [])

        run_cmd(["java", "-Xmx8G", "-Xms8G", "-jar",
                 "../GenomeAnalysisTK.jar",
                 "-T", "GenotypeGVCFs",                 # use GenotypeGVCFs to jointly call variants
                 "--num_threads", MAX_CPU_CORES,        # number of data threads to allocate to this analysis
                 "--includeNonVariantSites",            # include loci found to be non-variant after genotyping
                 "-R", "fasta/" + self.genome + ".fa",  # the indexed reference genome
                 "-o", "vcf/" + self.population + ".vcf"]
                + vcf_files)

        print "===== Joint genotyped population ======="

class Site_Frequency_Spectrum(luigi.Task):
    """
    Produce the site frequency spectrum, based on genotype calls from GATK GenotypeGVCFs
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        for pop in self.populations:
            yield GATK_GenotypeGVCFs(pop, self.populations[pop], self.genome, self.targets)

    def output(self):
        return luigi.LocalTarget("fsdata/" + self.genome + ".data")

    def run(self):

        # log everything to file
        logging.basicConfig(filename="fsdata/" + self.genome + ".log", level=logging.DEBUG)

        # generate the frequency spectrum
        fsdata = generate_frequency_spectrum(self.populations)

        # save the fsdata file
        with self.output().open('w') as fout:
            fout.write(fsdata)

        print "===== Generated the site frequency spectrum  ======="

class GATK_SelectVariants(luigi.Task):
    """
    Select only biallelic variant sites for downstream analysis
    """
    population = luigi.Parameter()
    samples = luigi.ListParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        yield GATK_GenotypeGVCFs(self.population, self.samples, self.genome, self.targets)

    def output(self):
        return luigi.LocalTarget("vcf/" + self.population + ".variant.vcf")

    def run(self):

        run_cmd(["java", "-Xmx2g", "-jar",
                 "../GenomeAnalysisTK.jar",
                 "-T", "SelectVariants",
                 "-R", "fasta/" + self.genome + ".fa",
                 "-V", "vcf/" + self.population + ".vcf",
                 "-o", "vcf/" + self.population + ".variant.vcf",
                 "--selectTypeToInclude", "SNP",
                 "--restrictAllelesTo", "BIALLELIC"])

        print "===== Selected variant sites ======="


class Plink_Make_Bed(luigi.Task):
    """
    Convert a population gVCF into a BED file (binary pedigree)
    """
    population = luigi.Parameter()
    samples = luigi.ListParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        yield GATK_SelectVariants(self.population, self.samples, self.genome, self.targets)

    def output(self):
        extensions = ['bed', 'bim', 'fam']
        return [luigi.LocalTarget("bed/" + self.population + "." + ext) for ext in extensions]

    def run(self):

        # convert VCF to PED
        run_cmd(["plink",
                 "--recode",
                 "--const-fid", self.population,
                 "--biallelic-only", "strict",
                 "--vcf-min-qual", MIN_GENOTYPE_QUAL,
                 "--vcf-require-gt",
                 "--vcf", "vcf/" + self.population + ".variant.vcf",
                 "--out", "ped/" + self.population])

        # use awk to add variant IDs, so we can identify pollyallelic sites during merge
        map = run_cmd(["awk", "$2=$1\"-\"$4", "ped/" + self.population + ".map"])

        # replace the old map file
        with open("ped/" + self.population + ".map", 'w') as fout:
            fout.write(map)

        # convert PED to BED
        run_cmd(["plink",
                 "--make-bed",
                 "--file", "ped/" + self.population,
                 "--out", "bed/" + self.population])

        print "===== Created population BED file ======="


class Plink_Merge_Beds(luigi.Task):
    """
    Merge multiple BED files into one
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        for pop in self.populations:
            yield Plink_Make_Bed(pop, self.populations[pop], self.genome, self.targets)

    def output(self):
        extensions = ['bed', 'bim', 'fam']
        return [luigi.LocalTarget("bed/" + self.group + "." + ext) for ext in extensions]

    def run(self):

        # generate a unique suffix for temporary files
        suffix = 'tmp' + str(random.getrandbits(100))

        # make a copy of the bed files because we'll need to filter them
        for pop in self.populations:
            for ext in ['bed', 'bim', 'fam']:
                copyfile("bed/" + pop + "." + ext,
                         "bed/" + pop + "." + suffix + "." + ext)

        # merge requires the first bed file to be named in the command
        bed_files = ["bed/" + pop + "." + suffix for pop in self.populations]
        bed_file1 = bed_files.pop(0)

        # make the merge-list with the remaining BED files
        with open("bed/" + self.group + ".list", 'w') as fout:
            fout.write("\n".join(bed_files))

        # compose the merge command, because we are going to need it twice
        merge = ["plink",
                 "--make-bed",
                 "--bfile", bed_file1,
                 "--merge-list", "bed/" + self.group + ".list",
                 "--out", "bed/" + self.group]

        try:
            # attempt the merge
            run_cmd(merge)

        except Exception as e:

            # handle merge errors
            if os.path.isfile("bed/" + self.group + "-merge.missnp") :

                # filter all the BED files, using the missnp file created by the failed merge
                for pop in self.populations:
                    run_cmd(["plink",
                             "--make-bed",
                             "--exclude", "bed/" + self.group + "-merge.missnp",
                             "--bfile", "bed/" + pop + "." + suffix,
                             "--out", "bed/" + pop + "." + suffix])

                # reattempt the merge
                run_cmd(merge)

            else:
                raise Exception(e)

        # tidy up all the temporary intermediate files
        for tmp in glob.glob("./bed/*{}*".format(suffix)):
            os.remove(tmp)

        print "===== Merged BED files ======="


class Plink_Prune_Bed(luigi.Task):
    """
    Prune the genome BED file for linkage disequilibrium
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return Plink_Merge_Beds(self.populations, self.genome, self.targets, self.group)

    def output(self):
        extensions = ['bed', 'bim', 'fam']
        return [luigi.LocalTarget("bed/" + self.group + ".pruned." + ext) for ext in extensions]

    def run(self):

        # calcualte the prune list (prune.in / prune.out)
        run_cmd(["plink",
                 "--indep-pairwise", 50, 10, 0.1,
                 "--bfile", "bed/" + self.group,
                 "--out", "bed/" + self.group])

        # apply the prune list
        run_cmd(["plink",
                 "--make-bed",
                 "--extract", "bed/" + self.group + ".prune.in",
                 "--bfile", "bed/" + self.group,
                 "--out", "bed/" + self.group + ".pruned"])

        print "===== Pruned BED file ======="


class Admixture_K(luigi.Task):
    """
    Run admixture, with K ancentral populations, on the pruned BED files
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()
    k = luigi.Parameter()

    def requires(self):
        return Plink_Prune_Bed(self.populations, self.genome, self.targets, self.group)

    def output(self):
        extensions = ['P', 'Q', 'log']
        return [luigi.LocalTarget("admix/" + self.group + ".pruned." + str(self.k) + "." + ext) for ext in extensions]

    def run(self):

        # admixture only outputs to the current directory
        os.chdir('./admix')

        log = run_cmd(["admixture",
                       "-j{}".format(MAX_CPU_CORES),           # use multi-threading
                       "--cv",                                 # include cross-validation standard errors
                       "../bed/" + self.group + ".pruned.bed", # using this input file
                       self.k])                                # for K ancestral populations

        # restore previous working directory
        os.chdir('..')

        # save the log file
        with self.output()[2].open('w') as fout:
            fout.write(log)

        # TODO plot the admixture stacked column chart

        print "===== Admixture ======="

class Flashpca(luigi.Task):
    """
    Run flashpca on the pruned BED files
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return Plink_Prune_Bed(self.populations, self.genome, self.targets, self.group)

    def output(self):
        prefixes = ['eigenvalues', 'eigenvectors', 'pcs', 'pve']
        return [luigi.LocalTarget("flashpca/" + prefix + "_" + self.group + ".txt") for prefix in prefixes]

    def run(self):

        # flashpca only outputs to the current directory
        os.chdir('./flashpca')

        log = run_cmd(["flashpca",
                       "--bfile", "../bed/" + self.group + ".pruned",
                       "--numthreads", MAX_CPU_CORES,
                       "--suffix", "_" + self.group + ".txt"])

        # restore previous working directory
        os.chdir('..')

        print "===== Flashpca ======="

class RScript_Ggplot_PCA(luigi.Task):
    """
    Use ggplot to plot the PCA
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return Flashpca(self.populations, self.genome, self.targets, self.group)

    def output(self):
        extensions = ['data', 'pdf']
        return [luigi.LocalTarget("flashpca/pca_" + self.group + "." + ext) for ext in extensions]

    def run(self):

        # use awk to add population and sample names, needed for the PCA plot
        data = run_cmd(["awk '{print $1\"\t\"$2}' bed/" + self.group + ".fam | paste - flashpca/pcs_" + self.group + ".txt"], True)

        # save the labeled file
        with self.output()[0].open('w') as fout:
            fout.write(data)

        # generate a PNG of the PCA plot
        run_cmd(["Rscript",
                 "flashpca.R",
                 "flashpca/pca_" + self.group])

        print "===== RScript ggplot ======="

class RScript_Tree(luigi.Task):
    """
    Create a phylogenetic tree from a pruned BED file
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return Plink_Prune_Bed(self.populations, self.genome, self.targets, self.group)

    def output(self):
        return luigi.LocalTarget("tree/" + self.group + ".tree")

    def run(self):

        # make the distance matrix
        run_cmd(["plink",
                 "--distance", "square", "1-ibs",
                 "--bfile", "bed/" + self.group + ".pruned",
                 "--out", "tree/" + self.group])

        # use awk to extract short versions of the pop and sample names
        awk = "awk '{ print substr($1, length($1)-2, 3) \"-\" substr($2, length($2)-2, 3) }' ./bed/" + self.group + ".fam"

        # fetch the names as a row
        header = run_cmd([awk + " | xargs"], True)

        # add the samples names as a column to the mdist data
        data = run_cmd([awk + " | paste - ./bed/" + self.group + ".mdist"], True)

        # save the labeled file
        with open("tree/" + self.group + ".data", 'w') as fout:
            fout.write("\t"+header)
            fout.write(data)

        # generate a tree from the labeled data
        run_cmd(["Rscript",
                 "tree.R",
                 "tree/" + self.group + ".data",
                 "tree/" + self.group + ".tree"])

        # TODO plot the tree with either ggtree or ape, reusing the colour scheme from the PCA

        print "===== RScript tree ======="

class Custom_Genome_Pipeline(luigi.Task):
    """
    Run all the samples through the pipeline
    """

    def requires(self):

        # make the SFS for dadi
        yield Site_Frequency_Spectrum(POPULATIONS, GENOME, TARGETS)

        # setup some population groups
        groups = dict()

        # one with all the populations
        groups['all-pops'] = POPULATIONS.copy()

        # and one without the outgroup
        groups['no-outgroup'] = POPULATIONS.copy()
        del groups['no-outgroup']['OUT']

        # and one without either the outgroup or the most divergent wild species (O. cuniculus algirus)
        groups['no-out-ib1'] = POPULATIONS.copy()
        del groups['no-out-ib1']['OUT']
        del groups['no-out-ib1']['WLD-IB1']

        for group in groups:
            # run admixture or each population and each value of K
            for k in range(1, MAX_ANCESTRAL_K + 1):
                yield Admixture_K(groups[group], GENOME, TARGETS, group, k)

            # run flashpca for each population
            yield RScript_Ggplot_PCA(groups[group], GENOME, TARGETS, group)

        yield RScript_Tree(POPULATIONS, GENOME, TARGETS, 'all-pops')

if __name__=='__main__':
    luigi.run()