#!/usr/bin/python
import luigi
import multiprocessing
import os
import os.path
import subprocess
from shutil import copyfile
import glob
import random
import re

# import the custom vcf parser
from vcfparser import *

# Based on example pipeline from http://coderscrowd.com/app/public/codes/view/229

# the reference genome
GENOME = "OryCun2.0"
GENOME_URL = "ftp://ftp.ensembl.org/pub/release-84/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0" \
             ".dna.toplevel.fa.gz"

# file containing the list of sequence capture regions, format <chr>:<start>-<stop>
TARGETS = './targets.interval_list'

# populations and sample accession codes (n=28)
POPULATIONS = {

    # Wild mountain hare / Lepus timidus (n=1)
    'OUT': ['SRR824842'],

    # Domestic breeds (n=8)
    'DOM': ['SRR997325', 'SRR997320', 'SRR997321', 'SRR997327', 'SRR997323', 'SRR997326', 'SRR997324', 'SRR997322'],

    # Wild French (n=7)
    'WLD-FRE': ['SRR997319', 'SRR997317', 'SRR997304', 'SRR997303', 'SRR997318', 'SRR997316', 'SRR997305'],

    # Wild Iberian / Oryctolagus cuniculus algirus (n=6)
    'WLD-IB1': ['SRR827758', 'SRR827761', 'SRR827762', 'SRR827763', 'SRR827764', 'SRR827765'],

    # Wild Iberian / Oryctolagus cuniculus cuniculus (n=6)
    'WLD-IB2': ['SRR827759', 'SRR827760', 'SRR827766', 'SRR827767', 'SRR827768', 'SRR827769']
}

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
    :param shell: Use the native shell
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
                        "-c",                 # output to stdout
                        "-p", MAX_CPU_CORES,  # use the maximum cores
                        gzip])                # gzip file

    except OSError as e:
        # if it's not installed
        if e.errno == os.errno.ENOENT:

            # use single-threaded gunzip
            return run_cmd(["gunzip",
                            "-c",             # output to stdout
                            gzip])            # gzip file
        else:
            # escalate the exception
            raise e


def curl_download(url, filename):
    """
    Downloads a remote url to a local file path using cURL
    """

    # download the file
    run_cmd(["curl",
             "-s",                  # download silently
             "--output", filename,  # output path
             url])                  # from this url


class SamplePairedEndFastq(luigi.Task):
    """
    Fetches the two paired-end FASTQ files for a given sample ID
    """
    sample = luigi.Parameter()

    def output(self):
        return [luigi.LocalTarget("fastq/{0}_{1}.fastq.gz".format(self.sample, i)) for i in [1, 2]]

    def run(self):

        # use the NCBI SRA toolkit to fetch the fastq files
        run_cmd(["fastq-dump",
                 "--gzip",               # output gzipped files
                 "--split-files",        # split paired end files
                 "--outdir", "./fastq",  # output directory
                 self.sample])


class ReferenceGenomeFasta(luigi.Task):
    """
    Fetches the reference genome
    """
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("fasta/{0}.fa".format(self.genome))

    def run(self):

        zipfile = "fasta/{0}.fa.gz".format(self.genome)

        # TODO install ensembl perl API for a generic method of downloading fasta files
        # download the gzipped fasta file of the reference genome
        curl_download(GENOME_URL, zipfile)

        # unzip the file
        fasta = unzip_file(zipfile)

        with self.output().open('w') as fout:
            fout.write(fasta)


class SamtoolsFaidx(luigi.Task):
    """
    Create the samtools fasta index for the reference genome, needed for GATK
    """
    genome = luigi.Parameter()

    def requires(self):
        return ReferenceGenomeFasta(self.genome)

    def output(self):
        return luigi.LocalTarget("fasta/{0}.fa.fai".format(self.genome))

    def run(self):
        run_cmd(["samtools1.3",
                 "faidx",                              # index needed for CRAM files
                 "fasta/{0}.fa".format(self.genome)])  # input file


class PicardCreateSequenceDictionary(luigi.Task):
    """
    Create the sequence dictionary for the reference genome, needed for GATK
    """
    genome = luigi.Parameter()

    def requires(self):
        return ReferenceGenomeFasta(self.genome)

    def output(self):
        return luigi.LocalTarget("fasta/{0}.fa".format(self.genome))

    def run(self):
        run_cmd(["java", "-jar", "$PICARD",
                 "CreateSequenceDictionary",
                 "R=fasta/{0}.fa".format(self.genome),     # reference fasta file
                 "O=fasta/{0}.dict".format(self.genome)])  # dictionary file


class BwaIndexBwtsw(luigi.Task):
    """
    Builds the BWA index for the reference genome, needed for performing alignments
    """
    genome = luigi.Parameter()

    def requires(self):
        return ReferenceGenomeFasta(self.genome)

    def output(self):
        extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
        return [luigi.LocalTarget("fasta/{0}.fa.{1}".format(self.genome, ext)) for ext in extensions]

    def run(self):

        run_cmd(["bwa",
                 "index",                              # index needed for bwa alignment
                 "-a", "bwtsw",                        # algorithm suitable for mammals
                 "fasta/{0}.fa".format(self.genome)])  # input file


class BwaMem(luigi.Task):
    """
    Align the fastq files to the reference genome
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return [SamplePairedEndFastq(self.sample),
                BwaIndexBwtsw(self.genome)]

    def output(self):
        return luigi.LocalTarget("sam/{0}.sam".format(self.sample))

    def run(self):
        read_group = "@RG\\tID:{sample}\\tSM:{sample}".format(sample=self.sample)

        # perform the alignment
        sam = run_cmd(["bwa",
                       "mem",                                        # align using the mem algorithm
                       "-t", MAX_CPU_CORES,                          # number of cores
                       "-R", read_group,                             # read group metadata
                       "fasta/{0}.fa".format(self.genome),           # reference genome
                       "fastq/{0}_1.fastq.gz".format(self.sample),   # pair 1
                       "fastq/{0}_2.fastq.gz".format(self.sample)])  # pair 2

        # save the SAM file
        with self.output().open('w') as fout:
            fout.write(sam)


class SamtoolsSortBam(luigi.Task):
    """
    Convert an uncompressed SAM file into BAM format, and sort it
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return BwaMem(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bam/{0}.bam".format(self.sample))

    def run(self):
        # perform the SAM -> BAM conversion and sorting
        bam = run_cmd(["samtools1.3",
                       "sort",                              # sort the reads
                       "-l", DEFAULT_COMPRESSION,           # level of compression
                       "-@", MAX_CPU_CORES,                 # number of cores
                       "-O", "bam",                         # output a BAM file
                       "sam/{0}.sam".format(self.sample)])  # input a SAM file

        # save the BAM file
        with self.output().open('w') as fout:
            fout.write(bam)


class PicardMarkDuplicates(luigi.Task):
    """
    Remove PCR duplicates, so we don't overestimate coverage
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return SamtoolsSortBam(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bam/{0}.rmdup.bam".format(self.sample))

    def run(self):

        # https://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation
        run_cmd(["java", "-jar", "$PICARD",
                 "FixMateInformation",
                 "INPUT=bam/{0}.bam".format(self.sample)])

        # https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
        run_cmd(["java", "-jar", "$PICARD",
                 "MarkDuplicates",
                 "INPUT=bam/{0}.bam".format(self.sample),
                 "OUTPUT=bam/{0}.rmdup.bam".format(self.sample),
                 "METRICS_FILE=bam/{0}.rmdup.txt".format(self.sample),
                 "REMOVE_DUPLICATES=true",
                 "QUIET=true"])


class SamtoolsIndexBam(luigi.Task):
    """
    Create an index for the BAM file, for fast random access
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("bam/{0}.rmdup.bam.bai".format(self.sample))

    def requires(self):
        return PicardMarkDuplicates(self.sample, self.genome)

    def run(self):

        run_cmd(["samtools1.3",
                 "index",
                 "-b",                                      # create a BAI index
                 "bam/{0}.rmdup.bam".format(self.sample)])  # file to index


class GatkHaplotypeCaller(luigi.Task):
    """
    Convert the BAM file into a gVCF
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        return [PicardCreateSequenceDictionary(self.genome),
                SamtoolsFaidx(self.genome),
                SamtoolsIndexBam(self.sample, self.genome)]

    def output(self):
        return luigi.LocalTarget("vcf/{0}.vcf".format(self.sample))

    def run(self):

        run_cmd(["java", "-Xmx8G", "-jar", "$GATK",
                 "-T", "HaplotypeCaller",                   # use the HaplotypeCaller to call variants
                 "-R", "fasta/{0}.fa".format(self.genome),  # the indexed reference genome
                 "--genotyping_mode", "DISCOVERY",          # variant discovery
                 "--emitRefConfidence", "GVCF",             # reference model emitted with condensed non-variant blocks
                 "-variant_index_type", "LINEAR",
                 "-variant_index_parameter", "128000",
                 "--output_mode", "EMIT_ALL_SITES",         # produces calls at any callable site
                 "-L", self.targets,                        # limit to the list of regions defined in the targets file
                 "-stand_emit_conf", "10",                  # min confidence threshold
                 "-stand_call_conf", MIN_GENOTYPE_QUAL,     # min call threshold
                 "-I", "bam/{0}.rmdup.bam".format(self.sample),
                 "-o", "vcf/{0}.vcf".format(self.sample)])


class GatkGenotypeGVCFs(luigi.Task):
    """
    Joint genotype all the samples in a population
    """
    population = luigi.Parameter()
    samples = luigi.ListParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        for sample in self.samples:
            yield GatkHaplotypeCaller(sample, self.genome, self.targets)

    def output(self):
        return luigi.LocalTarget("vcf/{0}.vcf".format(self.population))

    def run(self):
        # make a list of input files
        vcf_files = sum([["-V", "vcf/{0}.vcf".format(sample)] for sample in self.samples], [])

        run_cmd(["java", "-Xmx8G", "-jar", "$GATK",
                 "-T", "GenotypeGVCFs",                     # use GenotypeGVCFs to jointly call variants
                 "--num_threads", MAX_CPU_CORES,            # number of data threads to allocate to this analysis
                 "--includeNonVariantSites",                # include loci found to be non-variant after genotyping
                 "-R", "fasta/{0}.fa".format(self.genome),  # the indexed reference genome
                 "-o", "vcf/{0}.vcf".format(self.population)]
                + vcf_files)


class SiteFrequencySpectrum(luigi.Task):
    """
    Produce the site frequency spectrum, based on genotype calls from GATK GenotypeGVCFs
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        for pop in self.populations:
            yield GatkGenotypeGVCFs(pop, self.populations[pop], self.genome, self.targets)

    def output(self):
        return luigi.LocalTarget("fsdata/{0}.data".format(self.genome))

    def run(self):

        # log everything to file
        logging.basicConfig(filename="fsdata/{0}.log".format(self.genome), level=logging.DEBUG)

        # generate the frequency spectrum
        fsdata = generate_frequency_spectrum(self.populations)

        # save the fsdata file
        with self.output().open('w') as fout:
            fout.write(fsdata)


class GatkSelectVariants(luigi.Task):
    """
    Select only biallelic variant sites for downstream analysis
    """
    population = luigi.Parameter()
    samples = luigi.ListParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        yield GatkGenotypeGVCFs(self.population, self.samples, self.genome, self.targets)

    def output(self):
        return luigi.LocalTarget("vcf/{0}.variant.vcf".format(self.population))

    def run(self):

        run_cmd(["java", "-Xmx2g", "-jar", "$GATK",
                 "-T", "SelectVariants",
                 "-R", "fasta/{0}.fa".format(self.genome),
                 "-V", "vcf/{0}.vcf".format(self.population),
                 "-o", "vcf/{0}.variant.vcf".format(self.population),
                 "--selectTypeToInclude", "SNP",
                 "--restrictAllelesTo", "BIALLELIC"])


class PlinkMakeBed(luigi.Task):
    """
    Convert a population gVCF into a BED file (binary pedigree)
    """
    population = luigi.Parameter()
    samples = luigi.ListParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()

    def requires(self):
        yield GatkSelectVariants(self.population, self.samples, self.genome, self.targets)

    def output(self):
        extensions = ['bed', 'bim', 'fam']
        return [luigi.LocalTarget("bed/{0}.{1}".format(self.population, ext)) for ext in extensions]

    def run(self):

        # convert VCF to PED
        run_cmd(["plink",
                 "--recode",
                 "--const-fid", self.population,
                 "--biallelic-only", "strict",
                 "--vcf-min-qual", MIN_GENOTYPE_QUAL,
                 "--vcf-require-gt",
                 "--vcf", "vcf/{0}.variant.vcf".format(self.population),
                 "--out", "ped/{0}".format(self.population)])

        # use awk to add variant IDs, so we can identify polyallelic sites during merge
        map = run_cmd(["awk '$2=$1\"-\"$4' ped/" + str(self.population) + ".map", True])

        # replace the old map file
        with open("ped/{0}.map".format(self.population), 'w') as fout:
            fout.write(map)

        # convert PED to BED
        run_cmd(["plink",
                 "--make-bed",
                 "--file", "ped/{0}".format(self.population),
                 "--out", "bed/{0}".format(self.population)])


class PlinkMergeBeds(luigi.Task):
    """
    Merge multiple BED files into one
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        for pop in self.populations:
            yield PlinkMakeBed(pop, self.populations[pop], self.genome, self.targets)

    def output(self):
        extensions = ['bed', 'bim', 'fam']
        return [luigi.LocalTarget("bed/{0}.{1}".format(self.group, ext)) for ext in extensions]

    def run(self):

        # TODO should this use the group name instead?
        # generate a unique suffix for temporary files
        suffix = 'tmp' + str(random.getrandbits(100))

        # make a copy of the bed files because we'll need to filter them
        for pop in self.populations:
            for ext in ['bed', 'bim', 'fam']:
                copyfile("bed/{0}.{1}".format(pop, ext),
                         "bed/{0}.{1}.{2}".format(pop, suffix, ext))

        # merge requires the first bed file to be named in the command
        bed_files = ["bed/{0}.{1}".format(pop, suffix) for pop in self.populations]
        bed_file1 = bed_files.pop(0)

        # make the merge-list with the remaining BED files
        with open("bed/{0}.list".format(self.group), 'w') as fout:
            fout.write("\n".join(bed_files))

        # compose the merge command, because we are going to need it twice
        merge = ["plink",
                 "--make-bed",
                 "--bfile", bed_file1,
                 "--merge-list", "bed/{0}.list".format(self.group),
                 "--out", "bed/{0}".format(self.group)]

        try:
            # attempt the merge
            run_cmd(merge)

        except Exception as e:

            # handle merge errors
            if os.path.isfile("bed/{0}-merge.missnp".format(self.group)) :

                # filter all the BED files, using the missnp file created by the failed merge
                for pop in self.populations:
                    run_cmd(["plink",
                             "--make-bed",
                             "--exclude", "bed/{0}-merge.missnp".format(self.group),
                             "--bfile", "bed/{0}.{1}".format(pop, suffix),
                             "--out", "bed/{0}.{1}".format(pop, suffix)])

                # reattempt the merge
                run_cmd(merge)

            else:
                raise Exception(e)

        # tidy up all the temporary intermediate files
        for tmp in glob.glob("./bed/*{}*".format(suffix)):
            os.remove(tmp)


class PlinkPruneBed(luigi.Task):
    """
    Prune the genome BED file for linkage disequilibrium
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return PlinkMergeBeds(self.populations, self.genome, self.targets, self.group)

    def output(self):
        extensions = ['bed', 'bim', 'fam']
        return [luigi.LocalTarget("bed/{0}.pruned.{1}".format(self.group, ext)) for ext in extensions]

    def run(self):

        # calcualte the prune list (prune.in / prune.out)
        run_cmd(["plink",
                 "--indep-pairwise", 50, 10, 0.1,
                 "--bfile", "bed/{0}".format(self.group),
                 "--out", "bed/{0}".format(self.group)])

        # apply the prune list
        run_cmd(["plink",
                 "--make-bed",
                 "--extract", "bed/{0}.prune.in".format(self.group),
                 "--bfile", "bed/{0}".format(self.group),
                 "--out", "bed/{0}.pruned".format(self.group)])


class AdmixtureK(luigi.Task):
    """
    Run admixture, with K ancestral populations, on the pruned BED files
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()
    k = luigi.IntParameter()

    def requires(self):
        return PlinkPruneBed(self.populations, self.genome, self.targets, self.group)

    def output(self):
        extensions = ['P', 'Q', 'log']
        return [luigi.LocalTarget("admix/{0}.pruned.{1}.{2}".format(self.group, self.k, ext)) for ext in extensions]

    def run(self):

        # admixture only outputs to the current directory
        os.chdir('./admix')

        # TODO what about bootstrapping (-B / -B2000)
        log = run_cmd(["admixture",
                       "-j{}".format(MAX_CPU_CORES),           # use multi-threading
                       "--cv",                                 # include cross-validation standard errors
                       "../bed/{0}.pruned.bed".format(self.group), # using this input file
                       self.k])                                # for K ancestral populations

        # restore previous working directory
        os.chdir('..')

        # save the log file
        with self.output()[2].open('w') as fout:
            fout.write(log)


class PlotAdmixtureK(luigi.Task):
    """
    Use ggplot to plot the admixture Q stats
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()
    k = luigi.IntParameter()

    def requires(self):
        return AdmixtureK(self.populations, self.genome, self.targets, self.group, self.k)

    def output(self):
        extensions = ['data', 'pdf']
        return [luigi.LocalTarget("admix/{0}.pruned.{1}.data".format(self.group, self.k)),
                luigi.LocalTarget("pdf/{0}.K.{1}.pdf".format(self.group, self.k))]

    def run(self):

        # use awk and paste to add population and sample names, needed for the plot
        awk = "awk '{ print substr($1, length($1)-2, 3) \"-\" substr($2, length($2)-2, 3) }' bed/" + str(self.group) + ".fam | " \
              "paste - admix/" + str(self.group) + ".pruned." + str(self.k) + ".Q"

        data = run_cmd([awk], True)

        # compose the header row
        header = ["Pop{}".format(i) for i in range(1, int(self.k) + 1)]
        header.insert(0, "Samples")

        # save the labeled file
        with self.output()[0].open('w') as fout:
            fout.write("\t".join(header)+"\n")
            fout.write(data)

        # generate a PDF of the admixture stacked column chart
        run_cmd(["Rscript",
                 "plot-admix-k.R",
                 self.output()[0].path,
                 self.output()[1].path])


class AdmixtureCV(luigi.Task):
    """
    Run admixture for the given population, determine the optimal K value, and plot the graphs
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        # run admixture or each population and each value of K
        for k in range(1, MAX_ANCESTRAL_K + 1):
            yield AdmixtureK(self.populations, self.genome, self.targets, self.group, k)

    def output(self):
        return [luigi.LocalTarget("admix/{0}.CV.data".format(self.group)),
                luigi.LocalTarget("pdf/{0}.CV.pdf".format(self.group)),]

    def run(self):

        # use grep to extract the cross-validation scores from all the log files
        cvs = run_cmd(["grep -h CV admix/{0}*.log".format(self.group)], True)

        # extract the K value and CV score
        # e.g. "CV error (K=1): 1.20340"
        data = [tuple([eval(re.sub(r'[^\d.]+', "", val)) for val in cv.split(":")]) for cv in cvs.splitlines()]

        # write the scores to a data file
        with self.output()[0].open('w') as fout:
            fout.write("\t".join(["K", "CV"]) + "\n")
            for row in data:
                fout.write("\t".join(str(datum) for datum in row) + "\n")

        # plot the CV values as a line graph
        run_cmd(["Rscript",
                 "plot-admix-cv.R",
                 self.output()[0].path,
                 self.output()[1].path])

        # get the three lowest CV scores
        bestfit = sorted(data, key=lambda x: x[1])[0:3]

        # plot the admixture percentages for the 3 best fitting values of k
        for k, cv in bestfit:
            yield PlotAdmixtureK(self.populations, self.genome, self.targets, self.group, int(k))


class FlashPCA(luigi.Task):
    """
    Run flashpca on the pruned BED files
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return PlinkPruneBed(self.populations, self.genome, self.targets, self.group)

    def output(self):
        prefixes = ['eigenvalues', 'eigenvectors', 'pcs', 'pve']
        return [luigi.LocalTarget("flashpca/{0}_{1}.txt".format(prefix, self.group)) for prefix in prefixes]

    def run(self):

        # flashpca only outputs to the current directory
        os.chdir('./flashpca')

        log = run_cmd(["flashpca",
                       "--bfile", "../bed/{0}.pruned".format(self.group),
                       "--numthreads", MAX_CPU_CORES,
                       "--suffix", "_{0}.txt".format(self.group)])

        # restore previous working directory
        os.chdir('..')


class PlotFlashPCA(luigi.Task):
    """
    Use ggplot to plot the PCA
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()
    pcs1 = luigi.IntParameter()
    pcs2 = luigi.IntParameter()
    labeled = luigi.BoolParameter(default=False)

    def requires(self):
        return FlashPCA(self.populations, self.genome, self.targets, self.group)

    def output(self):
        return [luigi.LocalTarget("flashpca/pca_{0}.data".format(self.group)),
                luigi.LocalTarget("flashpca/pve_{0}.txt".format(self.group)),
                luigi.LocalTarget("pdf/{0}.PCA.{1}.{2}.pdf".format(self.group, self.pcs1, self.pcs2))]

    def run(self):

        # use awk to add population and sample names, needed for the plot
        data = run_cmd(["awk '{print $1\"\t\"$2}' bed/" + str(self.group) + ".fam | "
                        "paste - flashpca/pcs_" + str(self.group) + ".txt"], True)

        # save the labeled file
        with self.output()[0].open('w') as fout:
            fout.write(data)

        # generate a PDF of the PCA plot
        run_cmd(["Rscript",
                 "plot-flashpca.R",
                 self.output()[0].path,      # pcs data
                 self.output()[1].path,      # pve data
                 self.output()[2].path,      # pdf location
                 self.pcs1,                  # component for x-axis
                 self.pcs2,                  # component for y-axis
                 1 if self.labeled else 0])  # label points (yes/no)


class PlotPhyloTree(luigi.Task):
    """
    Create a phylogenetic tree from a pruned BED file
    """
    populations = luigi.DictParameter()
    genome = luigi.Parameter()
    targets = luigi.Parameter()
    group = luigi.Parameter()

    def requires(self):
        return PlinkPruneBed(self.populations, self.genome, self.targets, self.group)

    def output(self):
        return [luigi.LocalTarget("tree/{0}.data".format(self.group)),
                luigi.LocalTarget("tree/{0}.tree".format(self.group)),
                luigi.LocalTarget("pdf/{0}.Phylo.pdf".format(self.group))]

    def run(self):

        # TODO what about bootstrapping?
        # make the distance matrix
        run_cmd(["plink",
                 "--distance", "square", "1-ibs",
                 "--bfile", "bed/{0}.pruned".format(self.group),
                 "--out", "tree/{0}".format(self.group)])

        # use awk to extract short versions of the pop and sample names
        awk = "awk '{ print substr($1, length($1)-2, 3) \"-\" substr($2, length($2)-2, 3) }' ./bed/" + str(self.group) + ".fam"

        # fetch the names as a row
        head = run_cmd([awk + " | xargs"], True)

        # add the samples names as a column to the mdist data
        data = run_cmd([awk + " | paste - ./bed/{0}.mdist".format(self.group)], True)

        # save the labeled file
        with self.output()[0].open('w') as fout:
            fout.write("\t"+head)
            fout.write(data)

        # generate a tree from the labeled data
        run_cmd(["Rscript",
                 "plot-phylo-tree.R",
                 self.output()[0].path,
                 self.output()[1].path,
                 self.output()[2].path])


class CustomGenomePipeline(luigi.Task):
    """
    Run all the samples through the pipeline
    """

    def complete(self):
        # always run the pipeline
        return False

    def requires(self):

        # make the SFS for dadi
        yield SiteFrequencySpectrum(POPULATIONS, GENOME, TARGETS)

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
            # run flashpca for each population (for the top 5 components)
            for pcs1, pcs2 in [(1,2), (3,4), (5,6)]:
                yield PlotFlashPCA(groups[group], GENOME, TARGETS, group, pcs1, pcs2, True)

        yield PlotPhyloTree(POPULATIONS, GENOME, TARGETS, 'all-pops')

        yield AdmixtureCV(groups['no-outgroup'], GENOME, TARGETS, 'no-outgroup')

if __name__=='__main__':
    luigi.run()