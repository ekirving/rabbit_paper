#!/usr/bin/python
import luigi, os, os.path, glob, random, re, shutil

# import my custom modules
from pipeline_dadi import *
from pipeline_utils import *


class SraToolsFastqDump(SetupTask):
    """
    Fetches the paired-end and single end FASTQ files for a given sample ID, using SRA Tools
    """
    sample = luigi.Parameter()

    def output(self):
        return [luigi.LocalTarget("fastq/{0}.fastq.gz".format(self.sample)),    # unpaired
                luigi.LocalTarget("fastq/{0}_1.fastq.gz".format(self.sample)),  # forward pair
                luigi.LocalTarget("fastq/{0}_2.fastq.gz".format(self.sample))]  # reverse pair

    def run(self):

        # use the NCBI SRA toolkit to fetch the fastq files
        run_cmd(["fastq-dump",
                 "--gzip",               # output gzipped files
                 "--split-3",            # split into two paired end fastq files + one unpaired fastq
                 "--outdir", "./fastq",  # output directory
                 self.sample])


class EnsemblReferenceFasta(SetupTask):
    """
    Fetches the reference genome FASTA file from Ensembl
    """
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("fasta/{0}.fa".format(self.genome))

    def run(self):

        # get the name of the zip file
        zipfile = os.path.basename(GENOME_URL)

        # TODO install ensembl perl API for a generic method of downloading fasta files

        # download the gzipped fasta file of the reference genome
        curl_download(GENOME_URL, zipfile)

        # unzip the file
        fasta = unzip_file(zipfile)

        with self.output().open('w') as fout:
            fout.write(fasta)


class SamtoolsFaidx(SetupTask):
    """
    Create the samtools fasta index for the reference genome, needed for GATK
    """
    genome = luigi.Parameter()

    def requires(self):
        return EnsemblReferenceFasta(self.genome)

    def output(self):
        return luigi.LocalTarget("{0}.fai".format(self.input().path))

    def run(self):
        run_cmd([SAMTOOLS,
                 "faidx",
                 self.input().path])  # fasta file


class PicardCreateSequenceDictionary(SetupTask):
    """
    Create the sequence dictionary for the reference genome, needed for GATK
    """
    genome = luigi.Parameter()

    def requires(self):
        return EnsemblReferenceFasta(self.genome)

    def output(self):
        filepath = trim_ext(self.input().path)
        return luigi.LocalTarget("{0}.dict".format(filepath))

    def run(self):
        run_cmd(["java", "-jar", PICARD,
                 "CreateSequenceDictionary",
                 "R={0}".format(self.input().path),    # reference fasta file
                 "O={0}".format(self.output().path)])  # dictionary file


class BwaIndexBwtsw(SetupTask):
    """
    Builds the BWA index for the reference genome, needed for performing alignments
    """
    genome = luigi.Parameter()

    def requires(self):
        return EnsemblReferenceFasta(self.genome)

    def output(self):
        extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
        return [luigi.LocalTarget("{0}.{1}".format(self.input().path, ext)) for ext in extensions]

    def run(self):

        run_cmd(["bwa",
                 "index",                              # index needed for bwa alignment
                 "-a", "bwtsw",                        # algorithm suitable for mammals
                 self.input().path])  # input file


class BwaMem(SetupTask):
    """
    Align the paired end fastq files to the reference genome
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()
    paired = luigi.BoolParameter()

    # TODO add resource constrants to bwa mem processes so the server doesn't crash
    # http://stackoverflow.com/questions/33429847/how-to-avoid-running-a-specific-task-simultaneously-in-luigi-with-multiple-worke

    def requires(self):
        return [SraToolsFastqDump(self.sample),
                BwaIndexBwtsw(self.genome)]

    def output(self):
        suffix = "paired" if self.paired else "single"
        return luigi.LocalTarget("sam/{0}.{1}.sam".format(self.sample, suffix))

    def run(self):

        read_group = "@RG\\tID:{sample}\\tSM:{sample}".format(sample=self.sample)

        # TODO how to resolve self.input() here
        if self.paired:
            fastq = ["fastq/{0}_1.fastq.gz".format(self.sample),  # pair 1
                     "fastq/{0}_2.fastq.gz".format(self.sample)]  # pair 2
        else :
            fastq = ["fastq/{0}.fastq.gz".format(self.sample)]    # unapired

        # perform the alignment
        sam = run_cmd(["bwa",
                       "mem",                               # align using the mem algorithm
                       "-t", MAX_CPU_CORES,                 # number of cores
                       "-R", read_group,                    # read group metadata
                       "fasta/{0}.fa".format(self.genome)]  # reference genome
                       + fastq,                             # input files
                      returnout=True)

        # save the SAM file
        with self.output().open('w') as fout:
            fout.write(sam)


class SamtoolsSortBam(SetupTask):
    """
    Convert an uncompressed SAM file into BAM format, and sort it
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()
    paired = luigi.BoolParameter()

    def requires(self):
        return BwaMem(self.sample, self.genome, self.paired)

    def output(self):
        filename = trim_path_ext(self.input().path)
        return luigi.LocalTarget("bam/{0}.bam".format(filename))

    def run(self):

        # perform the SAM -> BAM conversion and sorting
        bam = run_cmd([SAMTOOLS,
                       "sort",                    # sort the reads
                       "-l", DEFAULT_COMPRESSION, # level of compression
                       "-@", MAX_CPU_CORES,       # number of cores
                       "-O", "bam",               # output a BAM file
                       self.input().path],
                      returnout=True)

        # save the BAM file
        with self.output().open('w') as fout:
            fout.write(bam)


class SamtoolsMergeBam(SetupTask):
    """
    Merge the two SAM files from the paired and unpaired reads.
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        # we require both the paired and unpaired BAM files
        return [SamtoolsSortBam(self.sample, self.genome, True),
                SamtoolsSortBam(self.sample, self.genome, False)]

    def output(self):
        return luigi.LocalTarget("bam/{0}.bam".format(self.sample))

    def run(self):

        # get the bam files to merge
        in_files = [input.path for input in self.input()]

        # perform the SAM -> BAM conversion and sorting
        bam = run_cmd([SAMTOOLS,
                       "merge",              # merge the files
                       "-c",                 # combine the @RG headers
                       "-@", MAX_CPU_CORES,  # number of cores
                       self.output().path]   # output file
                       + in_files,           # input files
                      returnout=True)


class PicardMarkDuplicates(SetupTask):
    """
    Remove PCR duplicates, so we don't overestimate coverage
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return SamtoolsMergeBam(self.sample, self.genome)

    def output(self):
        filepath = trim_ext(self.input().path)
        return luigi.LocalTarget("{0}.rmdup.bam".format(filepath))

    def run(self):

        # TODO remove when confirmed unnecessary
        # # https://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation
        # run_cmd(["java", "-jar", PICARD,
        #          "FixMateInformation",
        #          "INPUT={0}".format(self.input().path)])

        run_cmd(["java", "-jar", PICARD,
                 "MarkDuplicates",
                 "INPUT={0}".format(self.input().path),
                 "OUTPUT={0}".format(self.output().path),
                 "METRICS_FILE={0}.txt".format(trim_ext(self.output().path)),
                 "REMOVE_DUPLICATES=true",
                 "QUIET=true"])


class SamtoolsIndexBam(SetupTask):
    """
    Create an index for the BAM file, for fast random access
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget("{0}.bai".format(self.input().path))

    def requires(self):
        return PicardMarkDuplicates(self.sample, self.genome)

    def run(self):

        run_cmd([SAMTOOLS,
                 "index",
                 "-b",                # create a BAI index
                 self.input().path])  # file to index


class GatkHaplotypeCaller(SetupTask):
    """
    Convert the BAM file into a gVCF
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        # reference must be indexed properly
        yield PicardCreateSequenceDictionary(self.genome)
        yield SamtoolsFaidx(self.genome)

        # bam file must also be indexed
        yield SamtoolsIndexBam(self.sample, self.genome)

    def output(self):
        # filename = trim_path_ext(self.input().path)
        return luigi.LocalTarget("vcf/{0}.g.vcf".format(self.sample))

    def run(self):

        # TODO how to resolve self.input()

        run_cmd(["java", "-Xmx8G", "-jar", GATK,
                 "-T", "HaplotypeCaller",                   # use the HaplotypeCaller to call variants
                 "-R", "fasta/{0}.fa".format(self.genome),  # the indexed reference genome
                 "--genotyping_mode", "DISCOVERY",          # variant discovery
                 "--emitRefConfidence", "GVCF",             # reference model emitted with condensed non-variant blocks
                 "--output_mode", "EMIT_ALL_SITES",         # produces calls at any callable site
                 "-L", TARGETS_LIST,                        # limit to the list of regions defined in the targets file
                 "-stand_emit_conf", "10",                  # min confidence threshold
                 "-stand_call_conf", MIN_GENOTYPE_QUAL,     # min call threshold
                 "-I", "bam/{0}.rmdup.bam".format(self.sample),
                 "-o", self.output().path])


class GatkGenotypeGVCFs(SetupTask):
    """
    Joint genotype all the samples in a population
    """
    population = luigi.Parameter()
    samples = luigi.ListParameter()
    genome = luigi.Parameter()

    def requires(self):
        # samples must be individually called
        for sample in self.samples:
            yield GatkHaplotypeCaller(sample, self.genome)

    def output(self):
        return luigi.LocalTarget("vcf/{0}.vcf".format(self.population))

    def run(self):
        # TODO how to resolve self.input()

        # make a list of input files
        vcf_files = sum([["-V", "vcf/{0}.g.vcf".format(sample)] for sample in self.samples], [])

        run_cmd(["java", "-Xmx8G", "-jar", GATK,
                 "-T", "GenotypeGVCFs",                     # use GenotypeGVCFs to jointly call variants
                 "--num_threads", MAX_CPU_CORES,            # number of data threads to allocate to this analysis
                 "--includeNonVariantSites",                # include loci found to be non-variant after genotyping
                 "-R", "fasta/{0}.fa".format(self.genome),  # the indexed reference genome
                 "-o", self.output().path]
                + vcf_files)


class GatkSelectVariants(SetupTask):
    """
    Select only biallelic variant sites for downstream analysis
    """
    population = luigi.Parameter()
    samples = luigi.ListParameter()
    genome = luigi.Parameter()

    def requires(self):
        # population must have been joint called
        return GatkGenotypeGVCFs(self.population, self.samples, self.genome)

    def output(self):
        filepath = trim_ext(self.input().path)
        return luigi.LocalTarget("{0}.variant.vcf".format(filepath))

    def run(self):

        run_cmd(["java", "-Xmx2g", "-jar", GATK,
                 "-T", "SelectVariants",
                 "--selectTypeToInclude", "SNP",
                 "--restrictAllelesTo", "BIALLELIC",
                 "-R", "fasta/{0}.fa".format(self.genome), # TODO how can I fetch this from downstream
                 "-V", self.input().path,
                 "-o", self.output().path])

        # because SelectVariants doesn't handle having "*" in the ALT column we need to grep all those sites out
        vcf = run_cmd(["grep -Pv '\t\*\t' " + str(self.output().path)], shell=True, returnout=True)

        # overwrite the vcf with the filtered data
        with self.output().open('w') as fout:
            fout.write(vcf)


class PlinkMakeBed(SetupTask):
    """
    Convert a population gVCF into a BED file (binary pedigree)
    """
    population = luigi.Parameter()
    samples = luigi.ListParameter()
    genome = luigi.Parameter()

    def requires(self):
        yield GatkSelectVariants(self.population, self.samples, self.genome)

    def output(self):
        filename = trim_path_ext(self.input().path)
        return [luigi.LocalTarget("ped/{0}.ped".format(filename)),
                luigi.LocalTarget("bed/{0}.bed".format(filename))]

    def run(self):

        # plink requires extensionless filenames
        pedname = trim_ext(self.output()[0])
        bedname = trim_ext(self.output()[1])

        # convert VCF to PED
        run_cmd(["plink",
                 "--recode",
                 "--const-fid", self.population,  # use population name as "family" id
                 "--biallelic-only", "strict",
                 "--vcf-min-qual", MIN_GENOTYPE_QUAL,
                 "--vcf-require-gt",
                 "--vcf", self.input().path,
                 "--out", pedname])

        # plink requres that all sites have a unique variant name, we we'll have to bodge the map file a little
        mapfile = "{0}.map".format(pedname)

        # use awk to add variant IDs, so we can identify polyallelic sites during merge
        data = run_cmd(["awk '$2=$1\"-\"$4' " + mapfile], returnout=True, shell=True)

        # replace the old map file
        with open(mapfile, 'w') as fout:
            fout.write(data)

        # convert PED to BED
        run_cmd(["plink",
                 "--make-bed",
                 "--file", pedname,
                 "--out", bedname])


class PlinkMergeBeds(SetupTask):
    """
    Merge multiple BED files into one
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        for population, samples in GROUPS[self.group].iteritems():
            yield PlinkMakeBed(population, samples, self.genome)

    def output(self):
        return luigi.LocalTarget("bed/{0}.bed".format(self.group))

    def run(self):

        # generate a unique suffix for temporary files
        suffix = 'tmp' + str(random.getrandbits(100))

        # make a copy of the bed files because we'll need to filter them
        for population in GROUPS[self.group]:
            for ext in ['bed', 'bim', 'fam']:
                shutil.copyfile("bed/{0}.{1}".format(population, ext),
                                "bed/{0}.{1}.{2}".format(population, suffix, ext))

        # merge requires the first bed file to be named in the command
        bed_files = ["bed/{0}.{1}".format(population, suffix) for population in GROUPS[self.group]]
        bed_file1 = bed_files.pop(0)

        # make the merge-list with the remaining BED files
        with open("bed/{0}.list".format(self.group), 'w') as fout:
            fout.write("\n".join(bed_files))

        # TODO make --geno a separate task and add to filename

        # compose the merge command, because we are going to need it twice (retains ~8.6k variants)
        merge = ["plink",
                 "--make-bed",
                 "--geno", '0.1',       # remove all sites with more than 10% of genotypes missing (i.e. force 90% geno)
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
                for population in GROUPS[self.group]:
                    run_cmd(["plink",
                             "--make-bed",
                             "--exclude", "bed/{0}-merge.missnp".format(self.group),
                             "--bfile", "bed/{0}.{1}".format(population, suffix),
                             "--out", "bed/{0}.{1}".format(population, suffix)])

                # reattempt the merge
                run_cmd(merge)

            else:
                raise Exception(e)

        # tidy up all the temporary intermediate files
        for tmp in glob.glob("./bed/*{}*".format(suffix)):
            os.remove(tmp)


class PlinkIndepPairwise(SetupTask):
    """
    Produce a list of SNPs with high discriminating power, by filtering out minor allele frequency and sites under
    linkage disequilibrium
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return PlinkMergeBeds(self.group, self.genome)

    def output(self):
        extensions = ['in', 'out']
        return [luigi.LocalTarget("bed/{0}.prune.{1}".format(self.group, ext)) for ext in extensions]

    def run(self):

        # TODO make --maf a separate task and add to filename
        # apply a filter on minor allele frequencies
        # retains ~94k SNPs (from ~123k)
        run_cmd(["plink",
                 "--make-bed",
                 "--maf", "0.10",  # allele must be observed > 10% of the time
                 "--bfile", "bed/{0}".format(self.group),
                 "--out", "bed/{0}.maf-0.1".format(self.group)])

        # TODO make --indep-pairwise a separate task and add to filename
        # calculate the prune list (prune.in / prune.out)
        # retains ~12k SNPs (from ~94k)
        run_cmd(["plink",
                 "--indep-pairwise", 50, 10, 0.5,  # accept R^2 coefficient of up to 0.5
                 "--bfile", "bed/{0}.maf-0.1".format(self.group),
                 "--out", "bed/{0}".format(self.group)])


class PlinkPruneBed(SetupTask):
    """
    Prune the merged group BED file using the prune.in list from PlinkIndepPairwise()
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return PlinkIndepPairwise(self.group, self.genome)

    def output(self):
        filepath = trim_ext(self.input().path)
        return luigi.LocalTarget("{0}.pruned.bed".format(filepath))

    def run(self):

        # TODO fix inputs
        # apply the prune list
        run_cmd(["plink",
                 "--make-bed",
                 "--extract", "bed/{0}.prune.in".format(self.group),
                 "--bfile", "bed/{0}".format(self.group),
                 "--out", trim_ext(self.output().path)])


class AdmixtureK(SetupTask):
    """
    Run admixture, with K ancestral populations, on the pruned BED files
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()
    k = luigi.IntParameter()

    def requires(self):
        return PlinkPruneBed(self.group, self.genome)

    def output(self):
        filename = trim_path_ext(self.input().path)
        return luigi.LocalTarget("admix/{0}.{1}.Q".format(filename, self.k))

    def run(self):

        # admixture only outputs to the current directory
        os.chdir('./admix')

        # TODO what about bootstrapping (-B / -B2000)
        log = run_cmd(["admixture",
                       "-j{0}".format(MAX_CPU_CORES),       # use multi-threading
                       "--cv",                              # include cross-validation standard errors
                       "../{0}".format(self.input().path),  # using this input file
                       self.k],                             # for K ancestral populations
                      returnout=True, pwd='../')

        # restore previous working directory
        os.chdir('..')

        # TODO log everything like this!
        log_file = replace_ext(self.output().path, "log")

        # save the log file
        with open(log_file, 'w') as fout:
            fout.write(log)


class AdmixturePlotK(SetupTask):
    """
    Use ggplot to plot the admixture Q stats
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()
    k = luigi.IntParameter()

    def requires(self):
        return AdmixtureK(self.group, self.genome, self.k)

    def output(self):
        filename = trim_path_ext(self.input().path)
        return luigi.LocalTarget("pdf/admix.{0}.pdf".format(filename))

    def run(self):

        # TODO what about the ".fam" file?
        # use awk and paste to add population and sample names, needed for the plot
        awk = "awk '{ print substr($1, length($1)-2, 3) \"-\" substr($2, length($2)-2, 3) }' bed/" + str(self.group) + ".fam | " \
              "paste - " + self.input().path

        data = run_cmd([awk], returnout=True, shell=True)

        # compose the header row
        header = ["Pop{}".format(i) for i in range(1, int(self.k) + 1)]
        header.insert(0, "Samples")

        labeled_file = replace_ext(self.input().path, "data")

        # save the labeled file
        with open(labeled_file, 'w') as fout:
            fout.write("\t".join(header)+"\n")
            fout.write(data)

        # generate a PDF of the admixture stacked column chart
        run_cmd(["Rscript",
                 "rscript/plot-admix-k.R",
                 labeled_file,
                 self.output().path])


class AdmixtureCV(SetupTask):
    """
    Run admixture for the given population, determine the optimal K value, and plot the graphs
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        # run admixture or each population and each value of K
        for k in range(1, MAX_ANCESTRAL_K + 1):
            yield AdmixturePlotK(self.group, self.genome, k)

    def output(self):
        filename = trim_path_ext(trim_ext(self.input()[0].path))
        return luigi.LocalTarget("pdf/admix.{0}.CV.pdf".format(filename))

    def run(self):

        # TODO replace dependence on filename pattern
        # use grep to extract the cross-validation scores from all the log files
        cvs = run_cmd(["grep -h CV admix/{0}*.log".format(self.group)], returnout=True, shell=True)

        # extract the K value and CV score
        # e.g. "CV error (K=1): 1.20340"
        data = [tuple([re.sub(r'[^\d.]+', "", val) for val in cv.split(":")]) for cv in cvs.splitlines()]

        # get the three lowest CV scores
        # bestfit = sorted(data, key=lambda x: x[1])[0:3]

        cvdata_file = "{0}.CV.data".format(trim_ext(trim_ext(self.input()[0].path)))

        # write the scores to a data file
        with open(cvdata_file, 'w') as fout:
            fout.write("\t".join(["K", "CV"]) + "\n")
            for row in data:
                fout.write("\t".join(str(datum) for datum in row) + "\n")

        # plot the CV values as a line graph
        run_cmd(["Rscript",
                 "rscript/plot-line-graph.R",
                 cvdata_file,
                 self.output().path,
                 "Ancestral populations (K)",
                 "Cross-validation Error"])


class sNMF_Ped2Geno(SetupTask):
    """
    Convert from BED to PED to Geno file type
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return PlinkPruneBed(self.group, self.genome)

    def output(self):
        filename = trim_path_ext(self.input().path)
        return luigi.LocalTarget("snmf/{0}.geno".format(filename))

    def run(self):

        # plink requires extensionless filenames
        filename = trim_path_ext(self.input().path)
        bedname = "bed/{0}".format(filename)
        pedname = "ped/{0}".format(filename)

        # convert from BED to PED (because sNMF can't handle BED files)
        run_cmd(["plink",
                 "--recode", "12",
                 "--bfile", bedname,
                 "--out", pedname])

        # convert from PED to GENO
        run_cmd(["ped2geno",
                 "{0}.ped".format(pedname),
                 self.output().path])


class sNMF_K(SetupTask):
    """
    Run sNMF, with K ancestral populations, on the pruned BED files
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()
    k = luigi.IntParameter()

    def requires(self):
        return sNMF_Ped2Geno(self.group, self.genome)

    def output(self):
        filename = trim_path_ext(self.input().path)
        return luigi.LocalTarget("snmf/{0}.{1}.Q".format(filename, self.k))

    def run(self):

        # run sNMF on the geno file, and compute the cross-entropy criterion
        log = run_cmd(["sNMF",
                       "-p", MAX_CPU_CORES,
                       "-x", self.input().path,
                       "-K", self.k,
                       "-c"],
                      returnout=True)

        log_file = replace_ext(self.input().path, "log")

        # save the log file
        with open(log_file, 'w') as fout:
            fout.write(log)


class sNMF_PlotK(SetupTask):
    """
    Use ggplot to plot the admixture Q stats
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()
    k = luigi.IntParameter()

    def requires(self):
        return sNMF_K(self.group, self.genome, self.k)

    def output(self):
        filename = trim_path_ext(self.input().path)
        return luigi.LocalTarget("pdf/snmf.{0}.pdf".format(filename))

    def run(self):

        # TODO fam file dependency
        # use awk and paste to add population and sample names, needed for the plot
        awk = "awk '{ print substr($1, length($1)-2, 3) \"-\" substr($2, length($2)-2, 3) }' bed/" + str(self.group) + ".fam | " \
              "paste - " + self.input().path

        data = run_cmd([awk], returnout=True, shell=True)

        # compose the header row
        header = ["Pop{}".format(i) for i in range(1, int(self.k) + 1)]
        header.insert(0, "Samples")

        labeled_file = replace_ext(self.input().path, "data")

        # save the labeled file
        with open(labeled_file, 'w') as fout:
            fout.write("\t".join(header)+"\n")
            fout.write(data)

        # generate a PDF of the admixture stacked column chart
        run_cmd(["Rscript",
                 "rscript/plot-admix-k.R",
                 labeled_file,
                 self.output().path])


class sNMF_CE(SetupTask):
    """
    Run sNMF for the given population, determine the optimal K value, and plot the graphs
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        # run admixture or each population and each value of K
        for k in range(1, MAX_ANCESTRAL_K + 1):
            yield sNMF_PlotK(self.group, self.genome, k)

    def output(self):
        filename = trim_path_ext(trim_ext(self.input()[0].path))
        return luigi.LocalTarget("pdf/sNMF.{0}.CV.pdf".format(filename))

    def run(self):

        # key phrase to look for in log files
        keyphrase = "Cross-Entropy (masked data):"

        # TODO replace dependence on filename pattern
        # use grep to find the lines containing the cross-entropy scores from all the log files
        ces = run_cmd(["grep '" + keyphrase + "' snmf/{0}*.log".format(self.group)], returnout=True, shell=True)

        # extract the K value and cross-entropy score
        data = [tuple([re.sub(r'[^\d.]+', "", val).strip('.') for val in ce.split(keyphrase)]) for ce in ces.splitlines()]

        cedata_file = "{0}.CV.data".format(trim_ext(trim_ext(self.input()[0].path)))

        # write the scores to a data file
        with open(cedata_file, 'w') as fout:
            fout.write("\t".join(["K", "Cross-Entropy"]) + "\n")
            for row in data:
                fout.write("\t".join(str(datum) for datum in row) + "\n")

        # plot the CV values as a line graph
        run_cmd(["Rscript",
                 "rscript/plot-line-graph.R",
                 cedata_file,
                 self.output().path,
                 "Ancestral populations (K)",
                 "Cross-Entropy (masked data)"])


class FlashPCA(SetupTask):
    """
    Run flashpca on the pruned BED files
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return PlinkMergeBeds(self.group, self.genome)

    def output(self):
        prefixes = ['pcs', 'pve']
        filename = trim_path_ext(self.input().path)
        return [luigi.LocalTarget("flashpca/{0}_{1}.txt".format(prefix, filename)) for prefix in prefixes]

    def run(self):

        # flashpca only outputs to the current directory
        os.chdir('./flashpca')

        filename = trim_path_ext(self.input().path)

        run_cmd(["flashpca",
                 "--bfile", "../{0}".format(self.input().path),
                 "--numthreads", MAX_CPU_CORES,
                 "--suffix", "_{0}.txt".format(filename)],
                pwd='../')

        # restore previous working directory
        os.chdir('..')


class FlashPCAPlot(SetupTask):
    """
    Use ggplot to plot the PCA
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return FlashPCA(self.group, self.genome)

    def output(self):

        yield luigi.LocalTarget("flashpca/pca_{0}.data".format(self.group))
        yield luigi.LocalTarget("flashpca/pve_{0}.txt".format(self.group))

        for pcs1, pcs2 in [(1, 2), (3, 4), (5, 6)]:
            yield luigi.LocalTarget("pdf/{0}.PCA.{1}.{2}.pdf".format(self.group, pcs1, pcs2))

    def run(self):

        # use awk to add population and sample names, needed for the plot
        data = run_cmd(["awk '{print $1\"\t\"$2}' bed/" + str(self.group) + ".fam | "
                        "paste - flashpca/pcs_" + str(self.group) + ".txt"], returnout=True, shell=True)

        # save the labeled file
        with open("flashpca/pca_{0}.data".format(self.group), 'w') as fout:
            fout.write(data)

        # plot the first 6 components
        for pcs1, pcs2 in [(1, 2), (3, 4), (5, 6)]:

            # generate both labeled and unlabeled PDFs
            pdfs = {
                0: "pdf/{0}.PCA.{1}.{2}.pdf".format(self.group, pcs1, pcs2),
                1: "pdf/{0}.PCA.{1}.{2}.labeled.pdf".format(self.group, pcs1, pcs2)
            }

            for labeled, pdf_path in pdfs.iteritems():
                # generate a PDF of the PCA plot
                run_cmd(["Rscript",
                         "rscript/plot-flashpca.R",
                         "flashpca/pca_{0}.data".format(self.group),  # pca data
                         "flashpca/pve_{0}.txt".format(self.group),   # pve data
                         pdf_path,                                    # pdf location
                         pcs1,                                        # component for x-axis
                         pcs2,                                        # component for y-axis
                         labeled])                                    # show point labels (0/1)


class APE_PlotTree(SetupTask):
    """
    Create a phylogenetic tree from a pruned BED file
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return PlinkMergeBeds(self.group, self.genome)

    def output(self):
        return [luigi.LocalTarget("tree/{0}.data".format(self.group)),
                luigi.LocalTarget("tree/{0}.tree".format(self.group)),
                luigi.LocalTarget("pdf/{0}.Phylo.pdf".format(self.group))]

    def run(self):

        # TODO what about bootstrapping?
        # make the distance matrix
        run_cmd(["plink",
                 "--distance", "square", "1-ibs",
                 "--bfile", "bed/{0}".format(self.group),
                 "--out", "tree/{0}".format(self.group)])

        # use awk to extract short versions of the pop and sample names
        awk = "awk '{ print substr($1, length($1)-2, 3) \"-\" substr($2, length($2)-2, 3) }' ./bed/" + str(self.group) + ".fam"

        # fetch the names as a row
        head = run_cmd([awk + " | xargs"], returnout=True, shell=True)

        # add the samples names as a column to the mdist data
        data = run_cmd([awk + " | paste - ./tree/{0}.mdist".format(self.group)], returnout=True, shell=True)

        # save the labeled file
        with self.output()[0].open('w') as fout:
            fout.write("\t"+head)
            fout.write(data)

        # generate a tree from the labeled data
        run_cmd(["Rscript",
                 "rscript/plot-phylo-tree.R",
                 self.output()[0].path,
                 self.output()[1].path,
                 self.output()[2].path])


class CustomGenomePipeline(luigi.WrapperTask):
    """
    Run all the samples through the pipeline
    """

    def requires(self):

        # run admixture (on pruned data)
        yield AdmixtureCV('no-outgroup', GENOME)

        # run sNMF (on pruned data)
        yield sNMF_CE('no-outgroup', GENOME)

        # plot a phylogenetic tree
        yield APE_PlotTree('all-pops', GENOME)

        # run flashpca for each population (for the top 6 components)
        for group in GROUPS:
            yield FlashPCAPlot(group, GENOME)


if __name__=='__main__':
    luigi.run()