from subprocess import check_output
import subprocess
import luigi
import os
import random

# downloaded from http://coderscrowd.com/app/public/codes/view/229

# the URL for the reference genome
genome = "OryCun2.0"
genome_url = "ftp://ftp.ensembl.org/pub/release-84/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz"

# the URLs for the paried end fastq files
pair1_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_1.fastq.gz"
pair2_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_2.fastq.gz"

def run_cmd(cmd):
    """
    Executes the given command in a system subprocess
    """
    print cmd
    proc = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    (output, error) = proc.communicate()
    if error:
        print error
        exit()
    return output

class Curl_Download(luigi.Task):
    """
    Downloads a remote url to a local file path using cURL
    """
    url = luigi.Parameter()
    file = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.file)

    def run(self):
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
        return [Curl_Download(pair1_url.format(sample=self.sample), "fastq/%s_1.fastq.gz" % self.sample),
                Curl_Download(pair2_url.format(sample=self.sample), "fastq/%s_2.fastq.gz" % self.sample)]

    def output(self):
        return [luigi.LocalTarget("fastq/%s_1.fastq" % self.sample),
                luigi.LocalTarget("fastq/%s_2.fastq" % self.sample)]

    def run(self):
        # unzip the files
        (file1, file2) = self.output()
        run_cmd(["gunzip", "-k", file1.path, file2.path])

        print "====== Downloaded Paired-end FASTQ ======"

class Genome_Fasta(luigi.Task):
    """
    Fetches the reference genome
    """
    genome = luigi.Parameter()
    # url = luigi.Parameter()

    def requires(self):
        global genome_url
        # download the gzipped fasta file of the reference genome
        return Curl_Download(genome_url, "fasta/%s.fa.gz" % self.genome)

    def output(self):
        return luigi.LocalTarget("fasta/%s.fa" % self.genome)

    def run(self):
        # unzip the file
        run_cmd(["gunzip", "-k", self.output().path])

        print "====== Downloaded genome FASTA ======"

class Bwa_Index(luigi.Task):
    """
    Builds the BWA index for the reference genome
    """
    genome = luigi.Parameter()

    def requires(self):
        return Genome_Fasta(self.genome)

    def output(self):
        files = []
        for ext in ['amb', 'ann', 'bwt', 'pac', 'sa']: # 'fai',
            files.append(luigi.LocalTarget("fasta/{genome}.fa.{ext}".format(genome=self.genome, ext=ext)))
        return files

    def run(self):
        run_cmd(["bwa",
                 "index",
                 "-a bwtsw",
                 "fasta/%s.fa" % self.genome])
        print "====== Built the BWA index ======"

class Bwa_Mem(luigi.Task):
    """
    Align the fastq files to the reference genome
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return [PairedEnd_Fastq(self.sample), Bwa_Index(self.genome)]

    def output(self):
        return luigi.LocalTarget("sam/%s.sam" % self.sample)

    def run(self):
        # TODO add the sample name to the readgroup header
        # TODO paramaterise the number of cores
        # perform the alignment
        sam = run_cmd(["bwa",
                       "mem",
                       "-t 8",
                       '-R "@RG\tID:{sample}\tSM:{sample}\t"'.format(sample=self.sample),
                       "fasta/%s.fa" % self.genome,
                       "fastq/%s_1.fastq" % self.sample,
                       "fastq/%s_2.fastq" % self.sample])

        # save the SAM file
        with self.output().open('w') as file:
            file.write(sam)

        print "====== Aligned the FASTQ using BWA ======"


class Convert_Sam_Bam(luigi.Task):
    """
    Convert an uncompressed SAM file into BAM format
    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Bwa_Mem(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("bam/%s.bam" % self.sample)

    def run(self):
        # perform the SAM -> BAM conversion
        bam = run_cmd(["samtools", "view", "-bS", "sam/%s.sam" % self.sample])

        # save the BAM file
        with self.output().open('w') as bam_out:
            bam_out.write(bam)

        print "===== Converting Sam file to Bam file ======"

class Sort_Bam(luigi.Task):

    sample = luigi.Parameter()

    def requires(self):
        return Convert_Sam_Bam(self.sample)

    def run(self):
        tmp = run_cmd(["samtools",
                       "sort",
                       "bam/"+self.sample+".bam",
                       "bam/"+self.sample+".sorted"])

        print "===== Sorting Bam ======="

class Index_Bam(luigi.Task):

    sample = luigi.Parameter()

    def output(self):
        print "Indexing ... Done"

    def requires(self):
        return Sort_Bam(self.sample)

    def run(self):
        tmp = run_cmd(["samtools",
                       "index",
                       "bam/"+self.sample+".sorted.bam"])

        print "==== Indexing Bam ===="


class Call_Variant(luigi.Task):

    sample = luigi.Parameter()

    def requires(self):
        return Index_Bam(self.sample)

    def output(self):
        return luigi.LocalTarget("bcf/%s.bcf" % self.sample)

    def run(self):
        tmp = run_cmd(["samtools",
                       "mpileup",
                       "-u",
                       "-f",
                       genome,
                       "bam/"+self.sample+".sorted.bam"])

        with self.output().open("w") as bcf_file:
            bcf_file.write(tmp)

        print "=====  Calling Variants ======="


class Convert_Bcf_Vcf(luigi.Task):

    sample = luigi.Parameter()

    def requires(self):
        return Call_Variant(self.sample)

    def output(self):
        return luigi.LocalTarget("vcf/%s.vcf" % self.sample)

    def run(self):
        tmp = run_cmd(["bcftools",
                       "view",
                       "-v",
                       "-c",
                       "-g",
                       "bcf/"+self.sample+".bcf"])

        with self.output().open("w") as vcf_file:
            vcf_file.write(tmp)

        print "===== Creating VCF ======="


class Custom_Genome_Pipeline(luigi.Task):

    def requires(self):
        return [Convert_Bcf_Vcf(sample) for sample in sample_list]

    def output(self):
        return luigi.LocalTarget("log.txt")

    def run(self):
        print "running..."
        #for sample in sample_list:
            #print sample+"\n"
            #return Convert_Bcf_Vcf(sample)


if __name__=='__main__':
    luigi.run()