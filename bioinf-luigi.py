from subprocess import check_output
import subprocess
import luigi
import os
import random

# downloade from http://coderscrowd.com/app/public/codes/view/229

#faking the pipline for now

sample_list = ["GP140-10A_S19","GP140-10B_S20","GP140-11A_S21","GP140-11B_S22","GP140-1A_S1","GP140-1B_S2","GP140-2A_S3","GP140-2B_S4","GP140-3A_S5","GP140-3B_S6","GP140-4A_S7","GP140-4B_S8","GP140-5A_S9","GP140-5B_S10","GP140-6A_S11","GP140-6B_S12","GP140-7A_S13","GP140-7B_S14","GP140-8A_S15","GP140-8B_S16","GP140-9A_S17","GP140-9B_S18"]

genome = "reference/gp140.fa"


def run_cmd(cmd):
    """
    Executes the given command in a system subprocess
    """
    proc = subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=subprocess.PIPE)
    (output, error) = proc.communicate()

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
        # download the gzipped fastq files from the EBI ftp server
        url1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_1.fastq.gz".format(sample=self.sample)
        url2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR997/{sample}/{sample}_2.fastq.gz".format(sample=self.sample)

        return [Curl_Download(url1, "fastq/%s_1.fastq.gz" % self.sample),
                Curl_Download(url2, "fastq/%s_2.fastq.gz" % self.sample)]

    def output(self):
        return [luigi.LocalTarget("fastq/%s_1.fastq" % self.sample),
                luigi.LocalTarget("fastq/%s_2.fastq" % self.sample)]

    def run(self):
        # unzip the files
        (file1, file2) = self.output()
        run_cmd(["gunzip", "-k", file1.path, file2.path])

        print "====== Downloaded Paired-end FASTQ ======"

class Bwa_Index(luigi.Task):
    """
    Builds the BWA index for the given fasta file
    """
    sample = luigi.Parameter()

    def requires(self):
        return

    def output(self):
        return luigi.LocalTarget("fasta/%s.fa" % self.genome),

    def run(self):
        return

class Bwa_Mem(luigi.Task):

    # Ideally we have to show this but for now I am gonna fake it by creating a list of sample IDs
    # fastq_path = luigi.Parameter()

    sample = luigi.Parameter()

    def requires(self):
       return PairedEnd_Fastq(self.sample)

    def output(self):
        return luigi.LocalTarget("sam/%s.sam" % self.sample)

    def run(self):
        tmp = run_cmd(["bwa",
                     "mem",
                     genome,
                     "fastq/%s_1.fastq" % self.sample,
                     "fastq/%s_2.fastq" % self.sample])
        with self.output().open('w') as sam_out:
            sam_out.write(tmp)
        print "====== Running BWA ======"


class Convert_Sam_Bam(luigi.Task):

    sample = luigi.Parameter()
    #genome = luigi.Parameter()

    def requires(self):
        return Bwa_Mem(self.sample)

    def output(self):
        return luigi.LocalTarget("bam/%s.bam" % self.sample)

    def run(self):
        tmp = run_cmd(["samtools",
                       "view",
                       "-bS",
                       "sam/"+self.sample+".sam"])
        with self.output().open('w') as bam_out:
            bam_out.write(tmp)
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