#!/usr/bin/python
import luigi, dadi, numpy

from pipeline_gatk import *
from pipeline_utils import *


class SiteFrequencySpectrum(luigi.Task):
    """
    Produce the site frequency spectrum, based on genotype calls from GATK GenotypeGVCFs
    """
    group = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        for population, samples in GROUPS[self.group].iteritems():
            yield GatkGenotypeGVCFs(population, samples, self.genome)

    def output(self):
        return luigi.LocalTarget("fsdata/{0}.data".format(self.genome))

    def run(self):

        # log everything to file
        logging.basicConfig(filename="fsdata/{0}.log".format(self.genome), level=logging.DEBUG)

        # generate the frequency spectrum
        fsdata = generate_frequency_spectrum(GROUPS[self.group])

        # save the fsdata file
        with self.output().open('w') as fout:
            fout.write(fsdata)


class DadiOptimizeLog(luigi.Task):
    """

    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget("".format())

    def run(self):

        data = run_cmd([""], returnout=True)

        # save the output
        with self.output().open('w') as fout:
            fout.write(data)


class CustomDadiPipeline(luigi.WrapperTask):
    """
    Run the dadi models
    """

    def requires(self):

        # make the SFS for dadi
        yield SiteFrequencySpectrum('all-pops', GENOME)

