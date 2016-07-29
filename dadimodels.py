#!/usr/bin/python
import luigi, dadi, numpy

from utilities import *


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
