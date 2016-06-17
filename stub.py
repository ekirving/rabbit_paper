import luigi

class Some_Task(luigi.Task):
    """

    """
    sample = luigi.Parameter()
    genome = luigi.Parameter()

    def requires(self):
        return Index_Cram(self.sample, self.genome)

    def output(self):
        return luigi.LocalTarget("pileup/"+self.sample+".pileup")

    def run(self):
        stdout = run_cmd([])          #

        # save the pileup file
        with self.output().open('w') as fout:
            fout.write(stdout)

        print "=====  ======="

