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
        return luigi.LocalTarget("fsdata/{0}.data".format(self.group))

    def run(self):

        # log everything to file
        logging.basicConfig(filename="fsdata/{0}.log".format(self.group), level=logging.DEBUG)

        # generate the frequency spectrum
        fsdata = generate_frequency_spectrum(GROUPS[self.group])

        # save the fsdata file
        with self.output().open('w') as fout:
            fout.write(fsdata)


class DadiSpectrum(luigi.Task):
    """
    Generate the dadi Spectrum file for the two populations
    """
    group = luigi.Parameter()
    pop1 = luigi.Parameter()
    pop2 = luigi.Parameter()

    def requires(self):
        return SiteFrequencySpectrum(self.group, GENOME)

    def output(self):
        return luigi.LocalTarget("fsdata/{0}_{1}.fs".format(self.pop1, self.pop2))

    def run(self):

        # parse the data file to generate the data dictionary
        dd = dadi.Misc.make_data_dict("fsdata/{0}.data".format(self.group))

        # get the two populations
        pops = [self.pop1, self.pop2]

        # project each population down to allow for incomplete coverage
        prj = [round(len(POPULATIONS[pop]) * DADI_PROJECT, 0) for pop in pops]

        # extract the spectrum for the two populations from the dictionary, and project down to 80% of sample size
        fs = dadi.Spectrum.from_data_dict(dd, pops, prj, polarized=True)

        # save it to a file
        fs.to_file(self.output().path)


class DadiOptimizeLogParams(luigi.Task):
    """
    Optimise the paramaters for the given model
    """
    pop1 = luigi.Parameter()
    pop2 = luigi.Parameter()
    # model = luigi.Parameter()
    n = luigi.IntParameter()

    def requires(self):
        return DadiSpectrum('all-pops', self.pop1, self.pop2)

    def output(self):
        # TODO make model a param
        return luigi.LocalTarget("fsdata/{0}.{1}.opt".format("split_mig", self.n))

    def run(self):

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file("fsdata/{0}_{1}.fs".format(self.pop1, self.pop2))

        # These are the grid point settings will use for extrapolation.
        pts_l = [10, 50, 60]

        # The Demographics1D and Demographics2D modules contain a few simple models,
        # mostly as examples. We could use one of those.
        func = dadi.Demographics2D.split_mig

        # The upper_bound and lower_bound lists are for use in optimization.
        upper_bound = [100, 100, 3, 10]
        lower_bound = [1e-4, 1e-4, 0, 0]

        # TODO randomise these starting values
        # This is our initial guess for the parameters, which is somewhat arbitrary.
        p0 = [2, 0.1, 0.2, 0.2]

        # Make the extrapolating version of our demographic model function.
        func_ex = dadi.Numerics.make_extrap_log_func(func)

        # Perturb our parameters before optimization. This does so by taking each
        # parameter a up to a factor of two up or down.
        p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

        # Do the optimization...
        popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                           lower_bound=lower_bound,
                                           upper_bound=upper_bound,
                                           # verbose=len(p0),
                                           maxiter=DADI_MAX_ITER)

        ns = fs.sample_sizes

        # Calculate the best-fit model AFS.
        model = func_ex(popt, ns, pts_l)

        # Likelihood of the data given the model AFS.
        ll_model = dadi.Inference.ll_multinom(model, fs)

        # The optimal value of theta given the model.
        theta = dadi.Inference.optimal_sfs_scaling(model, fs)

        # save the output
        with self.output().open('w') as fout:
            fout.write('Maximum log composite likelihood: {0}\n'.format(ll_model))
            fout.write('Best-fit parameters: {0}\n'.format(popt))
            fout.write('Optimal value of theta: {0}\n'.format(theta))



class CustomDadiPipeline(luigi.WrapperTask):
    """
    Run the dadi models
    """

    def requires(self):

        for n in range(0, 10):
            yield DadiOptimizeLogParams('DOM', 'WLD-FRE', n)
