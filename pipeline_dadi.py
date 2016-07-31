#!/usr/bin/python

# load matplotlib before dadi so we can disable the screen
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import luigi, dadi, numpy, pylab, random

# import the custom pipelines
from pipeline_gatk import *
from pipeline_utils import *


class SiteFrequencySpectrum(SetupTask):
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


class DadiSpectrum(SetupTask):
    """
    Generate the dadi Spectrum file for the two populations
    """
    group = luigi.Parameter()
    pop1 = luigi.Parameter()
    pop2 = luigi.Parameter()

    def setup(self):
        # get the two populations
        self.pops = [self.pop1, self.pop2]

        # project each population down by one diploid sample (i.e. 2) to maximise the data in the frequency spectrum
        self.prjt = [(len(POPULATIONS[pop]) - 1) * 2 for pop in self.pops]

    def requires(self):
        return SiteFrequencySpectrum(self.group, GENOME)

    def output(self):
        return luigi.LocalTarget("fsdata/{0}_{1}_{2}_{3}.fs".format(self.pops[0], self.prjt[0],
                                                                    self.pops[1], self.prjt[1]))

    def run(self):

        # parse the data file to generate the data dictionary
        dd = dadi.Misc.make_data_dict(self.input().path)

        # extract the spectrum for the two populations from the dictionary and project down
        fs = dadi.Spectrum.from_data_dict(dd, self.pops, self.prjt, polarized=True)

        # save it to a file
        fs.to_file(self.output().path)


class DadiOptimizeLogParams(SetupTask):
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
        return [luigi.LocalTarget("fsdata/{0}.{1}.opt".format("split_mig", self.n)),
                luigi.LocalTarget("pdf/{0}_{1}_{2}.dadi.fs.pdf".format(self.pop1, self.pop2, self.n))]

    def run(self):

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file("fsdata/{0}_{1}.fs".format(self.pop1, self.pop2))

        # TODO what effect do these have
        # These are the grid point settings will use for extrapolation.
        pts_l = [10, 50, 60]

        # The Demographics1D and Demographics2D modules contain a few simple models,
        # mostly as examples. We could use one of those.
        func = dadi.Demographics2D.split_mig

        # The upper_bound and lower_bound lists are for use in optimization.
        upper_bound = [100, 100, 3, 10]
        lower_bound = [1e-4, 1e-4, 0, 0]

        # randomly generated starting values within the bounding ranges
        p0 = [random.uniform(lower_bound[i], upper_bound[i]) for i in range(0, len(upper_bound))]

        # This is our initial guess for the parameters, which is somewhat arbitrary.
        # p0 = [2, 0.1, 0.2, 0.2]

        # p0 = [0.0006,  # nu1: Size of population 1 after split.
        #       0.1,     # nu2: Size of population 2 after split.
        #       0.002,   # T: Time in the past of split (in units of 2*Na generations)
        #       0.001]   # m: Migration rate between populations (2*Na*m)

        # Perturb our parameters before optimization. This does so by taking each
        # parameter a up to a factor of two up or down.
        p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

        # Make the extrapolating version of our demographic model function.
        func_ex = dadi.Numerics.make_extrap_log_func(func)

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
        with self.output()[0].open('w') as fout:
            fout.write('Starting params: {0}\n'.format(p0))
            fout.write('Maximum log composite likelihood: {0}\n'.format(ll_model))
            fout.write('Best-fit parameters: {0}\n'.format(popt))
            fout.write('Optimal value of theta: {0}\n'.format(theta))

        # # save the figure as a PDF
        # fig = plt.figure(1)
        # dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=1, resid_range=3, fig_num=1)
        # fig.savefig(self.output()[1].path)
        # plt.close(fig)

        # # estimate of mutation rate... probably totally wrong
        # mu=1.25e-8
        #
        # # theta=4*Ne*mu
        # # Ne=theta/(4*mu)
        # # reference effective population size of the ancestral population (kind of)
        # Ne = theta / (4 * mu)
        #
        # # Ne / theta for population 1
        # theta1=popt[1]*theta
        # N1=popt[1]*Ne
        #
        # # Ne / theta for population 2
        # theta2=popt[2]*theta
        # N2 = popt[2] * Ne
        #
        # #T=2*Ne*t
        # #t=T/(2*Ne)
        # # This the time in generations
        # time=popt[3]*(2*Ne)
        #
        # print "mu={}".format(mu)
        # print "Ne={}".format(Ne)
        # print "theta1={}".format(theta1)
        # print "N1={}".format(N1)
        # print "theta2={}".format(theta2)
        # print "N2={}".format(N2)
        # print "time={}".format(time)


class CustomDadiPipeline(luigi.WrapperTask):
    """
    Run the dadi models
    """

    def requires(self):

        for n in range(0, 100000):
            yield DadiOptimizeLogParams('DOM', 'WLD-FRE', n)
