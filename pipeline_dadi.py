#!/usr/bin/env python
# -*- coding: utf-8 -*-

# load matplotlib before dadi so we can disable the screen
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import luigi, dadi, numpy, pylab, random

# import the custom pipelines
from pipeline_gatk import *
from pipeline_utils import *

# buffer the logs created by dadi so we can inspect them
strm = LogBuffer()
hndl = logging.StreamHandler(strm)
logger = logging.getLogger('Inference')
logger.addHandler(hndl)

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

        # project each population down by one sample to allow for a little missing coverage
        prj = [(len(POPULATIONS[pop]) - 1)* 2 for pop in pops]

        # extract the spectrum for the two populations from the dictionary and project down
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
        return luigi.LocalTarget("pdf/dadi.{0}_{1}_{2}.fs.pdf".format(self.pop1, self.pop2, self.n))

    def run(self):

        # reset the buffered log (just in case)
        strm.log = []

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file("fsdata/{0}_{1}.fs".format(self.pop1, self.pop2))

        # TODO what effect do these have
        # These are the grid point settings will use for extrapolation.
        pts_l = [10, 50, 60]

        # TODO make model a param
        func = dadi.Demographics2D.split_mig

        model_params = ["nu1",  # nu1: Size of population 1 after split.
                        "nu2",  # nu2: Size of population 2 after split.
                        "T",    # T: Time in the past of split (in units of 2*Na generations)
                        "m"]    # m: Migration rate between populations (2*Na*m)

        # The upper_bound and lower_bound lists are for use in optimization.
        upper_bound = [100, 100, 3, 10]
        lower_bound = [1e-4, 1e-4, 0, 0]

        # randomly generated starting values within the bounding ranges
        p0 = [random.uniform(lower_bound[i], upper_bound[i]) for i in range(0, len(model_params))]

        # This is our initial guess for the parameters, which is somewhat arbitrary.
        # p0 = [2, 0.1, 0.2, 0.2]

        # Perturb our parameters before optimization. This does so by taking each
        # parameter a up to a factor of two up or down.
        # p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

        # Make the extrapolating version of our demographic model function.
        func_ex = dadi.Numerics.make_extrap_log_func(func)

        # Do the optimization...
        popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                           lower_bound=lower_bound,
                                           upper_bound=upper_bound,
                                           maxiter=DADI_MAX_ITER)

        # Calculate the best-fit model AFS.
        model = func_ex(popt, fs.sample_sizes, pts_l)

        # Likelihood of the data given the model AFS.
        ll_model = dadi.Inference.ll_multinom(model, fs)

        # The optimal value of theta given the model.
        theta = dadi.Inference.optimal_sfs_scaling(model, fs)

        # collate all the data for logging
        data = [self.pop1, self.pop2, self.n] + popt + [ll_model, theta] + strm.log

        # add an entry to to log
        with open("fsdata/{0}_{1}_{2}.tsv".format("split_mig", self.pop1, self.pop2), "a") as tsv:
            tsv.write("\t".join(str(x).strip("\n") for x in data) + "\n")

        # # save the figure as a PDF
        # fig = plt.figure(1)
        # dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin=1, resid_range=3, fig_num=1)
        # fig.savefig(self.output()[1].path)
        # plt.close(fig)

        # # estimate of mutation rate... probably totally wrong
        # mu = 1.25e-8
        # mu = 1e-3         # rate for microsatelites in rabbits (Surridge et al., 1999)
        # mu = 1.622e-9     # rate for rabbits based on assumed Oryctolagus-Lepus split time of 11.8 MYA (Carneiro et al 2009:596)
        # mu = 1.74e-9      # (Carneiro et al., 2011)
        # mu = 2.02-2.35e-9 # mutations per site per generation, which is similar to but slightly lower than estimates in mice (3.4-4.1e-9; (Carneiro et al., 2012)
        # mu = 1.18e-6      # (Sousa et al., 2013)  μ

        # # theta=4*Ne*mu
        # # Ne=theta/(4*mu)
        # # reference effective population size of the ancestral population (kind of)
        # Ne = theta / (4 * mu)
        #
        # # Ne / theta for population 1
        # theta1=popt[0]*theta
        # N1=popt[0]*Ne
        #
        # # Ne / theta for population 2
        # theta2=popt[1]*theta
        # N2 = popt[1] * Ne
        #
        # #T=2*Ne*t
        # #t=T/(2*Ne)
        # # This the time in generations
        # time=popt[2]*(2*Ne)
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

        for n in range(0, 10000):
            yield DadiOptimizeLogParams('DOM', 'WLD-FRE', n)
