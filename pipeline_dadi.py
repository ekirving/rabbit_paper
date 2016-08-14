#!/usr/bin/env python
# -*- coding: utf-8 -*-

# load matplotlib before dadi so we can disable the screen
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import luigi, dadi, numpy, pylab, random
import pickle
import csv

# import the custom pipelines
from pipeline_gatk import *
from pipeline_utils import *

# buffer the logs created by dadi so we can inspect them
log_buffer = LogBuffer()
logger = logging.getLogger('Inference')
logger.addHandler(logging.StreamHandler(log_buffer))

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
        return luigi.LocalTarget("fsdata/{0}_{1}_{2}.fs".format(self.group, self.pop1, self.pop2))

    def run(self):

        # parse the data file to generate the data dictionary
        dd = dadi.Misc.make_data_dict(self.input().path)

        # get the two populations
        pops = [self.pop1, self.pop2]

        # project each population down by one sample to allow for a little missing coverage
        prj = [(len(GROUPS[self.group][pop]) - 1) * 2 for pop in pops]

        # extract the spectrum for the two populations from the dictionary and project down
        fs = dadi.Spectrum.from_data_dict(dd, pops, prj, polarized=True)

        # save it to a file
        fs.to_file(self.output().path)


class DadiModelOptimizeParams(luigi.Task):
    """
    Optimise the log likelihood of the model paramaters for the given frequency spectrum
    """
    group = luigi.Parameter()
    pop1 = luigi.Parameter()
    pop2 = luigi.Parameter()

    model = luigi.Parameter()
    scenario = luigi.Parameter()
    grid_size = luigi.ListParameter()
    upper_bound = luigi.ListParameter()
    lower_bound = luigi.ListParameter()
    param_start = luigi.ListParameter(significant=False)

    n = luigi.IntParameter()

    def requires(self):
        return DadiSpectrum(self.group, self.pop1, self.pop2)

    def output(self):
        return luigi.LocalTarget("fsdata/opt/{0}_{1}_{2}_{3}_{4}_{5}.opt".format(self.group, self.pop1, self.pop2,
                                                                                 self.model, self.scenario, self.n))

    def run(self):

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file(self.input().path)
        ns = fs.sample_sizes

        # get the demographic model to test
        func = getattr(dadi.Demographics2D, self.model)

        # Make the extrapolating version of our demographic model function.
        func_ex = dadi.Numerics.make_extrap_log_func(func)

        # keep a list of the optimal params
        p_best = []

        # start with the random values passed to this task
        p_start = self.param_start

        # run the optimisation many times
        for i in range(0, DADI_MAX_ITER):

            # Perturb our parameters before optimization. This does so by taking each
            # parameter a up to a factor of two up or down.
            p_perturb = dadi.Misc.perturb_params(p_start,
                                                 fold=1,
                                                 upper_bound=self.upper_bound,
                                                 lower_bound=self.lower_bound)

            print('Beginning optimization: {} | {} | {} **** n={:>3} **** i={:>3} **** '.format(self.group,
                                                                                                self.model,
                                                                                                self.scenario,
                                                                                                self.n,
                                                                                                i+1))
            start = datetime.datetime.now()

            # do the optimization...
            p_opt = dadi.Inference.optimize_log(p_perturb, fs, func_ex, self.grid_size,
                                                lower_bound=self.lower_bound,
                                                upper_bound=self.upper_bound,
                                                verbose=20,
                                                maxiter=DADI_MAX_ITER)

            end = datetime.datetime.now()
            diff = int((end - start).total_seconds() / 60)

            print('Finshed optimization: {} | {} | {} **** n={:>3} **** i={:>3} **** t={:>3} mins'.format(self.group,
                                                                                                          self.model,
                                                                                                          self.scenario,
                                                                                                          self.n,
                                                                                                          i+1,
                                                                                                          diff))

            # reset the log buffer
            log_buffer.log = []

            # Calculate the best-fit model AFS.
            model = func_ex(p_opt, ns, self.grid_size)

            # Likelihood of the data given the model AFS.
            ll_model = dadi.Inference.ll_multinom(model, fs)

            # The optimal value of theta given the model.
            theta = dadi.Inference.optimal_sfs_scaling(model, fs)

            # get the buffered warnings
            warnings = " ".join(log_buffer.log)

            # we only care about non-masked data
            if "Model is masked" not in warnings:

                print('Maximum log composite likelihood: {0}'.format(ll_model))
                print('Optimal value of theta: {0}'.format(theta))

                # record the best fitting params for this run
                p_best.append([ll_model, theta] + list(p_opt))

                # sort the params (largest ll_model first)
                p_best.sort(reverse=True)

            try:
                # run the optimisation again, starting from the best fit we've seen so far
                p_start = p_best[0][2:]
            except IndexError:
                # otherwise, generate a set of new random starting params (so we don't get stuck in bad param space)
                p_start = [random.uniform(self.lower_bound[j], self.upper_bound[j])
                           for j in range(0, len(self.upper_bound))]

        # if we've run the iteration DADI_MAX_ITER times and not found any non-masked params then we've failed
        if not p_best:
            raise Exception("{} | {} | {} | n={}: FAILED to find any non-masked params".format(self.group,
                                                                                               self.model,
                                                                                               self.scenario,
                                                                                               self.n))
        # save the list of optimal params by pickeling them in a file
        with self.output().open('w') as fout:
            pickle.dump(p_best, fout)


class DadiModelMaximumLikelihood(luigi.Task):
    """
    Find the maximum log likelihood parameters for a given model, and plot the resulting spectrum
    """
    group = luigi.Parameter()
    pop1 = luigi.Parameter()
    pop2 = luigi.Parameter()

    model = luigi.Parameter()
    scenario = luigi.Parameter()
    param_names = luigi.ListParameter()
    grid_size = luigi.ListParameter()
    upper_bound = luigi.ListParameter()
    lower_bound = luigi.ListParameter()

    def requires(self):
        # TODO restore DADI_MAX_ITER when done testing
        for n in range(0, 10):
            # randomly generate starting params, within the bounding ranges
            param_start = [random.uniform(self.lower_bound[i], self.upper_bound[i])
                           for i in range(0, len(self.upper_bound))]

            # find the optimal params
            yield DadiModelOptimizeParams(self.group, self.pop1, self.pop2, self.model, self.scenario, self.grid_size,
                                          self.upper_bound, self.lower_bound, param_start, n)

    def output(self):
        return [luigi.LocalTarget("fsdata/{0}_{1}_{2}_{3}_{4}.csv".format(self.group, self.pop1, self.pop2, self.model,
                                                                          self.scenario)),
                luigi.LocalTarget("fsdata/{0}_{1}_{2}_{3}_{4}.pdf".format(self.group, self.pop1, self.pop2, self.model,
                                                                          self.scenario))]

    def run(self):

        # make a master list of optimal params
        p_best = []

        for opt in self.input():
            with opt.open('r') as fin:
                # unpickle the optimal params
                p_best += pickle.load(fin)

        # sort the params (largest ll_model first)
        p_best.sort(reverse=True)

        header = ["likelihood", "theta"] + list(self.param_names)

        # dump all the data to csv
        with self.output()[0].open("w") as fout:
            writer = csv.writer(fout)
            writer.writerow(header)
            writer.writerows(p_best)

        # get the params with the maximum log likelihood, from all ~1e6 iterations
        p_opt = p_best[0][2:]

        # load the frequency spectrum
        fs = dadi.Spectrum.from_file("fsdata/{0}_{1}_{2}.fs".format(self.group, self.pop1, self.pop2))
        ns = fs.sample_sizes

        # get the demographic model to test
        func = getattr(dadi.Demographics2D, self.model)

        # Make the extrapolating version of our demographic model function.
        func_ex = dadi.Numerics.make_extrap_log_func(func)

        # Calculate the best-fit model AFS.
        model = func_ex(p_opt, ns, self.grid_size)

        # plot the figure
        fig = plt.figure(1)
        dadi.Plotting.plot_2d_comp_multinom(model, fs, fig_num=1)
        fig.savefig(self.output()[1].path)
        plt.close(fig)


class CustomDadiPipeline(luigi.WrapperTask):
    """
    Run the dadi models
    """

    def requires(self):

        # run the whole analysis for multiple sets of pairwise comparisons
        for group, pop1, pop2 in [('all-pops',  'DOM',     'WLD-FRE'),
                                  ('all-pops',  'WLD-FRE', 'WLD-IB2'),
                                  ('all-pops',  'WLD-IB2', 'WLD-IB1'),
                                  ('no-admix',  'DOM',     'WLD-FRE2'),
                                  ('no-admix',  'DOM',     'WLD-FRE3'),
                                  ('split-fre', 'FRE-1',   'FRE-2')]:

            # TODO experiment to see how changing this effects run time and model fitting
            grid_size = [10, 50, 60]

            # ----------------------------------------------------------------------------------------------------------

            model = 'split_mig'
            param_names = ['nu1',  # Size of population 1 after split.
                           'nu2',  # Size of population 2 after split.
                           'T',    # Time in the past of split (in units of 2*Na generations)
                           'm']    # Migration rate between populations (2*Na*m)

            scenario = "best-fit"
            upper_bound = [100, 100, 3, 10]
            lower_bound = [1e-2, 1e-2, 0, 0]

            yield DadiModelMaximumLikelihood(group, pop1, pop2, model, scenario, param_names, grid_size,
                                             upper_bound, lower_bound)

            # ----------------------------------------------------------------------------------------------------------

            model = 'IM'
            param_names = ['s',    # Size of pop 1 after split. (Pop 2 has size 1-s.)
                           'nu1',  # Final size of pop 1.
                           'nu2',  # Final size of pop 2.
                           'T',    # Time in the past of split (in units of 2*Na generations)
                           'm12',  # Migration from pop 2 to pop 1 (2*Na*m12)
                           'm21']  # Migration from pop 1 to pop 2

            # -------------------
            scenario = "best-fit"
            upper_bound = [0.9999, 100, 100, 3, 10, 10]
            lower_bound = [0.0001, 1e-2, 1e-2, 0, 0, 0]

            yield DadiModelMaximumLikelihood(group, pop1, pop2, model, scenario, param_names, grid_size,
                                             upper_bound, lower_bound)

            # -------------------
            # now lets run the model where we fix certain params, so we can do a direct comparison between them...
            # simulate a simple model, by fixing S
            scenario = "fixed-S"
            upper_bound = [0.5, 100, 100, 3, 10, 10]
            lower_bound = [0.5, 1e-2, 1e-2, 0, 0, 0]

            yield DadiModelMaximumLikelihood(group, pop1, pop2, model, scenario, param_names, grid_size,
                                             upper_bound, lower_bound)

            # -------------------
            # simulate a simple model, by fixing S, m12, m21
            scenario = "fixed-S-m12-m21"
            upper_bound = [0.5, 100, 100, 3, 0, 0]
            lower_bound = [0.5, 1e-2, 1e-2, 0, 0, 0]

            yield DadiModelMaximumLikelihood(group, pop1, pop2, model, scenario, param_names, grid_size,
                                             upper_bound, lower_bound)

            # -------------------
            # simulate a simple model, by fixing m12, m21
            scenario = "fixed-m12-m21"
            upper_bound = [0.9999, 100, 100, 3, 0, 0]
            lower_bound = [0.0001, 1e-2, 1e-2, 0, 0, 0]

            yield DadiModelMaximumLikelihood(group, pop1, pop2, model, scenario, param_names, grid_size,
                                             upper_bound, lower_bound)

