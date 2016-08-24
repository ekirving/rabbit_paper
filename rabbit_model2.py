# Numpy is the numerical library dadi is built upon
from numpy import array

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import dadi

# Load the data
fs = dadi.Spectrum.from_file('fsdata/all-pops_DOM_WLD-FRE.fs')
ns = fs.sample_sizes

# These are the grid point settings will use for extrapolation.
# pts_l = [40,50,60]

pts_l = [ns[0]+10, ns[0]+20, ns[0]+30]

# Isolation-with-migration model with exponential pop growth.
func = dadi.Demographics2D.IM

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.

# crazy big bounds !!!
upper_bound = [1-1e-10,  100, 100, 10, 1e5, 1e5]
lower_bound = [  1e-10, 1e-10, 1e-10,  0,   0,   0]

# more realistic
# upper_bound = [0.9999, 100,  100, 10, 30, 30]
# lower_bound = [0.001, 1e-2, 1e-2, 0,  0,  0]


# upperIM = [0.99, 10, 10, 5, 20, 20]
# lowerIM = [0.01, 0.5, 0.5, 0.005, .1, .1]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [0.5,  # s: Size of pop 1 after split. (Pop 2 has size 1-s.)
      0.1,  # nu1: Final size of pop 1.
      0.1,  # nu2: Final size of pop 2.
      0.2,  # T: Time in the past of split (in units of 2*Na generations)
      0.2,  # m12: Migration from pop 2 to pop 1 (2*Na*m12)
      0.2]  # m21: Migration from pop 1 to pop 2

# p0 = [ 0.0537521  ,  0.0100264  ,  0.0684193  ,  0.0068241  ,  0.0434472  ,  0.0262704  ] # -2257.8
# p0 = [ 0.0491897  ,  0.0100047  ,  0.130444   ,  0.00764804 ,  0.0424375  ,  0.0175239  ] # -1986.96
# p0 = [ 0.0367591  ,  0.0100035  ,  0.168718   ,  0.00690335 ,  0.0250734  ,  0.0289296  ] # -1924.57
# p0 = [ 0.0303529  ,  0.0106935  ,  2.60008    ,  0.00660695 ,  1.93043e-06,  0.000138871] # -1866.61
# p0 = [ 0.0141832  ,  0.0100017  ,  0.482211   ,  0.00466663 ,  6.20664e-07,  0.000143181] # -1835.16
# p0 = [ 0.00865567 ,  0.0100043  ,  0.264837   ,  0.00364864 ,  1.23579e-06,  9.23756e-05] # -1821.73

# best fit from manual runs
# p0 = [ 0.00365723 ,  0.0100486  ,  0.385703   ,  0.00226085 ,  1.29317e-06,  5.33581e-05] # -1796.36

# p0 = [ 0.002396513,	0.010378533,	0.295010157,	0.001797582,	9.59E-07,	6.51E-05] # -1781.545499

# best fit from opt logs... really weird numbers >:(
# p0 = [ 0.021885814,	0.046397118,	0.682089788,	0.017228833,	9.696112236,	9.977892015]

# p0 = [ 0.022131   ,  0.0377687  ,  0.478812   ,  0.0169264  ,  9.78311    ,  19.8904]
# p0 = [ 0.0191467  ,  0.01       ,  0.0638444  ,  0.00983257 ,  37.1688    ,  66.0204]
# p0 = [ 0.00131067 ,  0.000101352,  0.000898061,  0.000221085,  3756.02    ,  4290.99    ]
# p0 = [ 0.000468924,  4.54758e-05,  0.000377352,  9.98896e-05,  9144.36    ,  9999.13    ]
p1 = [ 4.53703467e-05,   4.71541615e-06,   2.95542838e-05,   1.01754795e-05, 8.90730854e+04,   9.96531052e+04]
# p0 = [ 0.002 ,  0.1,  0.1,  0.000221085,  20    ,  25]

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# # Make the extrapolating version of our demographic model function.
# func_ex = dadi.Numerics.make_extrap_log_func(func)
#
# # Calculate the best-fit model AFS.
# model = func_ex(p0, ns, pts_l)
#
# # Likelihood of the data given the model AFS.
# ll_model = dadi.Inference.ll_multinom(model, fs)
#
# # The optimal value of theta given the model.
# theta = dadi.Inference.optimal_sfs_scaling(model, fs)
#
# print ('P best: {0}'.format(p0))
# print('Maximum log composite likelihood: {0}'.format(ll_model))
# print('Optimal value of theta: {0}'.format(theta))
#
# quit()
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

verbose = True

# make a list of the optimal params
popt = []

for i in range(0, 1):

    # Make the extrapolating version of our demographic model function.
    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # Perturb our parameters before optimization. This does so by taking each
    # parameter a up to a factor of two up or down.
    # p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
    #                               lower_bound=lower_bound)

    if verbose : print('Beginning optimization ************************************************')

    # # do the optimization...
    # p1 = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
    #                                  lower_bound=lower_bound,
    #                                  upper_bound=upper_bound,
    #                                  verbose=20 if verbose else 0,
    #                                  maxiter=100)

    if verbose: print('Finshed optimization **************************************************')

    # Calculate the best-fit model AFS.
    model = func_ex(p1, ns, pts_l)

    # Likelihood of the data given the model AFS.
    ll_model = dadi.Inference.ll_multinom(model, fs)

    # The optimal value of theta given the model.
    theta = dadi.Inference.optimal_sfs_scaling(model, fs)

    if verbose:
        print ('P best: {0}'.format(p1))
        print('Maximum log composite likelihood: {0}'.format(ll_model))
        print('Optimal value of theta: {0}'.format(theta))

    # record the best fitting params for this run
    popt.append([ll_model, theta] + list(p1))

    # sort the params (largest ll_model first)
    popt.sort(reverse=True)

    # run the optimisation again, starting from the best fit we've seen so far
    p0 = popt[0][2:]

    # plot the figure
    fig = plt.figure(1)
    dadi.Plotting.plot_2d_comp_multinom(model, fs, fig_num=1)
    fig.savefig('test.pdf')
    plt.close(fig)

# display the optimal params
print popt



    # # These are the actual best-fit model parameters, which we found through
    # # longer optimizations and confirmed by running multiple optimizations.
    # # We'll work with them through the rest of this script.
    # popt = [1.881, 0.0710, 1.845, 0.911, 0.355, 0.111]
    # print('Best-fit parameters: {0}'.format(popt))
    #
    # # Calculate the best-fit model AFS.
    # model = func_ex(popt, ns, pts_l)
    # # Likelihood of the data given the model AFS.
    # ll_model = dadi.Inference.ll_multinom(model, data)
    # print('Maximum log composite likelihood: {0}'.format(ll_model))
    # # The optimal value of theta given the model.
    # theta = dadi.Inference.optimal_sfs_scaling(model, data)
    # print('Optimal value of theta: {0}'.format(theta))
    #
    # # Plot a comparison of the resulting fs with the data.
    # import pylab
    # pylab.figure(1)
    # dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
    #                                     pop_ids =('YRI','CEU'))
    # # This ensures that the figure pops up. It may be unecessary if you are using
    # # ipython.
    # pylab.show()
    # # Save the figure
    # pylab.savefig('YRI_CEU.png', dpi=50)
    #
    # # Let's generate some data using ms, if you have it installed.
    # mscore = demographic_models.prior_onegrow_mig_mscore(popt)
    # # I find that it's most efficient to simulate with theta=1, average over many
    # # iterations, and then scale up.
    # mscommand = dadi.Misc.ms_command(1., ns, mscore, int(1e5))
    # # If you have ms installed, uncomment these lines to see the results.
    #
    # # We use Python's os module to call this command from within the script.
    # import os
    # return_code = os.system('{0} > test.msout'.format(mscommand))
    # # We check the return code, so the script doesn't crash if you don't have ms
    # # installed
    # if return_code == 0:
    #     msdata = dadi.Spectrum.from_ms_file('test.msout')
    #     pylab.figure(2)
    #     dadi.Plotting.plot_2d_comp_multinom(model, theta*msdata, vmin=1,
    #                                         pop_ids=('YRI','CEU'))
    #     pylab.show()
    #
    # # Estimate parameter uncertainties using the Godambe Information Matrix, to
    # # account for linkage in the data. To use the GIM approach, we need to have
    # # spectra from bootstrapping our data.  Let's load the ones we've provided for
    # # the example.
    # # (We're using Python list comprehension syntax to do this in one line.)
    # all_boot = [dadi.Spectrum.from_file('bootstraps/{0:02d}.fs'.format(ii))
    #             for ii in range(100)]
    # uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data,
    #                                   multinom=True)
    # # uncert contains the estimated standard deviations of each parameter, with
    # # theta as the final entry in the list.
    # print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))
    #
    # # For comparison, we can estimate uncertainties with the Fisher Information
    # # Matrix, which doesn't account for linkage in the data and thus underestimates
    # # uncertainty. (Although it's a fine approach if you think your data is truly
    # # unlinked.)
    # uncerts_fim = dadi.Godambe.FIM_uncert(func_ex, pts_l, popt, data, multinom=True)
    # print('Estimated parameter standard deviations from FIM: {0}'.format(uncerts_fim))
    #
    # print('Factors by which FIM underestimates parameter uncertainties: {0}'.format(uncerts/uncerts_fim))
    #
    # # What if we fold the data?
    # # These are the optimal parameters when the spectrum is folded. They can be
    # # found simply by passing data.fold() to the above call to optimize_log.
    # popt_fold =  array([1.907,  0.073,  1.830,  0.899,  0.425,  0.113])
    # uncerts_folded = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt_fold,
    #                                          data.fold(), multinom=True)
    # print('Folding increases parameter uncertainties by factors of: {0}'.format(uncerts_folded/uncerts))
    #
    # # Let's do a likelihood-ratio test comparing models with and without migration.
    # # The no migration model is implemented as
    # # demographic_models.prior_onegrow_nomig
    # func_nomig = demographic_models.prior_onegrow_nomig
    # func_ex_nomig = dadi.Numerics.make_extrap_log_func(func_nomig)
    # # These are the best-fit parameters, which we found by multiple optimizations
    # popt_nomig = array([ 1.897,  0.0388,  9.677,  0.395,  0.070])
    # model_nomig = func_ex_nomig(popt_nomig, ns, pts_l)
    # ll_nomig = dadi.Inference.ll_multinom(model_nomig, data)
    #
    # # Since LRT evaluates the complex model using the best-fit parameters from the
    # # simple model, we need to create list of parameters for the complex model
    # # using the simple (no-mig) best-fit params.  Since evalution is done with more
    # # complex model, need to insert zero migration value at corresponding migration
    # # parameter index in complex model. And we need to tell the LRT adjust function
    # # that the 3rd parameter (counting from 0) is the nested one.
    # p_lrt = [1.897,  0.0388,  9.677, 0, 0.395,  0.070]
    #
    # adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, p_lrt, data,
    #                               nested_indices=[3], multinom=True)
    # D_adj = adj*2*(ll_model - ll_nomig)
    # print('Adjusted D statistic: {0:.4f}'.format(D_adj))
    #
    # # Because this is test of a parameter on the boundary of parameter space
    # # (m cannot be less than zero), our null distribution is an even proportion
    # # of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
    # # point percent function for a weighted sum of chi^2 dists.
    # pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
    # print('p-value for rejecting no-migration model: {0:.4f}'.format(pval))
