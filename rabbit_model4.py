# Numpy is the numerical library dadi is built upon
from numpy import array

import dadi
import pylab
import os

# Load the data
fs = dadi.Spectrum.from_file('fsdata/all-pops_DOM_WLD-FRE.fs')
ns = fs.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [40,50,60]

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# Isolation-with-migration model with exponential pop growth.
func_complex = dadi.Demographics2D.IM

# best fit params from multiple optimisation runs
popt_complex = [0.021885814, 0.046397118, 0.682089788, 0.017228833, 9.696112236, 9.977892015]



# Split into two populations of specifed size, with migration.
func_simple = dadi.Demographics2D.split_mig

# best fit params from multiple optimisation runs
popt_simple = [ 0.040546754, 0.522623312, 0.024983701, 9.98911627]

# Since LRT evaluates the complex model using the best-fit parameters from the
# simple model, we need to create list of parameters for the complex model
# using the simple (no-mig) best-fit params.  Since evalution is done with more
# complex model, need to insert zero migration value at corresponding migration
# parameter index in complex model. And we need to tell the LRT adjust function
# that the 3rd parameter (counting from 0) is the nested one.

p_lrt = [ 0.5, popt_simple[1], popt_simple[2], popt_simple[3], popt_simple[4], popt_simple[4]]

nested_indices=[0, 5]

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func_complex)

# Calculate the best-fit model AFS.
model = func_ex(popt_complex, ns, pts_l)

# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)

# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

print('Maximum log composite likelihood: {0}'.format(ll_model))
print('Optimal value of theta: {0}'.format(theta))

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# Let's generate some data using ms
mscore = dadi.Demographics2D.IM_mscore(popt_complex)
mscommand = dadi.Misc.ms_command(1., ns, mscore, int(1e5))

# We use Python's os module to call this command from within the script.
return_code = os.system('{0} > test.msout'.format(mscommand))

all_boot = dadi.Spectrum.from_ms_file('test.msout', bootstrap_segments=100)

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt_complex, fs,
                                  multinom=True)

# uncert contains the estimated standard deviations of each parameter, with
# theta as the final entry in the list.
print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))

# For comparison, we can estimate uncertainties with the Fisher Information
# Matrix, which doesn't account for linkage in the data and thus underestimates
# uncertainty. (Although it's a fine approach if you think your data is truly
# unlinked.)
uncerts_fim = dadi.Godambe.FIM_uncert(func_ex, pts_l, popt_complex, fs, multinom=True)

print('Estimated parameter standard deviations from FIM: {0}'.format(uncerts_fim))
print('Factors by which FIM underestimates parameter uncertainties: {0}'.format(uncerts/uncerts_fim))

    # TODO look at folding...
    # # # What if we fold the data?
    # # # These are the optimal parameters when the spectrum is folded. They can be
    # # # found simply by passing data.fold() to the above call to optimize_log.
    # # popt_fold =  array([1.907,  0.073,  1.830,  0.899,  0.425,  0.113])
    # uncerts_folded = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt_fold,
    #                                          fs.fold(), multinom=True)
    # # print('Folding increases parameter uncertainties by factors of: {0}'.format(uncerts_folded/uncerts))

# Let's do a likelihood-ratio test comparing the complex and simple models
func_ex_simple = dadi.Numerics.make_extrap_log_func(func_simple)
model_simple = func_ex_simple(popt_simple, ns, pts_l)
ll_simple = dadi.Inference.ll_multinom(model_simple, fs)

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, p_lrt, fs, nested_indices=nested_indices, multinom=True)

D_adj = adj*2*(ll_model - ll_simple)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))

# Because this is test of a parameter on the boundary of parameter space
# (m cannot be less than zero), our null distribution is an even proportion
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))

print('p-value for rejecting simple model: {0:.4f}'.format(pval))