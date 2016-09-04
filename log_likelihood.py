# Numpy is the numerical library dadi is built upon
from numpy import array

import dadi
import pylab
import os


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# Isolation-with-migration model with exponential pop growth.
func_complex = dadi.Demographics2D.IM

# s: Size of pop 1 after split. (Pop 2 has size 1-s.)
# nu1: Final size of pop 1.
# nu2: Final size of pop 2.
# T: Time in the past of split (in units of 2*Na generations)
# m12: Migration from pop 2 to pop 1 (2*Na*m12)
# m21: Migration from pop 1 to pop 2

# best fit params from multiple optimisation runs

# DOM /FRE
# fs = dadi.Spectrum.from_file('fsdata/all-pops_DOM_WLD-FRE.fs')
# popt_complex = [0.021885814, 0.046397118, 0.682089788, 0.017228833, 9.696112236, 9.977892015]
# popt_simple = [ 0.040546754, 0.522623312, 0.024983701, 9.98911627]

    # Adjusted D statistic: -0.1258
    # p-value for rejecting simple model: 1.0000

# DOM / FRE2
# fs = dadi.Spectrum.from_file('fsdata/no-admix_DOM_WLD-FRE2.fs')
# popt_complex = [0.973950359, 0.036676267, 1.801709748, 0.119370203, 9.735916425, 9.961141793]
# popt_simple = [0.051099955, 0.308993758, 0.045069965, 9.976460257]

    # Adjusted D statistic: 0.0577
    # p-value for rejecting simple model: 0.4051

    # BAD!!!! -37327.96607
    # popt_simple= [60.13287855, 82.47864547, 0.179351516, 0.416954267]
    # popt_simple= [0.025804311, 0.472274087, 0.012950424, 8.25229803]

# DOM / FRE3
# fs = dadi.Spectrum.from_file('fsdata/no-admix_DOM_WLD-FRE3.fs')
# popt_complex = [0.905197319, 0.04151391, 0.325184361, 0.141967601, 9.900435858, 9.458934399]
# popt_simple = [0.058741985, 0.178016717, 0.059008035, 9.973357793]

    # Adjusted D statistic: 0.0132
    # p-value for rejecting simple model: 0.4543

# Split FRE1/2
# fs = dadi.Spectrum.from_file('fsdata/split-fre_FRE-1_FRE-2.fs')
# popt_complex = [0.013048138, 11.73178929, 1.918929748, 0.025678896, 7.545295451, 9.979874929]
# popt_simple = [0.046784111, 9.674163542, 0.008100691, 9.976127942]

    # Adjusted D statistic: 0.0122
    # p-value for rejecting simple model: 0.4561

# FRE / IB2
# fs = dadi.Spectrum.from_file('fsdata/all-pops_WLD-FRE_WLD-IB2.fs')
# popt_complex = [0.008112191, 0.256233251, 1.124115328, 0.638081552, 8.11299697, 2.115227161]
# popt_simple = [0.28917809, 1.269741392, 0.10457709, 2.836067032]

    # Adjusted D statistic: 0.0283
    # p-value for rejecting simple model: 0.4333

# IB2 / IB1
fs = dadi.Spectrum.from_file('fsdata/all-pops_WLD-IB2_WLD-IB1.fs')
popt_complex = [0.984309077, 0.275898384, 2.691461934, 0.834728827, 3.708014791, 2.07175748]
popt_simple = [0.719422317, 2.568937764, 0.477190374, 1.192350901]

    # Adjusted D statistic: 0.0118
    # p-value for rejecting simple model: 0.4567

# Split into two populations of specifed size, with migration.
func_simple = dadi.Demographics2D.split_mig

# nu1: Size of population 1 after split.
# nu2: Size of population 2 after split.
# T: Time in the past of split (in units of 2*Na generations)
# m: Migration rate between populations (2*Na*m)

# best fit params from multiple optimisation runs


# make 's' the final pop-size ratio
s = popt_simple[0]/popt_simple[1]

# and the second migration param the same as the first
m2 = popt_simple[3]

# now compose the param set for evaluating the simple model
# p_lrt = [s] + popt_simple + [m2]
p_lrt = [s] + popt_simple + [m2]

# tell the evaluation that the first and last params are 'nested'
nested_indices=[0, 5]

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
ns = fs.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [40,50,60]

# ----------------------------------------------------------------------------------------------------------------------

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func_complex)

# Calculate the best-fit model AFS.
model = func_ex(popt_complex, ns, pts_l)

# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)

# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

print('\n--- complex ---');
print('Maximum log composite likelihood: {0}'.format(ll_model))
print('Optimal value of theta: {0}'.format(theta))

# # Calculate the best-fit model AFS.
# model2 = func_ex(p_lrt, ns, pts_l)
#
# # Likelihood of the data given the model AFS.
# ll_model2 = dadi.Inference.ll_multinom(model2, fs)
#
# # The optimal value of theta given the model.
# theta2 = dadi.Inference.optimal_sfs_scaling(model2, fs)
#
# print('\n--- nested ---');
# print('Maximum log composite likelihood: {0}'.format(ll_model2))
# print('Optimal value of theta: {0}'.format(theta2))


# ----------------------------------------------------------------------------------------------------------------------

# Make the extrapolating version of our demographic model function.
func_ex_simple = dadi.Numerics.make_extrap_log_func(func_simple)

# Calculate the best-fit model AFS.
model_simple = func_ex_simple(popt_simple, ns, pts_l)

# Likelihood of the data given the model AFS.
ll_simple = dadi.Inference.ll_multinom(model_simple, fs)

# The optimal value of theta given the model.
theta_simple = dadi.Inference.optimal_sfs_scaling(model_simple, fs)

print('\n--- simple ---');
print('Maximum log composite likelihood: {0}'.format(ll_simple))
print('Optimal value of theta: {0}'.format(theta_simple))
print('\n')

# ----------------------------------------------------------------------------------------------------------------------

# Let's generate some data using ms
mscore = dadi.Demographics2D.IM_mscore(popt_complex)
mscommand = dadi.Misc.ms_command(1., ns, mscore, int(1e5))

# We use Python's os module to call this command from within the script.
return_code = os.system('{0} > test.msout'.format(mscommand))

all_boot = dadi.Spectrum.from_ms_file('test.msout', bootstrap_segments=100)

# uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt_complex, fs,
#                                   multinom=True)
#
# # uncert contains the estimated standard deviations of each parameter, with
# # theta as the final entry in the list.
# print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))
#
# # For comparison, we can estimate uncertainties with the Fisher Information
# # Matrix, which doesn't account for linkage in the data and thus underestimates
# # uncertainty. (Although it's a fine approach if you think your data is truly
# # unlinked.)
# uncerts_fim = dadi.Godambe.FIM_uncert(func_ex, pts_l, popt_complex, fs, multinom=True)
#
# print('Estimated parameter standard deviations from FIM: {0}'.format(uncerts_fim))
# print('Factors by which FIM underestimates parameter uncertainties: {0}'.format(uncerts/uncerts_fim))

# Let's do a likelihood-ratio test comparing the complex and simple models
func_ex_simple = dadi.Numerics.make_extrap_log_func(func_simple)
model_simple = func_ex_simple(popt_simple, ns, pts_l)
ll_simple = dadi.Inference.ll_multinom(model_simple, fs)

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, popt_complex, fs, nested_indices=nested_indices, multinom=True)

print "adj = {}".format(adj)
print "ll_model = {}".format(ll_model)
print "ll_simple = {}".format(ll_simple)
print "difference = {}".format((ll_model - ll_simple))

D_adj = adj*2*(ll_model - ll_simple)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))

# Because this is test of a parameter on the boundary of parameter space
# (m cannot be less than zero), our null distribution is an even proportion
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))

print('p-value for rejecting simple model: {0:.4f}'.format(pval))