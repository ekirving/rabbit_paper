# Numpy is the numerical library dadi is built upon
from numpy import array

import dadi
import pylab
import os


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------


# s: Size of pop 1 after split. (Pop 2 has size 1-s.)
# nu1: Final size of pop 1.
# nu2: Final size of pop 2.
# T: Time in the past of split (in units of 2*Na generations)
# m12: Migration from pop 2 to pop 1 (2*Na*m12)
# m21: Migration from pop 1 to pop 2

# best fit params from multiple optimisation runs

# DOM /FRE
# fs = dadi.Spectrum.from_file('fsdata/all-pops_DOM_WLD-FRE.fs')
# popt_complex = [0.5, 0.036586947, 0.584179277, 0.086225099, 9.703680107, 9.951360252] # fixed s
# popt_simple = [ 0.040546754, 0.522623312, 0.024983701, 9.98911627]

    # Adjusted D statistic: 0.0186
    # p-value for rejecting simple model: 0.4458

# DOM / FRE2
# fs = dadi.Spectrum.from_file('fsdata/no-admix_DOM_WLD-FRE2.fs')
# popt_complex = [0.973950359, 0.036676267, 1.801709748, 0.119370203, 9.735916425, 9.961141793]
# popt_simple = [0.051099955, 0.308993758, 0.045069965, 9.976460257]


# DOM / FRE3
# fs = dadi.Spectrum.from_file('fsdata/no-admix_DOM_WLD-FRE3.fs')
# popt_complex = [0.905197319, 0.04151391, 0.325184361, 0.141967601, 9.900435858, 9.458934399]
# popt_simple = [0.058741985, 0.178016717, 0.059008035, 9.973357793]


# Split FRE1/2
# fs = dadi.Spectrum.from_file('fsdata/split-fre_FRE-1_FRE-2.fs')
# popt_complex = [0.013048138, 11.73178929, 1.918929748, 0.025678896, 7.545295451, 9.979874929]
# popt_simple = [0.046784111, 9.674163542, 0.008100691, 9.976127942]

    # Adjusted D statistic: 0.0891
    # p-value for rejecting simple model: 0.3827


# FRE / IB2
# fs = dadi.Spectrum.from_file('fsdata/all-pops_WLD-FRE_WLD-IB2.fs')
# popt_complex = [0.008112191, 0.256233251, 1.124115328, 0.638081552, 8.11299697, 2.115227161]
# popt_simple = [0.28917809, 1.269741392, 0.10457709, 2.836067032]

    # Adjusted D statistic: 8.8396
    # p-value for rejecting simple model: 0.0015


# IB2 / IB1
fs = dadi.Spectrum.from_file('fsdata/all-pops_WLD-IB2_WLD-IB1.fs')
popt_complex = [0.984309077, 0.275898384, 2.691461934, 0.834728827, 3.708014791, 2.07175748]
popt_simple = [0.719422317, 2.568937764, 0.477190374, 1.192350901]

    # Adjusted D statistic: 1.3382
    # p-value for rejecting simple model: 0.1237

# fs = dadi.Spectrum.from_file('fsdata/all-pops_DOM_WLD-FRE_folded.fs')
# popt_complex = [0.0216, 0.03870, 1.57378, 0.01545, 9.79209, 9.91776]
# popt_simple = [0.03324, 1.44079, 0.01804, 9.98202]

# Isolation-with-migration model with exponential pop growth.
func_complex = dadi.Demographics2D.IM

# Split into two populations of specifed size, with migration.
func_simple = dadi.Demographics2D.split_mig

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

# Let's do a likelihood-ratio test comparing the complex and simple models
func_ex_simple = dadi.Numerics.make_extrap_log_func(func_simple)
model_simple = func_ex_simple(popt_simple, ns, pts_l)
ll_simple = dadi.Inference.ll_multinom(model_simple, fs)
theta_simple = dadi.Inference.optimal_sfs_scaling(model_simple, fs)

print('\n--- simple ---');
print('Maximum log composite likelihood: {0}'.format(ll_simple))
print('Optimal value of theta: {0}'.format(theta_simple))
print('\n')


# bootstrap the spectrum
all_boot = [fs. sample () for i in range(100)]

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, p_lrt, fs, nested_indices=nested_indices, multinom=True)

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