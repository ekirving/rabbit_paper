import dadi
import numpy

# Parse the data file to generate the data dictionary
dd = dadi.Misc.make_data_dict('./fsdata/OryCun2.0.data')

data = dadi.Spectrum.from_data_dict(dd, ['DOM', 'WLD-FRE'], [16, 14], polarized=True)


# import pylab
# # # dadi.Plotting.plot_1d_fs(data)
# dadi.Plotting.plot_single_2d_sfs(data)
# pylab.show()
# exit()

ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [10,50,60]

# The Demographics1D and Demographics2D modules contain a few simple models,
# mostly as examples. We could use one of those.
func = dadi.Demographics2D.split_mig

# Now let's optimize parameters for this model.

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.
# Parameters are: (nu1F, nu2B, nu2F, m, Tp, T)
# params = (nu1,nu2,T,m)

upper_bound = [100, 100, 3, 10]
lower_bound = [1e-2, 1e-2, 0, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [2,0.1,0.2,0.2]

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Perturb our parameters before optimization. This does so by taking each
# parameter a up to a factor of two up or down.
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                              lower_bound=lower_bound)
# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters. If you want to fix
# theta, add a multinom=False to the call.
# The maxiter argument restricts how long the optimizer will run. For real
# runs, you will want to set this value higher (at least 10), to encourage
# better convergence. You will also want to run optimization several times
# using multiple sets of intial parameters, to be confident you've actually
# found the true maximum likelihood parameters.
# print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=20)

# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finshed optimization **************************************************')

print('Best-fit parameters: {0}'.format(popt))

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

# Plot a comparison of the resulting fs with the data.
import pylab
pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,
                                    pop_ids =('W','D'))
# This ensures that the figure pops up. It may be unecessary if you are using
# ipython.
pylab.show()
# Save the figure

# # ----------------------------------------------------------------------------------------------------------------------
#
func_mig=dadi.Demographics2D.IM
#ns = (n1,n2)
#params = (s,nu1,nu2,T,m12,m21)

upper_bound_1 = [1, 100, 100, 10, 3, 3]
lower_bound_1 = [0, 1e-2, 1e-2, 0, 0, 0]
p1 = [0.5, 2, 0.1, 0.2, 0.2, 0.2]

func_ex_mig = dadi.Numerics.make_extrap_log_func(func_mig)

# for i in range(0,9):
p1 = dadi.Misc.perturb_params(p1, fold=1, upper_bound=upper_bound_1,
                              lower_bound=lower_bound_1)

print('Beginning optimization ************************************************')
popt1 = dadi.Inference.optimize_log(p1, data, func_ex_mig, pts_l,
                                   lower_bound=lower_bound_1,
                                   upper_bound=upper_bound_1,
                                   # verbose=len(p1),
                                    maxiter=10)
# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finshed optimization **************************************************')

print('Best-fit parameters: {0}'.format(popt1))

# Calculate the best-fit model AFS.
model_mig = func_ex_mig(popt1, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model_mig = dadi.Inference.ll_multinom(model_mig, data)
print('Maximum log composite likelihood: {0}'.format(ll_model_mig))
# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model_mig, data)
print('Optimal value of theta: {0}'.format(theta))

pylab.figure()
dadi.Plotting.plot_2d_comp_multinom(model_mig, data, vmin=1, resid_range=3,
                                    pop_ids =('W','D'))
# This ensures that the figure pops up. It may be unecessary if you are using
# ipython.
pylab.show()