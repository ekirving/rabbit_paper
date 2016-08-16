import dadi

# Load the data
fs = dadi.Spectrum.from_file('fsdata/all-pops_DOM_WLD-FRE.fs')
ns = fs.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [40,50,60]

# Isolation-with-migration model with exponential pop growth.
func = dadi.Demographics2D.IM

# grid to optimise over
grid = dadi.Inference.index_exp[0.0001:0.9999:100j,  # s: Size of pop 1 after split. (Pop 2 has size 1-s.)
                                0.01:100:100j,       # nu1: Final size of pop 1.
                                0.01:100:100j,       # nu2: Final size of pop 2.
                                0:3:100j,            # T: Time in the past of split (in units of 2*Na generations)
                                0:10:100j,           # m12: Migration from pop 2 to pop 1 (2*Na*m12)
                                0:10:100j]           # m21: Migration from pop 1 to pop 2

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

print('Beginning optimization ************************************************')

# do the optimization...
p1 = dadi.Inference.optimize_grid(fs, func_ex, pts_l, grid, verbose=1)

print('Finshed optimization **************************************************')

# Calculate the best-fit model AFS.
model = func_ex(p1, ns, pts_l)

# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)

# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

print('Maximum log composite likelihood: {0}'.format(ll_model))
print('Optimal value of theta: {0}'.format(theta))