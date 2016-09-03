# load matplotlib before dadi so we can disable the screen
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

import glob, pickle, csv, dadi

group = 'all-pops'
pop1 = 'DOM'
pop2 = 'WLD-FRE'
# model = 'IM'
# scenario = 'folded-best-fit'
scenario = 'best-fit'
# pop1 = 'WLD-IB2'
# pop2 = 'WLD-IB1'
model = 'split_mig'


polarised = True

grid_size = [10, 50, 60]

# make a master list of optimal params
p_best = []

for opt in glob.glob("fsdata/opt/{}_{}_{}_{}_{}_*.opt".format(group, pop1, pop2, model, scenario)):
    print opt
    with open(opt, 'r') as fin:
        # unpickle the optimal params
        p_best += pickle.load(fin)

# sort the params (largest ll_model first)
p_best.sort(reverse=True)

param_names = ['s',    # Size of pop 1 after split. (Pop 2 has size 1-s.)
               'nu1',  # Final size of pop 1.
               'nu2',  # Final size of pop 2.
               'T',    # Time in the past of split (in units of 2*Na generations)
               'm12',  # Migration from pop 2 to pop 1 (2*Na*m12)
               'm21']  # Migration from pop 1 to pop 2


header = ["likelihood", "theta"] + list(param_names)

# dump all the data to csv
with open("fsdata/{}_{}_{}_{}_{}.csv".format(group, pop1, pop2, model, scenario), "w") as fout:
    writer = csv.writer(fout)
    writer.writerow(header)
    writer.writerows(p_best)

# get the params with the maximum log likelihood, from all ~1e6 iterations
p_opt = p_best[0][2:]

polar = '_folded' if not polarised else ''

fs = dadi.Spectrum.from_file("fsdata/{}_{}_{}{}.fs".format(group, pop1, pop2, polar))
# load the frequency spectrum
ns = fs.sample_sizes

# get the demographic model to test
func = getattr(dadi.Demographics2D, model)

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Calculate the best-fit model AFS.
md = func_ex(p_opt, ns, grid_size)

# plot the figure
fig = plt.figure(1)
dadi.Plotting.plot_2d_comp_multinom(md, fs, fig_num=1)
fig.savefig("fsdata/{}_{}_{}_{}_{}.pdf".format(group, pop1, pop2, model, scenario))
plt.close(fig)