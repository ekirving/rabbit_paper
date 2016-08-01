import glob
import dadi

import logging
logger = logging.getLogger('Inference')

# load the frequency spectrum
fs = dadi.Spectrum.from_file("fsdata/DOM_14_WLD-FRE_12.fs")

ns = fs.sample_sizes

pts_l = [10,50,60]

func = dadi.Demographics2D.split_mig

upper_bound = [100, 100, 3, 10]
lower_bound = [1e-2, 1e-2, 0, 0]

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)


class writer(object):
    log = []

    def write(self, data):
        self.log.append(data)

strm = writer()
hndl = logging.StreamHandler(strm)
logger.addHandler(hndl)

with open("results.tsv", "w") as results:

    # write the header
    header = ["nu1", "nu2", "T", "m", "Maximum log composite likelihood", "Optimal value of theta", "WARNINGS"]
    results.write("\t".join(header)+"\n")

    for filename in glob.iglob('fsdata/split_mig/*.opt'):
        with open(filename) as fin:

            # reset the log
            strm.log = []

            line = fin.readline().split()

            while line[0] != "Best-fit":
                line = fin.readline().split()

            # get the optimal params
            popt = [num.rstrip(']') for num in line[3:]]

            # Calculate the best-fit model AFS.
            model = func_ex(popt, ns, pts_l)

            # Likelihood of the data given the model AFS.
            ll_model = dadi.Inference.ll_multinom(model, fs)

            # if "Model is masked" in " ".join(strm.log):
            #     print "{}: Model is masked".format(filename)
            #     continue

            # The optimal value of theta given the model.
            theta = dadi.Inference.optimal_sfs_scaling(model, fs)

            data = popt + [ll_model, theta] + strm.log
            results.write("\t".join(str(x).strip("\n") for x in data)+"\n")
