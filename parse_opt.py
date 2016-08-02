#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import dadi

from pipeline_utils import *

# load the frequency spectrum
fs = dadi.Spectrum.from_file("fsdata/DOM_14_WLD-FRE_12.fs")

ns = fs.sample_sizes

pts_l = [10,50,60]

func = dadi.Demographics2D.split_mig

upper_bound = [100, 100, 3, 10]
lower_bound = [1e-2, 1e-2, 0, 0]

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# buffer the logs created by dadi
strm = LogBuffer()
hndl = logging.StreamHandler(strm)
logger = logging.getLogger('Inference')
logger.addHandler(hndl)

with open("results.tsv", "w") as results:

    # write the header
    header = ["filename", "nu1", "nu2", "T", "m", "Maximum log composite likelihood", "Optimal value of theta", "WARNINGS"]
    results.write("\t".join(header)+"\n")

    for filename in glob.iglob('fsdata/split_mig/*.opt'):
        with open(filename) as fin:

            # reset the log
            strm.log = []

            line = fin.readline()
            while line.split()[0] != "Best-fit":
                line = fin.readline()

            line = line.replace("]", "").split()

            # get the optimal params
            popt = [eval(num) for num in line[3:]]

            # Calculate the best-fit model AFS.
            model = func_ex(popt, ns, pts_l)

            # Likelihood of the data given the model AFS.
            ll_model = dadi.Inference.ll_multinom(model, fs)

            # if "Model is masked" in " ".join(strm.log):
            #     print "{}: Model is masked".format(filename)
            #     continue

            # The optimal value of theta given the model.
            theta = dadi.Inference.optimal_sfs_scaling(model, fs)

            data = [filename] + popt + [ll_model, theta] + strm.log
            results.write("\t".join(str(x).strip("\n") for x in data)+"\n")
