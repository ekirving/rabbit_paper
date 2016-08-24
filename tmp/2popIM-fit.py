"""
	Python script to fit the isolation with migration (IM) model to simulated data.
	Argparse is used to specify the filename where the frequency spectrum is saved.
	Requires mydemos.py in the working directory to define demographic models.
"""
import numpy
import dadi
import argparse
from time import time

T1 = time()

parser = argparse.ArgumentParser(description='Specify filename for analysis')
parser.add_argument('file')
args = parser.parse_args()
file = args.file

PODS= dadi.Spectrum.from_file(file)

nsam = PODS._get_sample_sizes()
sample_sizes = [nsam[0],nsam[1]]
pts_l = [nsam[0]+10, nsam[0]+20, nsam[0]+30]

import mydemos

IMmod = mydemos.IM
IM_ex = dadi.Numerics.make_extrap_log_func(IMmod)


upperIM = [0.99, 10, 10, 5, 20, 20]
lowerIM = [0.01, 0.5, 0.5, 0.005, .1, .1]
# lowerIM = [0.0001, 1e-2, 1e-2, 0,  0,  0]
# startIM = [0.25, 5, 5, 0.75, 2, 2]

startIM = [0.002396513,	0.010378533,	0.295010157,	0.001797582,	9.59E-07,	6.51E-05]
# startIM = [0.001632845,	0.012452512,	0.274687755,	0.001469009,	1.56E-06,	1.62E-05]

p0IM = dadi.Misc.perturb_params(startIM,fold=1,upper_bound=upperIM)


poptIM = dadi.Inference.optimize_log(p0IM, PODS, IM_ex, pts_l, verbose = 10, lower_bound = lowerIM, upper_bound = upperIM, maxiter=20, epsilon = 1e-6)

IMmod = IM_ex(poptIM, sample_sizes, pts_l)
thetaIM = dadi.Inference.optimal_sfs_scaling(IMmod, PODS)
IMmod *= thetaIM
llIM = dadi.Inference.ll(IMmod, PODS)

T2 = time()

print llIM, thetaIM, poptIM[0], poptIM[1], poptIM[2], poptIM[3], poptIM[4], poptIM[5], T2-T1