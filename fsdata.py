from pipeline_utils import *

group = 'all-pops'

# log everything to file
logging.basicConfig(filename="fsdata/{0}.log".format(group), level=logging.DEBUG)

# generate the frequency spectrum
fsdata = generate_frequency_spectrum(GROUPS[group])

# save the fsdata file
with open('fsdata.tmp', 'w') as fout:
    fout.write(fsdata)