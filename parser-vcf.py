import itertools
import numpy
import re
from parser_exceptions import *
from collections import defaultdict

SAMPLES = ['SRR997317','SRR997321', 'SRR997322']
OUTGROUP = 'SRR997317'

# OUTGROUP = 'SRR997325'

# ensure that the outgroup is the first element in the list
SAMPLES.remove(OUTGROUP)
SAMPLES.insert(0, OUTGROUP)

# VCF column headers
CHROM = 0
POS = 1
ID = 2
REF = 3
ALT = 4
QUAL = 5

# the index of the first wild sample
THRESHOLD = 2

# open all files for reading
filehandles = [open("vcf/{name}.vcf.short".format(name=sample), 'r') for sample in SAMPLES]

# skip over the variable length block comments
for file in filehandles:
    while file.readline().startswith("##"):
        pass

# get the next line from all the files
for lines in itertools.izip(*filehandles):

    try:
        # convert strings to lists
        lines = [line.split() for line in lines]

        # do all the sequence positions match
        if (len(set(i[POS] for i in lines)) != 1):
            raise OutOfSyncError(lines)

        # TODO drop sites with coverage lower than 1st quartile or higer than 3rd quartile

        # get the outgroup
        outgroup = lines[0]

        # skip het sites in the outgroup
        if (len(outgroup[ALT]) > 1):
            raise HeterozygousException(outgroup)

        # get the reference and outgroup alleles
        ref_allele = outgroup[REF]
        out_allele = outgroup[ALT]

        # keep track of all the observed alleles at this site
        alleles = set([ref_allele, out_allele])

        # dictionary for counting observations
        freqs = {}
        freqs['wild'] = defaultdict(int)
        freqs['doms'] = defaultdict(int)

        for idx, line in enumerate(lines[1:]):
            # skip sites with indels
            if (len(line[REF]) != 1):
                raise InDelException(line)

            # get the alleles observed in this individual
            obs = set(line[ALT].replace('.', line[REF]))

            # add them to the alleles set
            alleles |= obs

            # skip sites with more than two alleles observed across all samples
            if (len(alleles) > 2):
                raise PolyallelicException(line)

            # which group does this sample belong to
            group = 'wild' if idx >= THRESHOLD  else 'doms'

            for ob in obs:
                # count the observations
                freqs[group][obs] +=  1


        print alleles

        # TODO what about the previous and subsequent positions



    except (InDelException, PolyallelicException, HeterozygousException) as e:
            # skip all sites containing indels,
            # polyallelic sites in ingroup samples, and
            # heterozygous sites in the outgroup
            print 'Skipping site {} {} because of {}'.format(outgroup[CHROM], outgroup[POS], type(e).__name__)