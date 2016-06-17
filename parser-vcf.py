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


def rephase_files(lines, filehandles):
    """
    Advances the file handles until they reach consensus on the chromosome position

    :param lines:
    :param filehandles:
    :return:
    """
    # get the unique set of positions
    original = positions = set([int(line[POS]) for line in lines])

    # assume no consensus
    consensus = False


    skipped = []

    while not consensus:
        # get the maximum position of the current lines
        maxpos = max(positions)

        for i, unused in enumerate(lines):
            # check each file to see which ones are behind
            while int(lines[i][POS]) < maxpos:
                # keep track of skipped sites
                skipped.append(lines[i][POS])
                # advance to the next line
                lines[i] = filehandles[i].next().split()

        # get the unique set of positions
        positions = set([int(line[POS]) for line in lines])

        # do we have consensus
        if len(positions) == 1:
            consensus = True

    print "Rephased from {} to {}, skipping {} sites ".format(original, positions, len(set(skipped)))

# skip over the variable length block comments
for file in filehandles:
    while file.readline().startswith("##"):
        pass

# TODO remove when done testing
count = 0

try :

    # get the next line from all the files
    for lines in itertools.izip(*filehandles):

        try:
            # convert strings to lists
            lines = [line.split() for line in lines]

            # repahse the files if not all the sequence positions match
            if (len(set(line[POS] for line in lines)) != 1):
                rephase_files(lines, filehandles)

            # TODO drop sites with coverage lower than 1st quartile or higer than 3rd quartile

            # get the outgroup
            outgroup = lines[0]

            # skip het sites in the outgroup
            if (len(outgroup[ALT]) > 1):
                raise HeterozygousException(outgroup[ALT])

            # get the reference and outgroup alleles
            ref_allele = outgroup[REF]
            out_allele = outgroup[ALT].replace('.', ref_allele)

            # keep track of all the observed alleles at this site
            alleles = set([ref_allele, out_allele])

            # dictionary for counting observations
            freqs = {}
            freqs['wild'] = defaultdict(int) # default all entries to integers
            freqs['doms'] = defaultdict(int)

            for idx, line in enumerate(lines[1:]):
                # skip sites with indels
                if (len(line[REF]) != 1):
                    raise InDelException(line[REF])

                # get the alleles observed in this individual
                obs = set(line[ALT].replace('.', line[REF]))

                # add them to the alleles set
                alleles |= obs

                # skip sites with more than two alleles observed across all samples
                if (len(alleles) > 2):
                    raise PolyallelicException(alleles)

                # which group does this sample belong to
                group = 'wild' if idx >= THRESHOLD  else 'doms'

                # count the observations of each allele for each group
                for ob in obs:
                    freqs[group][ob] +=  1

            if len(alleles) == 1:
                raise HomozygousException(alleles)

            print freqs

            # TODO what about the previous and subsequent positions

            # TODO remove when done testing
            count += 1
            if (count > 100):
                break

        except (InDelException, PolyallelicException, HeterozygousException, HomozygousException) as e:
                # skip all sites containing indels, polyallelic sites in ingroup samples, heterozygous sites in the outgroup,
                # or homozygous sites across the all populations
                # print 'Skipping site {} {} because of {} - {}'.format(outgroup[CHROM], outgroup[POS], type(e).__name__, e)
                pass

except StopIteration as e:
    print 'Reached the end of one of the files {}'.format(e)
    pass

print 'Finished!'