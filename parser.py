import itertools
import numpy
import re
from parser_exceptions import *

SAMPLES = ['SRR997303', 'SRR997325']

OUTGROUP = 'SRR997325'

# ensure that the outgroup is the first element in the list
SAMPLES.remove(OUTGROUP)
SAMPLES.insert(0, OUTGROUP)

file = 'pileup/test.pileup'

# open all files for reading
filehandles = [open("pileup/{name}.test".format(name=sample), 'r') for sample in SAMPLES]



coverage = []

with open("pileup/SRR997303.test") as fin:
    for line in fin.readlines():
        depth = line.split()[3]
        if (depth!='1')


print coverage
exit()

# define the positions of the columns in a pileup
SEQUENCE_ID = 0
POSITION = 1
REFERENCE = 2
COVERAGE = 3
BASES = 4
QUALITY = 5

# dictionary of all samples
populations = {}

# wild population, w/ accession codes and sample ID
populations['wildfrench'] = {}
populations['wildfrench']['SRR997319']='Avey36' # Aveyron
populations['wildfrench']['SRR997317']='Fos6'   # Fos-su-Mer
populations['wildfrench']['SRR997304']='Fos2'   # Fos-su-Mer
populations['wildfrench']['SRR997303']='Her65'  # Herauld
populations['wildfrench']['SRR997318']='Lan7'   # Lancon
populations['wildfrench']['SRR997316']='Lan8'   # Lancon
populations['wildfrench']['SRR997305']='Vau73'  # Vaucluse

# domestic population, w/ accession codes and sample ID
populations['domestic'] = {}
populations['domestic']['SRR997325']='BH23'      # Belgian hare
populations['domestic']['SRR997320']='FA801'     # Champagne d'argent
populations['domestic']['SRR997321']='AC100'     # Angora
populations['domestic']['SRR997327']='A93015'    # Angora
populations['domestic']['SRR997323']='FL920'     # French lop
populations['domestic']['SRR997326']='FG3'       # Flemish giant
populations['domestic']['SRR997324']='FG4'       # Flemish giant
populations['domestic']['SRR997322']='REX12'     # Rex


# for population, samples in populations.iteritems():
#     for sample, code in samples.iteritems():
#         print sample, code

# TODO remove me
print filehandles

def parse_alleles(bases, coverage, reference):
    """
    Extract the observed alleles from a pileup bases string

    :param bases: Bases string from a pileup file
    :param depth: The depth of coverage
    :return: List of observed alleles
    """
    if re.match('\+[0-9]+[ACGTNacgtn]+', bases):
        raise InsertionException(bases)

    if re.match('-[0-9]+[ACGTNacgtn]+', bases):
        raise DeletionException(bases)

    # remove end of read and beginning of read notations
    bases = re.sub('\$|\^.', '', bases)

    if len(bases) != coverage:
        raise CoverageError(bases)

    # repalce dots and commas with the reference nucleotide
    bases = re.sub('[,.]', reference, bases)

    # standardise text case
    bases = bases.capitalize()

    # TODO what is a reasonable number of observations to consider a site to be het?

    # get the unique alleles
    alleles = set(bases)

    if len(alleles) > 2:
        raise PolyallelicException(alleles)

    return alleles


# get the next line from all the files
for lines in itertools.izip(*filehandles):

    # convert strings to lists
    lines = [line.split() for line in lines]

    # do all the sequence positions match
    if (len(set(i[POSITION] for i in lines)) != 1):
        raise OutOfSyncPileupError(alleles)

    # TODO drop sites with coverage lower than 1st quartile or higer than 3rd quartile
    # TODO drop sites with indels
    # TODO drop het sites in the outgroup

    # get the outgroup
    outgroup = lines[0]

    # get the reference allele
    ref_allele = lines[0][2]


    try:
        # get the alleles for the outgroup
        alleles = parse_alleles(outgroup[BASES], outgroup[COVERAGE], outgroup[REFERENCE])

        if len(alleles) != 1:
            raise HeterozygousException(alleles)

        # process the other samples
        for sample in lines[1:]:

            # get the alleles for the ingroup sample
            alleles = parse_alleles(sample[BASES], sample[COVERAGE], sample[REFERENCE])

    except (InsertionException, DeletionException, PolyallelicException, HeterozygousException):
            # skip all sites containing indels,
            # polyallelic sites in ingroup samples, and
            # heterozygous sites in the outgroup
            pass

