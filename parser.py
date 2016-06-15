import itertools
import numpy
import re

SAMPLES = ['SRR997303', 'SRR997325']

OUTGROUP = 'SRR997325'

# ensure that the outgroup is the first element in the list
SAMPLES.remove(OUTGROUP)
SAMPLES.insert(0, OUTGROUP)

file = 'pileup/test.pileup'

# open all files for reading
filehandles = [open("pileup/{name}.test".format(name=sample), 'r') for sample in SAMPLES]

# TODO calculate 1st quartile and 3rd quartile
(q1, q2) = numpy.percentile(numpy.loadtxt('pileup/SRR997303.test', usecols=[3]), [25, 75])

# print q1, q2
# exit()

# define the positions of the columns in a pileup
SEQUENCE_ID = 0
POSITION = 1
REFERENCE = 2
COVERAGE = 3
BASES = 4
QUALITY = 5


# TODO remove me
print filehandles


class InsertionError(Exception):
    pass

class DeletionError(Exception):
    pass

class CoverageError(Exception):
    pass

class PolyallelicError(Exception):
    pass

class HetOutgroup(Exception):
    pass

def parse_alleles(bases, coverage, reference):
    """
    Extract the observed alleles from a pileup bases string

    :param bases: Bases string from a pileup file
    :param depth: The depth of coverage
    :return: List of observed alleles
    """
    if re.match('\+[0-9]+[ACGTNacgtn]+', bases):
        raise InsertionError(bases)

    if re.match('-[0-9]+[ACGTNacgtn]+', bases):
        raise DeletionError(bases)

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
        raise PolyallelicError(alleles)

    return alleles


# get the next line from all the files
for lines in itertools.izip(*filehandles):

    # convert strings to lists
    lines = [line.split() for line in lines]

    # do all the sequence positions match
    if (len(set(i[POSITION] for i in lines)) != 1):
        print "Positions don't match"

    # TODO drop sites with coverage lower than 1st quartile or higer than 3rd quartile
    # TODO drop sites with indels
    # TODO drop het sites in the outgroup

    # get the outgroup
    outgroup = lines[0]

    try:
        # get the alleles for the outgroup
        alleles = parse_alleles(outgroup[BASES], outgroup[COVERAGE], outgroup[REFERENCE])

        if len(alleles) != 1:
            raise HetOutgroup(alleles)

        # process the other samples
        for sample in lines[1:]:

            # get the alleles for the ingroup sample
            alleles = parse_alleles(sample[BASES], sample[COVERAGE], sample[REFERENCE])

    except (InsertionError, DeletionError, CoverageError, PolyallelicError):
            pass


    # print outgroup;

    # get the reference allele
    ref_allele = lines[0][2]

    # print ref_allele
    # out_allele = lines[0][3]
