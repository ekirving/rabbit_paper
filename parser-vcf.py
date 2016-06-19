import itertools
import numpy
import re
from parser_exceptions import *
from collections import defaultdict
import logging

# log to file
logging.basicConfig(filename='parse-vcf.log',level=logging.DEBUG)

SAMPLES = {}

# outgroup (must be first element in dictionary)
SAMPLES['SRR997325']='BH23'   # Belgian hare

# wild population, w/ accession codes and sample ID
SAMPLES['SRR997319']='Avey36' # Aveyron
SAMPLES['SRR997317']='Fos6'   # Fos-su-Mer
SAMPLES['SRR997304']='Fos2'   # Fos-su-Mer
SAMPLES['SRR997303']='Her65'  # Herauld
SAMPLES['SRR997318']='Lan7'   # Lancon
SAMPLES['SRR997316']='Lan8'   # Lancon
SAMPLES['SRR997305']='Vau73'  # Vaucluse

# record the index of the last wild sample
WILD_THRESHOLD = len(SAMPLES)-1

# domestic population, w/ accession codes and sample ID
SAMPLES['SRR997320']='FA801'   # Champagne d'argent
SAMPLES['SRR997321']='AC100'   # Angora
SAMPLES['SRR997327']='A93015'  # Angora
SAMPLES['SRR997323']='FL920'   # French lop
SAMPLES['SRR997326']='FG3'     # Flemish giant
SAMPLES['SRR997324']='FG4'     # Flemish giant
SAMPLES['SRR997322']='REX12'   # Rex

# VCF column headers
CHROM = 0
POS = 1
ID = 2
REF = 3
ALT = 4
QUAL = 5
FILTER = 6
INFO = 7
FORMAT = 8
GENOTYPE = 9

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

    logging.debug("Rephased from {} to {}, skipping {} sites ".format(original, positions, len(set(skipped))))

# skip over the variable length block comments
for file in filehandles:
    while file.readline().startswith("##"):
        pass

# keep count of SNP sites
snpcount = 0

try :

    # output the header row
    print 'Rabbit\tHare\tAllele1\tWLD\tDOM\tAllele2\tWLD\tDOM\tGene\tPosition'

    # get the next line from all the files
    for lines in itertools.izip(*filehandles):

        try:
            # convert each line from a string to a list
            lines = [line.split() for line in lines]

            # rephase the files, if not all the sequence positions match
            if len(set(line[POS] for line in lines)) != 1:
                rephase_files(lines, filehandles)

            # TODO drop sites with coverage lower than 1st quartile or higher than 3rd quartile

            # get the outgroup
            outgroup = lines[0]

            # get the chromosome number and position
            chr = outgroup[CHROM]
            pos = outgroup[POS]

            # skip all sites with indels
            if 'INDEL' in outgroup[INFO]:
                raise InDelException(chr, pos, outgroup[INFO])

            # get the reference and outgroup alleles
            ref_allele = outgroup[REF]
            out_allele = outgroup[ALT].replace('.', ref_allele)

            # get the genotype of the outgroup
            out_genotype = outgroup[GENOTYPE].split(':')[0]

            # skip het sites in the outgroup
            if out_genotype == '0/1':
                raise HeterozygousException(chr, pos, outgroup[GENOTYPE])

            # keep track of all the observed alleles at this site
            all_alleles = set([ref_allele, out_allele])

            # dictionary for counting observations
            frequencies = {}

            # process all the samples (omitting the outgroup)
            for idx, line in enumerate(lines[1:]):

                # skip all sites with indels
                if 'INDEL' in line[INFO]:
                    raise InDelException(chr, pos, line[REF])

                # get the alt allele for this sample
                alt_allele = line[ALT].replace('.', ref_allele)

                # get the genotype of the sample
                genotype = line[GENOTYPE].split(':')[0]

                # resolve the genotype
                if genotype == '0/0':
                    sample_alleles = [ref_allele, ref_allele] # 0/0 - the sample is homozygous reference
                elif genotype == '0/1':
                    sample_alleles = [ref_allele, alt_allele] # 0/1 - the sample is heterozygous
                elif genotype == '1/1':
                    sample_alleles = [alt_allele, alt_allele] # 1/1 - the sample is homozygous alternate

                # add them to the all alleles set
                all_alleles |= set(sample_alleles)

                # skip sites with more than two alleles observed across all samples
                if len(all_alleles) > 2:
                    raise PolyallelicException(chr, pos, all_alleles)

                # use the index threshold to determine which group this sample belongs to
                group = 'wild' if idx < WILD_THRESHOLD else 'doms'

                # count the observations of each allele for each group
                for allele in sample_alleles:
                    # initialise the counter, if necessary
                    if allele not in frequencies:
                        frequencies[allele] = {'wild': 0, 'doms':0}
                    # increment the counter
                    frequencies[allele][group] += 1

            if len(all_alleles) == 1:
                # skip homozygous sites, because there is nothing to coalesce
                raise HomozygousException(chr, pos, all_alleles)

            if len(frequencies) == 1:
                # deal with fixed allele sites by initilising the missing allele to 0
                for allele in all_alleles:
                    if allele not in frequencies:
                        frequencies[allele] = {'wild': 0, 'doms':0}

            # TODO what about the previous and subsequent positions

            # start composing the output line...
            # Ref | Out | Allele1 | WILD | DOMS | Allele2 | WILD | DOMS | Gene | Position
            output = '-{ref}-\t-{out}-\t'.format(ref=ref_allele,
                                                 out=out_allele)

            for allele, count in frequencies.iteritems():
                # output the allele counts
                output += '{alle}\t{wild}\t{doms}\t'.format(alle=allele,
                                                            wild=count['wild'],
                                                            doms=count['doms'])

            # add the chromosome name and position
            output += 'chr{chr}\t{pos}'.format(chr=chr,
                                               pos=pos)

            # print the output
            print output

            # increment the SNP count
            snpcount += 1

        except (InDelException, PolyallelicException, HeterozygousException, HomozygousException) as e:
                # skip all sites containing indels, polyallelic sites in ingroup samples, heterozygous sites in the
                # outgroup, or homozygous sites across all the populations
                logging.debug('Skipping site chr{} {} because of a {} - {}'.format(outgroup[CHROM],
                                                                                   outgroup[POS],
                                                                                   type(e).__name__,
                                                                                   e))

except StopIteration as e:
    logging.debug('Reached the end of one of the files {}'.format(e))
    pass

logging.debug('Finished! Found {} suitable SNP sites'.format(snpcount))