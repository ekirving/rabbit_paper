"""
Custom VCF parser for generating a site frequency spectrum data file for consumption by dadi

"""
import logging
from collections import defaultdict
import pprint

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

def extract_variant_sites(population, samples, variants):

    with open('./vcf/' + population + '.vcf', 'r') as infile:

        # get the first line
        header = infile.readline()

        # skip over the block comments, stop when we reach the column header
        while header.startswith("##"):
            header = infile.readline()

        columns = header.split()

        for line in infile:

            # convert string to list
            locus = line.split()

            # skip low quality sites
            if 'LowQual' in locus[FILTER]:
                # TODO logging
                continue

            # get the reference and alternatice alleles
            ref = locus[REF]
            alt = locus[ALT].replace('.', ref)  # TODO check output from haplotype caller
            # <NON_REF>
            # T,<NON_REF>

            # skip indels
            if len(ref) > 1 or len(alt) > 1:
                # TODO logging
                continue

            # index by chromosome and position
            site = (locus[CHROM], locus[POS])

            if site not in variants:
                # initialise site locus dictionary
                variants[site] = dict()

            # initialise
            for allele in [ref, alt]:
                if allele not in variants[site]:
                    variants[site][allele] = defaultdict(int)

            # populate the variant dictionary
            for sample in samples:

                # get the genotype for the sample
                genotype = locus[columns.index(sample)][:3]

                # count the observed alleles
                if genotype == '0/0':
                    # homozygous reference
                    variants[site][ref][population] += 2
                elif genotype == '0/1':
                    # heterozygous
                    variants[site][ref][population] += 1
                    variants[site][alt][population] += 1
                elif genotype == '1/1':
                    # homozygous alternate
                    variants[site][alt][population] += 2

def generate_frequency_spectrum(populations):
    """
    Generates the site frequency spectrum for a given set of samples

    :param populations: Dictionary of populations
    :return:
    """

    # create a dictionary of variant sites
    variants = dict()

    # parse each of the populations and extract the variant sites
    for pop in populations:
        extract_variant_sites(pop, populations[pop], variants)

    pp = pprint.PrettyPrinter(depth=6)
    pp.pprint(variants)

if __name__ == '__main__':
    # TODO remove when done testing
    # population, accession codes
    POPULATIONS = dict()

    # Wild mountain hare / Lepus timidus
    POPULATIONS['OUT'] = ['SRR824842']

    # # Domestic breeds
    POPULATIONS['DOM'] = ['SRR997325', 'SRR997320'] #, 'SRR997321', 'SRR997327', 'SRR997323', 'SRR997326', 'SRR997324', 'SRR997322']
    #
    # # Wild French
    # POPULATIONS['WLD-FRE'] = ['SRR997319', 'SRR997317', 'SRR997304', 'SRR997303', 'SRR997318', 'SRR997316', 'SRR997305']
    #
    # # Wild Iberian / Oryctolagus cuniculus algirus
    # POPULATIONS['WLD-IB1'] = ['SRR827758', 'SRR827761', 'SRR827762', 'SRR827763', 'SRR827764', 'SRR827765']
    #
    # # Wild Iberian / Oryctolagus cuniculus cuniculus
    # POPULATIONS['WLD-IB2'] = ['SRR827759', 'SRR827760', 'SRR827766', 'SRR827767', 'SRR827768', 'SRR827769']

    print generate_frequency_spectrum(POPULATIONS)
