#!/usr/bin/python
"""
Custom VCF parser for generating a site frequency spectrum data file for consumption by dadi

"""

import logging
from collections import defaultdict

# the minimum depth of coverage for a site
MIN_COVERAGE_DEPTH = 8

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

    # is this the outgroup population (because we don't quality filer the outgroup)
    is_outgroup = (population == 'OUT')

    # keep track of indels so we can filter SNPs within +/-10 bases
    indels = []

    with open('./vcf/' + population + '.vcf', 'r') as infile:

        # get the first line
        header = infile.readline()

        # skip over the block comments, stop when we reach the column header
        while header.startswith("##"):
            header = infile.readline()

        # get the column headers
        columns = header.split()

        for line in infile:

            # convert string to list
            locus = line.split()

            # index by chromosome and position
            site = (locus[CHROM], int(locus[POS]))

            # skip low quality sites
            if 'LowQual' in locus[FILTER] and not is_outgroup:
                logging.debug('{}\t{}\tLowQual\t{}'.format(population, site, locus[QUAL]))
                continue

            # extract the locus info
            info = dict(item.split("=") for item in locus[INFO].split(";"))

            # get the joint depth (handle sites with no coverage)
            joint_depth = int(info['DP']) if 'DP' in info else 0

            # skip low coverage sites (average depth must be > MIN_COVERAGE_DEPTH)
            if joint_depth/len(samples) < MIN_COVERAGE_DEPTH and not is_outgroup:
                logging.debug('{}\t{}\tLowDepth\t{}/{}'.format(population, site, joint_depth, len(samples)))
                continue

            # get the reference allele
            ref = locus[REF]

            # get the alternative allele(s)
            alt_list = locus[ALT].split(',')

            # deal with the stupid <NON_REF> notation used by BP_RESOLUTION... grrr...
            if '<NON_REF>' in alt_list:
                alt_list.remove('<NON_REF>')

            # skip polyallelic sites
            if len(alt_list) > 1:
                logging.debug('{}\t{}\tPolyAllelic\t{}/{}'.format(population, site, ref, alt_list))
                continue

            # resolve empty list issue
            alt = ref if not alt_list else alt[0]

            # skip indels
            if len(ref) > 1 or len(alt) > 1:
                logging.debug('{}\t{}\tInDel\t{}/{}'.format(population, site, ref, alt))

                # remember where we found the indel
                indels.append(site)
                continue

            if site not in variants:
                # initialise site locus dictionary
                variants[site] = dict()

            if 'ref' not in variants[site]:
                # remember the reference allele
                variants[site]['ref'] = ref

            for allele in [ref, alt]:
                # initialise the allele count dictionary
                if allele not in variants[site]:
                    variants[site][allele] = defaultdict(int)

            # get the indices for the genotype column (e.g. GT:AD:DP:GQ:PL:SB)
            format = locus[FORMAT].split(':')

            GT = format.index('GT') # Genotype (1/1, 0/0, 0/1)

            # populate the variant dictionary
            for sample in samples:

                # get the genotype for the sample
                genotype = locus[columns.index(sample)].split(':')

                # count the observed alleles
                if genotype[GT] == '0/0':
                    # homozygous reference
                    variants[site][ref][population] += 2

                elif genotype[GT] == '0/1':
                    # heterozygous
                    variants[site][ref][population] += 1
                    variants[site][alt][population] += 1

                elif genotype[GT] == '1/1':
                    # homozygous alternate
                    variants[site][alt][population] += 2

    # filter sites within +/- 10 bases of each idels
    for indel in indels:
        chrom, pos = indel

        # filter any site within 10 bases
        sites = [(chrom, pos + offset) for offset in range(-10, 11)]

        for site in sites:
            # remove the current population, but leave the others
            for allele in variants.get(site, []):
                if population in variants[site][allele]:
                    del variants[site][allele][population]
                logging.debug('{}\t{}\tInDelProximity\t{}'.format(population, site, indel))

def find_flanking_bases(variants):

    with open('./vcf/OUT.vcf', 'r') as infile:

        # skip over the block comments
        while infile.readline().startswith("##"):
            pass

        # check all the loci
        for line in infile:

            # convert string to list
            locus = line.split()

            # get the positions of the adjacent flanks
            flank1 = (locus[CHROM], int(locus[POS]) - 1) # right flank
            flank2 = (locus[CHROM], int(locus[POS]) + 1) # left flank

            ref = locus[REF]
            alt = locus[ALT].replace('<NON_REF>', ref).replace('*', ref)

            # skip indels and pollyallelic sites
            if len(ref) > 1 or len(alt) > 1:
                continue

            if flank1 in variants:
                variants[flank1]['ref_rgt'] = ref
                variants[flank1]['alt_rgt'] = alt

            if flank2 in variants:
                variants[flank2]['ref_lft'] = ref
                variants[flank2]['alt_lft'] = alt


def generate_frequency_spectrum(populations):
    """
    Generates the site frequency spectrum for a given set of populations

    :param populations: Dictionary of populations
    :return:
    """

    # create a dictionary of variant sites
    variants = dict()

    # parse each of the populations and extract the variant sites
    for pop in populations:
        extract_variant_sites(pop, populations[pop], variants)

    # delete all the sites that don't conform to our requirements
    for site in list(variants):

        # remove all non variant and pollyallelic sites (2 alleles + 1 reference = 3)
        if len(variants[site]) != 3:
            del variants[site]
            continue

        # find the ancestral allele(s)
        ancestral = [allele for allele in variants[site] if 'OUT' in variants[site][allele]]

        # skip all het sites in the outgroup and those without coverage
        if len(ancestral) != 1:
            del variants[site]
            continue

        # remember the ancestral allele
        variants[site]['alt'] = ancestral[0]

    # now we've whittled down the sites, lets find the flanking bases
    find_flanking_bases(variants)

    # convert to a list, because dictionaries are not sorted and we need the iteration to be consistent
    site_list = list(variants)
    site_list.sort()

    pop_list = list(populations)
    pop_list.sort()

    # remove the pseudo-population OUT
    pop_list.remove('OUT')

    # start composing the output file
    output = 'Rabbit\tHare\t'

    for num in [1,2]:
        output += 'Allele{}\t'.format(num)

        for pop in pop_list:
            output += '{}\t'.format(pop)

    output += 'Chrom\tPos\n'

    for site in site_list:

        # get the chrom and pos
        chrom, pos = site

        # get the ref and alt alleles
        ref = variants[site].pop('ref')
        ref_lft = variants[site].pop('ref_lft', '-')
        ref_rgt = variants[site].pop('ref_rgt', '-')

        alt = variants[site].pop('alt')
        alt_lft = variants[site].pop('alt_lft', '-')
        alt_rgt = variants[site].pop('alt_rgt', '-')

        # output the alleles and their flanking bases
        output += '{}{}{}\t'.format(ref_lft, ref, ref_rgt)
        output += '{}{}{}\t'.format(alt_lft, alt, alt_rgt)

        # make sure the ref allele is first in the list
        alleles = list(variants[site])
        alleles.insert(0, alleles.pop(alleles.index(ref)))

        # output each of the alleles
        for allele in alleles:
            output += '{}\t'.format(allele)

            # output the allele count for each population
            for pop in pop_list:
                # note that SNPs cannot be projected up, so SNPs without enough calls in any population will be ignored
                # https://bitbucket.org/gutenkunstlab/dadi/wiki/DataFormats
                count = variants[site][allele][pop] if pop in variants[site][allele] else 0
                output += '{}\t'.format(count)

        # output the chromosome and position of the SNP
        output += 'chr{}\t{}\n'.format(chrom, pos)

    return output
