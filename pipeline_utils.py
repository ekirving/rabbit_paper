#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess, datetime, hashlib, os, logging, random
from collections import defaultdict

# import all the constants
from pipeline_consts import *

def run_cmd(cmd, returnout=True, shell=False, pwd='./'):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # TODO remove when done testing
    print cmd

    # has the command so we can match the logs together
    m = hashlib.md5()
    m.update(str(cmd))

    # log the command
    with open(pwd + 'log/luigi.cmd.log', 'a+') as fout:
        fout.write("{time} {hash}# {actn}\n".format(time=str(datetime.datetime.now()),
                                                    hash=m.hexdigest(),
                                                    actn=" ".join(cmd)))

    # run the command
    proc = subprocess.Popen(cmd,
                            shell=shell,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # fetch the output and error
    (stdout, stderr) = proc.communicate()

    # bail if something went wrong
    if proc.returncode:
        raise Exception(stderr)

    if returnout:
        return stdout
    else:
        # TODO send output to a log file
        # with open('', 'w') as fout:
        pass


def unzip_file(gzip):
    """
    Unzip a gzipped file, using multi-threading when available

    :param gzip: The path to the gzip file
    :return: The uncompressed data stream
    """
    try:
        # use unpigz for multithreaded unzipping (if installed)
        return run_cmd(["unpigz",
                        "-c",                 # output to stdout
                        "-p", MAX_CPU_CORES,  # use the maximum cores
                        gzip])                # gzip file

    except OSError as e:
        # if it's not installed
        if e.errno == os.errno.ENOENT:

            # use single-threaded gunzip
            return run_cmd(["gunzip",
                            "-c",             # output to stdout
                            gzip])            # gzip file
        else:
            # escalate the exception
            raise e


def curl_download(url, filename):
    """
    Downloads a remote url to a local file path using cURL
    """

    # download the file
    run_cmd(["curl",
             "-s",                  # download silently
             "--output", filename,  # output path
             url])                  # from this url


def random_params(lower_bound, upper_bound, fixed_params):
    """
    Randomly generate starting params, within the bounding ranges, with fixed params
    """

    # make random params
    random_params = [random.uniform(lower_bound[i], upper_bound[i])
                     for i in range(0, len(upper_bound))]

    # enforce any fixed params
    if fixed_params:
        for i in range(0, len(fixed_params)):
            if fixed_params[i] is not None:
                random_params[i] = fixed_params[i]

    return random_params

def extract_variant_sites(population, samples, variants):

    # is this the outgroup population (because we don't quality filter the outgroup)
    is_outgroup = (population == 'OUT')

    # keep track of indels so we can filter SNPs within ±10 bases
    indels = []

    # parse the population VCF file
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

            # deal with the stupid <NON_REF> notation used by GATK... grrr...
            if '<NON_REF>' in alt_list:
                alt_list.remove('<NON_REF>')

            # also handle the "." notation
            if '.' in alt_list:
                alt_list.remove('.')

            # skip polyallelic sites
            if len(alt_list) > 1:
                logging.debug('{}\t{}\tPolyAllelic\t{}/{}'.format(population, site, ref, alt_list))
                continue

            # resolve empty list issue
            alt = ref if not alt_list else alt_list[0]

            # skip indels
            if len(ref) > 1 or len(alt) > 1:
                logging.debug('{}\t{}\tInDel\t{}/{}'.format(population, site, ref, alt))

                # get the size of the indel
                size = max(len(ref), len(alt))

                # remember where we found the indel and its size
                indels.append(site + (size,))
                continue

            # skip sites around indels
            if alt == "*":
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

    # filter sites within ±10 bases of each indel
    for indel in indels:
        chrom, pos, size = indel

        # filter any site within 10 bases
        sites = [(chrom, pos + offset) for offset in range(-INDEL_BUFFER, (size + INDEL_BUFFER + 1))]

        for site in sites:
            # check each allele for the current population
            for allele in list(variants.get(site, [])):
                if population in variants[site][allele]:
                    # remove the current population, but leave the others
                    del variants[site][allele][population]
                    # if this was the only population, then remove the allele entirely
                    if len(variants[site][allele]) == 0:
                        del variants[site][allele]
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
            flank1 = (locus[CHROM], int(locus[POS]) - 1) # right flank of previous site
            flank2 = (locus[CHROM], int(locus[POS]) + 1) # left flank of subsequent site

            ref = locus[REF]
            alt = locus[ALT].replace('<NON_REF>', ref).replace('.', ref)

            # skip indels and pollyallelic sites
            if len(ref) > 1 or len(alt) > 1 or alt == "*":
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

        # remove all pollyallelic sites (3 = 2 x observed alleles + 1 dict key for the ref allele)
        if len(variants[site]) > 3:
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

    # record the number of viable sites (we need this for scaling theta in the output from dadi)
    viable_sites = len(variants)

    # now lets drop all the non variant sites
    for site in list(variants):
        # (4 = 2 x observed alleles + 2 dict keys for the ref and ancestral alleles)
        if len(variants[site]) < 4:
            del variants[site]

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
    output = "# Number of viable sites: {}\n".format(viable_sites)
    output += 'Rabbit\tHare\t'

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


class LogBuffer(object):
    """
    Simple class to act as a buffer for the logging module so we can inspect warnings created by dadi
    """
    log = []

    def write(self, data):
        self.log.append(data)