"""
Custom VCF parser for generating a site frequency spectrum data file for consumption by dadi

"""
import itertools
import logging

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


class InDelException(Exception):
    """
    Detected a site which contains an insertion or deletion
    """
    pass


class PolyallelicException(Exception):
    """
    Detected a site which contains more than two alleles
    """
    pass


class HeterozygousException(Exception):
    """
    Detected a site which is heterozygous
    """
    pass


class HomozygousException(Exception):
    """
    Detected a site which is homozygous
    """
    pass


def rephase_files(lines, filehandles):
    """
    Advances the file handles until they reach consensus on the chromosome position

    :param lines: List of the current lines in each file
    :param filehandles: List of the file handles
    :return:
    """
    # get the unique set of chromosome positions
    original = positions = set([int(line[POS]) for line in lines])

    # assume no consensus
    consensus = False

    skipped = []

    while not consensus:
        # get the maximum position of the current lines
        maxpos = max(positions)

        # TODO handle case where chrom numbers are unphased

        # check each file to see which ones are behind
        for i, unused in enumerate(lines):
            # advance until we catch up
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


def generate_frequency_spectrum(samples, wild_threshold):
    """
    Generates the site frequency spectrum for a given set of samples

    :param samples: List of sample accession codes
    :param wild_threshold: The index position of the last wild sample (used for resolving group membership)
    :return:
    """

    # open all files for reading
    filehandles = [open("vcf/{}.vcf.short".format(sample), 'r') for sample in samples]

    # skip over the block comments (which are variable length)
    for fin in filehandles:
        while fin.readline().startswith("##"):
            pass

    # keep count of SNP sites
    snpcount = 0

    # initialise
    ref_allele, out_allele = None, None

    chrm, pos = 0, 0

    # initialise look-behind variables
    prev_ref_allele, prev_out_allele, prev_chrm, prev_pos = '-', '-', None, None

    # output the header row
    output = 'Rabbit\tHare\tAllele1\tWLD\tDOM\tAllele2\tWLD\tDOM\tGene\tPosition\n'

    try:

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
                chrm = int(outgroup[CHROM])
                pos = int(outgroup[POS])

                # skip all sites with indels
                if 'INDEL' in outgroup[INFO]:
                    raise InDelException(chrm, pos, outgroup[INFO])

                # get the reference and outgroup alleles
                ref_allele = outgroup[REF]
                out_allele = outgroup[ALT].replace('.', ref_allele)

                # get the genotype of the outgroup
                out_genotype = outgroup[GENOTYPE].split(':')[0]

                # skip het sites in the outgroup
                if out_genotype == '0/1':
                    raise HeterozygousException(chrm, pos, outgroup[GENOTYPE])

                # keep track of all the observed alleles at this site
                all_alleles = {ref_allele, out_allele}

                # dictionary for counting observations
                frequencies = {}

                # process all the samples (omitting the outgroup)
                for idx, line in enumerate(lines[1:]):

                    # skip all sites with indels
                    if 'INDEL' in line[INFO]:
                        raise InDelException(chrm, pos, line[REF])

                    # get the alt allele for this sample
                    alt_allele = line[ALT].replace('.', ref_allele)

                    # get the genotype of the sample
                    genotype = line[GENOTYPE].split(':')[0]

                    # resolve the genotype
                    if genotype == '0/0':
                        sample_alleles = [ref_allele, ref_allele]  # 0/0 - the sample is homozygous reference
                    elif genotype == '0/1':
                        sample_alleles = [ref_allele, alt_allele]  # 0/1 - the sample is heterozygous
                    elif genotype == '1/1':
                        sample_alleles = [alt_allele, alt_allele]  # 1/1 - the sample is homozygous alternate

                    # add them to the all alleles set
                    all_alleles |= set(sample_alleles)

                    # skip sites with more than two alleles observed across all samples
                    if len(all_alleles) > 2:
                        raise PolyallelicException(chrm, pos, all_alleles)

                    # use the index threshold to determine which group this sample belongs to
                    group = 'wild' if idx < wild_threshold else 'doms'

                    # count the observations of each allele for each group
                    for allele in sample_alleles:
                        # initialise the counter, if necessary
                        if allele not in frequencies:
                            frequencies[allele] = {'wild': 0, 'doms': 0}
                        # increment the counter
                        frequencies[allele][group] += 1

                if len(all_alleles) == 1:
                    # skip homozygous sites, because there is nothing to coalesce
                    raise HomozygousException(chrm, pos, all_alleles)

                if len(frequencies) == 1:
                    # deal with fixed allele sites by initilising the missing allele to 0
                    for allele in all_alleles:
                        if allele not in frequencies:
                            frequencies[allele] = {'wild': 0, 'doms': 0}

                # TODO what about the previous and subsequent positions

                # start composing the output line...
                # Ref | Out | Allele1 | WILD | DOMS | Allele2 | WILD | DOMS | Gene | Position
                output += '{pref}{ref}-\t{pout}{out}-\t'.format(pref=prev_ref_allele,
                                                                ref=ref_allele,
                                                                pout=prev_out_allele,
                                                                out=out_allele)

                for allele, count in frequencies.iteritems():
                    # output the allele counts
                    output += '{alle}\t{wild}\t{doms}\t'.format(alle=allele,
                                                                wild=count['wild'],
                                                                doms=count['doms'])

                # add the chromosome name and position
                output += 'chr{chrm}\t{pos}\n'.format(chrm=chrm,
                                                      pos=pos)

                # increment the SNP count
                snpcount += 1

            except (InDelException, PolyallelicException, HeterozygousException, HomozygousException) as e:
                # skip all sites containing indels, polyallelic sites in ingroup samples, heterozygous sites in the
                # outgroup, or homozygous sites across all the populations
                logging.debug('Skipping site chr{} {} because of a {} - {}'.format(outgroup[CHROM],
                                                                                   outgroup[POS],
                                                                                   type(e).__name__,
                                                                                   e))

            # update the look-behind variables if the positions are sequential,
            # otherwise we don't know what the left flanking allele is
            prev_ref_allele, prev_out_allele = (ref_allele, out_allele) if (prev_chrm, prev_pos) == (chrm, pos - 1) else ('-', '-')

            # remember the previous positions
            prev_chrm, prev_pos = chrm, pos

    except StopIteration as e:
        logging.debug('Reached the end of one of the files {}'.format(e))
        pass

    logging.debug('Finished! Found {} suitable SNP sites'.format(snpcount))

    return output

if __name__=='__main__':

    # TODO remove when done testing
    SAMPLES = {}

    # the outgroup (must be first element in the dictionary)
    SAMPLES['SRR997325'] = 'BH23'  # Belgian hare

    # wild population, w/ accession codes and sample ID
    SAMPLES['SRR997319'] = 'Avey36'  # Aveyron
    SAMPLES['SRR997317'] = 'Fos6'  # Fos-su-Mer
    SAMPLES['SRR997304'] = 'Fos2'  # Fos-su-Mer
    SAMPLES['SRR997303'] = 'Her65'  # Herauld
    SAMPLES['SRR997318'] = 'Lan7'  # Lancon
    SAMPLES['SRR997316'] = 'Lan8'  # Lancon
    SAMPLES['SRR997305'] = 'Vau73'  # Vaucluse

    # record the index of the last wild sample
    WILD_THRESHOLD = len(SAMPLES) - 1

    # domestic population, w/ accession codes and sample ID
    SAMPLES['SRR997320'] = 'FA801'  # Champagne d'argent
    SAMPLES['SRR997321'] = 'AC100'  # Angora
    SAMPLES['SRR997327'] = 'A93015'  # Angora
    SAMPLES['SRR997323'] = 'FL920'  # French lop
    SAMPLES['SRR997326'] = 'FG3'  # Flemish giant
    SAMPLES['SRR997324'] = 'FG4'  # Flemish giant
    SAMPLES['SRR997322'] = 'REX12'  # Rex

    print generate_frequency_spectrum(SAMPLES, WILD_THRESHOLD)