#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing

# the reference genome
GENOME = "OryCun2.0"
GENOME_URL = "ftp://ftp.ensembl.org/pub/release-84/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0" \
             ".dna.toplevel.fa.gz"

# file containing the list of sequence capture regions, format <chr>:<start>-<stop>
TARGETS_LIST = './targets.interval_list'

# populations and sample accession codes (n=28)
POPULATIONS = {

    # Wild mountain hare / Lepus timidus (n=1)
    'OUT': ['SRR824842'],

    # Domestic breeds (n=8)
    'DOM': ['SRR997325', 'SRR997320', 'SRR997321', 'SRR997327', 'SRR997323', 'SRR997326', 'SRR997324', 'SRR997322'],

    # Wild French (n=7)
    'WLD-FRE': ['SRR997319', 'SRR997317', 'SRR997304', 'SRR997303', 'SRR997318', 'SRR997316', 'SRR997305'],

    # Wild Iberian / Oryctolagus cuniculus algirus (n=6)
    'WLD-IB1': ['SRR827758', 'SRR827761', 'SRR827762', 'SRR827763', 'SRR827764', 'SRR827765'],

    # Wild Iberian / Oryctolagus cuniculus cuniculus (n=6)
    'WLD-IB2': ['SRR827759', 'SRR827760', 'SRR827766', 'SRR827767', 'SRR827768', 'SRR827769']
}

# accession codes and sample IDs
SAMPLES = {

    # Wild hare, outgroup
    'SRR824842': 'LT2012',  # Lepus timidus

    # domestic population, w/ accession codes and sample ID
    'SRR997320': 'FA801',   # Champagne d'argent
    'SRR997321': 'AC100',   # Angora
    'SRR997322': 'REX12',   # Rex
    'SRR997323': 'FL920',   # French lop
    'SRR997324': 'FG4',     # Flemish giant
    'SRR997325': 'BH23',    # Belgian hare
    'SRR997326': 'FG3',     # Flemish giant
    'SRR997327': 'A93015',  # Angora

    # wild French population, w/ accession codes and sample ID
    'SRR997303': 'Her65',   # Herauld
    'SRR997304': 'Fos2',    # Fos-su-Mer
    'SRR997305': 'Vau73',   # Vaucluse
    'SRR997316': 'Lan8',    # Lancon
    'SRR997317': 'Fos6',    # Fos-su-Mer
    'SRR997318': 'Lan7',    # Lancon
    'SRR997319': 'Avey36',  # Aveyron

    # wild iberian 1, Oryctolagus cuniculus algirus
    'SRR827758': 'Hue52',   # Huelva
    'SRR827761': 'Pan2',    # Pancas
    'SRR827762': 'Ped16',   # Sevilla
    'SRR827763': 'Ped8',    # Sevilla
    'SRR827764': 'Pfr11',   # Sevilla
    'SRR827765': 'Pfr15',   # Sevilla

    # wild iberian 2, Oryctolagus cuniculus cuniculus
    'SRR827759': 'Lle10',   # Lleida
    'SRR827760': 'Lle11',   # Lleida
    'SRR827766': 'Zrg13',   # Zaragoza
    'SRR827767': 'Zrg2',    # Zaragoza
    'SRR827768': 'Zrg3',    # Zaragoza
    'SRR827769': 'Zrg7',    # Zaragoza
}

# setup some population groups
GROUPS = dict()

# one with all the populations
GROUPS['all-pops'] = POPULATIONS.copy()

# one without the outgroup
GROUPS['no-outgroup'] = POPULATIONS.copy()
del GROUPS['no-outgroup']['OUT']

# and one without either the outgroup or the most divergent wild species (O. cuniculus algirus)
GROUPS['no-out-ib1'] = POPULATIONS.copy()
del GROUPS['no-out-ib1']['OUT']
del GROUPS['no-out-ib1']['WLD-IB1']

# the samtools flag for BAM file comression
DEFAULT_COMPRESSION = 6

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_GENOTYPE_QUAL = 30

# use a lower threshold for the outgroup to maintain sufficient coverage
OUTGROUP_MIN_GQ = 20

# the minimum depth of coverage for a site
MIN_COVERAGE_DEPTH = 8

# the maximum number of ancestral populatons to run admiture for
MAX_ANCESTRAL_K = 10

# no single worker should use more than 50% of the available cores
MAX_CPU_CORES = int(multiprocessing.cpu_count() * 0.5)

# size of buffer around indels in which to drop sites
INDEL_BUFFER = 10

# location of software tools
PICARD = "/usr/local/picard-tools-2.5.0/picard.jar"
GATK = "/usr/local/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
SAMTOOLS = "/usr/local/bin/samtools1.3"

# how many iterations to use to optimise params
DADI_MAX_ITER = 100

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