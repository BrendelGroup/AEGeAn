#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import division
from __future__ import print_function
import argparse
import pandas
import re


def cli():
    """Define the command-line interface of the program."""
    desc = 'Summarize iLocus content of the specified genome(s)'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', '--cfgdir', default=None, metavar='DIR',
                        help='directory (or comma-separated list of '
                        'directories) from which to load user-supplied genome '
                        'configuration files')
    parser.add_argument('-w', '--workdir', metavar='WD', default='./species',
                        help='working directory for data files; default is '
                        '"./species"')
    parser.add_argument('--outfmt', metavar='FMT', choices=['tsv', 'tex'],
                        default='tsv', help='output format; "tsv" for machine '
                        'readability, "tex" for typesetting')
    parser.add_argument('-s', '--shuffled', action='store_true',
                        help='load input from shuffled iLocus data')
    parser.add_argument('species', nargs='+', help='species label(s)')
    return parser


def count_seqs(data):
    """Count sequences from iLocus positions."""
    seqids = dict()
    for locuspos in data['LocusPos']:
        posmatch = re.search(r'(\S+)_(\d+)-(\d+)', locuspos)
        assert posmatch, 'error parsing iLocus position: ' + locuspos
        seqid = posmatch.group(1)
        seqids[seqid] = True
    return len(seqids)


def get_row(ilocus_data, milocus_data, fmt):
    """Calculate the summary for a row of the table."""
    assert fmt in ['tsv', 'tex']

    species = ilocus_data['Species'][0]
    miloci = milocus_data.loc[milocus_data.LocusClass == 'miLocus']
    milocus_count = len(miloci)
    effective_genome = milocus_data.loc[milocus_data.LocusClass != 'fiLocus']
    effective_genome_size = effective_genome['EffectiveLength'].sum()
    milocus_occ = miloci['EffectiveLength'].sum()
    milocus_perc = milocus_occ / effective_genome_size * 100 
    gene_count = miloci['GeneCount'].quantile([0.25, 0.50, 0.75])
    gilocus_types = ['siLocus', 'ciLocus', 'niLocus']
    singletons = milocus_data.loc[milocus_data.LocusClass.isin(gilocus_types)]
    giloci = ilocus_data.loc[ilocus_data.LocusClass.isin(gilocus_types)]
    singleton_perc = len(singletons) / len(giloci) * 100

    if fmt == 'tsv':
        genecounts = ','.join(['{:.0f}'.format(gc) for gc in gene_count])
        row = [species, milocus_count, milocus_occ, milocus_perc,
               genecounts, len(singletons), singleton_perc       ]
    elif fmt == 'tex':
        count = '{:,d}'.format(milocus_count)
        occupancy = '{:,.1f} Mb ({:.1f}\\%)'.format(milocus_occ / 1000000,
                                                    milocus_perc          )
        genecounts = ', '.join(['{:.0f}'.format(gc) for gc in gene_count])
        singles = '{:,d} ({:.1f}\\%)'.format(len(singletons),
                                             singleton_perc  )
        row = [species, count, occupancy, genecounts, singles]

    return row


def print_row(values, fmt):
    assert fmt in ['tsv', 'tex']

    if fmt == 'tsv':
        print(*values, sep='\t')
    elif fmt == 'tex':
        vals = ['{:>20}'.format(v) for v in values]
        vals[-1] += '  \\\\'
        print(*vals, sep=' & ')


def main(args):
    column_names = ['Species', 'miLoci', 'Occupancy', 'GenomePercent',
                    'GeneCountQuartiles', 'Singletons', 'SingletonPercent']
    if args.outfmt == 'tex':
        column_names = ['Species', 'miLoci', 'Occupancy', 'Gene Count',
                        'Singletons']
    print_row(column_names, args.outfmt)

    for species in args.species:
        dtype = 'loci.shuffled' if args.shuffled else 'loci'
        ilocustable = '{wd:s}/{spec:s}/{spec:s}.i{dt:s}.tsv'.format(
            wd=args.workdir, spec=species, dt=dtype,
        )
        milocustable = '{wd:s}/{spec:s}/{spec:s}.mi{dt:s}.tsv'.format(
            wd=args.workdir, spec=species, dt=dtype,
        )
        iloci = pandas.read_table(ilocustable)
        miloci = pandas.read_table(milocustable)
        row = get_row(iloci, miloci, args.outfmt)
        print_row(row, args.outfmt)


if __name__ == '__main__':
    main(args=cli().parse_args())
