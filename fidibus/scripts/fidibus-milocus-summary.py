#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2016   Indiana University
#
# This file is part of fidibus (http://github.com/standage/fidibus) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from __future__ import division
from __future__ import print_function
import argparse
import pandas
import re
import fidibus


def cli():
    """Define the command-line interface of the program."""
    desc = 'Summarize iLocus content of the specified genome(s)'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-v', '--version', action='version',
                        version='fidibus v%s' % fidibus.__version__)
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
    milocus_perc = milocus_occ / effective_genome_size
    gene_count = miloci['GeneCount'].quantile([0.25, 0.50, 0.75])
    gilocus_types = ['siLocus', 'ciLocus', 'niLocus']
    singletons = milocus_data.loc[milocus_data.LocusClass.isin(gilocus_types)]
    giloci = ilocus_data.loc[ilocus_data.LocusClass.isin(gilocus_types)]
    single_frac = len(singletons) / len(giloci)

    if fmt == 'tsv':
        genecounts = ','.join(['{:.0f}'.format(gc) for gc in gene_count])
        row = [species, milocus_count, milocus_occ, milocus_perc,
               genecounts, len(singletons), len(giloci)]
    elif fmt == 'tex':
        count = '{:,d}'.format(milocus_count)
        occupancy = '{:,.1f} Mb ({:.1f}\\%)'.format(milocus_occ / 1000000,
                                                    milocus_perc * 100)
        genecounts = ', '.join(['{:.0f}'.format(gc) for gc in gene_count])
        singles = '{:,d} ({:.1f}\\%)'.format(len(singletons),
                                             single_frac * 100)
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
    column_names = ['Species', 'miLoci', 'Occupancy', 'GenomeFraction',
                    'GeneCountQuartiles', 'Singletons', 'SingletonFraction']
    if args.outfmt == 'tex':
        column_names = ['Species', 'miLoci', 'Occupancy', 'Gene Count',
                        'Singletons']
    print_row(column_names, args.outfmt)

    registry = fidibus.registry.Registry()
    if args.cfgdir:
        for cfgdirpath in args.cfgdir.split(','):
            registry.update(cfgdirpath)

    for species in args.species:
        db = registry.genome(species, workdir=args.workdir)
        if args.shuffled:
            iloci = pandas.read_table(db.ilocustableshuf)
            miloci = pandas.read_table(db.milocustableshuf)
        else:
            iloci = pandas.read_table(db.ilocustable)
            miloci = pandas.read_table(db.milocustable)
        row = get_row(iloci, miloci, args.outfmt)
        print_row(row, args.outfmt)


if __name__ == '__main__':
    main(args=cli().parse_args())
