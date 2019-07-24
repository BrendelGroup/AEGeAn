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
import sys
import fidibus


def cli():
    """Define the command-line interface of the program."""
    desc = 'Summarize piLocus content of the specified genome(s)'
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
    parser.add_argument('species', nargs='+', help='species label(s)')
    return parser


def get_row(iloci, premrnas, fmt):
    """Calculate the summary for a row of the table."""
    assert fmt in ['tsv', 'tex']

    species = iloci['Species'][0]
    piloci = iloci.loc[iloci.LocusClass.isin(['siLocus', 'ciLocus'])]
    pilocus_count = len(piloci)
    effective_genome = iloci.loc[iloci.LocusClass != 'fiLocus']
    effective_genome_size = effective_genome['EffectiveLength'].sum()
    pilocus_occ = piloci['EffectiveLength'].sum()
    pilocus_occ_perc = pilocus_occ / effective_genome_size
    single_exon_piloci = len(premrnas[premrnas['ExonCount'] == 1])
    single_exon_perc = single_exon_piloci / len(premrnas)

    if fmt == 'tsv':
        row = [species, pilocus_count, pilocus_occ, pilocus_occ_perc,
               single_exon_piloci, single_exon_perc]
    elif fmt == 'tex':
        count = '{:,d}'.format(pilocus_count)
        occupancy = '{:,.1f} Mb ({:.1f}\\%)'.format(pilocus_occ / 1000000,
                                                    pilocus_occ_perc * 100)
        sepiloci = '{:,d} ({:.1f}\\%)'.format(single_exon_piloci,
                                              single_exon_perc * 100)
        row = [species, count, occupancy, sepiloci]

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
    column_names = ['Species', 'piLoci', 'Occupancy', 'GenomeFraction',
                    'SingleExon_piLoci', 'SingleExonFraction']
    if args.outfmt == 'tex':
        column_names = ['Species', 'piLoci', 'Occupancy',
                        'Single Exon piLoci']
    print_row(column_names, args.outfmt)

    registry = fidibus.registry.Registry()
    if args.cfgdir:
        for cfgdirpath in args.cfgdir.split(','):
            registry.update(cfgdirpath)

    for species in args.species:
        db = registry.genome(species, workdir=args.workdir)
        iloci = pandas.read_table(db.ilocustable)
        premrnas = pandas.read_table(db.premrnatable)
        row = get_row(iloci, premrnas, args.outfmt)
        print_row(row, args.outfmt)


if __name__ == '__main__':
    main(args=cli().parse_args())
