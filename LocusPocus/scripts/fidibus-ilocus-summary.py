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
    parser.add_argument('-w', '--workdir', metavar='WD', default='./species',
                        help='working directory for data files; default is '
                        '"./species"')
    parser.add_argument('--outfmt', metavar='FMT', choices=['tsv', 'tex'],
                        default='tsv', help='output format; "tsv" for machine '
                        'readability, "tex" for typesetting')
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


def get_row(data, fmt):
    """Calculate the summary for a row of the table."""
    assert fmt in ['tsv', 'tex']

    row = [
        data['Species'][0],
        data['EffectiveLength'].sum() / 1000000,
        count_seqs(data),
        len(data.loc[data.LocusClass == 'fiLocus']),
        len(data.loc[data.LocusClass == 'iiLocus']),
        len(data.loc[data.LocusClass == 'niLocus']),
        len(data.loc[data.LocusClass == 'siLocus']),
        len(data.loc[data.LocusClass == 'ciLocus']),
    ]

    if fmt == 'tex':
        row[1] = '{:.1f}'.format(row[1])
        for i in range(2, len(row)):
            row[i] = '{:,d}'.format(row[i])

    return row


def print_row(values, fmt):
    assert fmt in ['tsv', 'tex']

    if fmt == 'tsv':
        print(*values, sep='\t')
    elif fmt == 'tex':
        vals = ['{:>12}'.format(v) for v in values] + [' \\\\']
        print(*vals, sep=' & ')


def main(args):
    column_names = ['Species', 'Mb', '#Seq', 'fiLoci', 'iiLoci', 'niLoci',
                    'siLoci', 'ciLoci']
    print_row(column_names, args.outfmt)

    for species in args.species:
        ilocustable = '{wd:s}/{spec:s}/{spec:s}.iloci.tsv'.format(
            wd=args.workdir, spec=species
        )
        data = pandas.read_table(ilocustable)
        row = get_row(data, args.outfmt)
        print_row(row, args.outfmt)


if __name__ == '__main__':
    main(args=cli().parse_args())
