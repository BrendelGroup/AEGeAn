#!/usr/bin/env python

# Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS
#
# The AEGeAn Toolkit is distributed under the ISC License. See
# the 'LICENSE' file in the AEGeAn source code distribution or
# online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

from __future__ import print_function
import re
import sys


def merge_iloci(loci):
    """Merge ajacent or overlapping gene-containing iLoci."""
    assert len(loci) > 0
    if len(loci) == 1:
        line = re.sub('ID=[^;\n]+;*', '', loci[0])
        line = re.sub('Name=[^;\n]+;*', '', line)
        return line

    seqid = None
    start, end = -1, -1
    attrs = {}
    for locus in loci:
        fields = locus.split('\t')
        assert len(fields) == 9
        if seqid:
            assert fields[0] == seqid
        seqid = fields[0]
        lstart = int(fields[3])
        lend = int(fields[4])
        if start == -1 or lstart < start:
            start = lstart
        end = max(end, lend)
        numeric_attrs = re.findall('([^;=]+=\d+)', fields[8])
        for key_value_pair in numeric_attrs:
            assert '=' in key_value_pair, \
                'malformed key/value pair %s' % key_value_pair
            key, value = key_value_pair.split('=')
            if key in ['left_overlap', 'right_overlap']:
                continue
            value = int(value)
            if key not in attrs:
                attrs[key] = 0
            attrs[key] += value

    attrstring = 'merged=true;iLocus_type=miLocus'
    for key in attrs:
        attrstring += ';%s=%d' % (key, attrs[key])
    gff3 = [seqid, 'AEGeAn::miloci.py', 'locus', str(start), str(end),
            '%d' % len(loci), '.', '.',    attrstring]
    return '\t'.join(gff3)


def parse_iloci(fp):
    """
    Input: a GFF3 file containing iLoci (LocusPocus output)
    Output: merged iLoci; gene-containing iLoci that are adjacent or
            overlapping are combined
    """
    seqid = None
    prev_loci = []
    for line in fp:
        line = line.rstrip()
        if '\tlocus\t' not in line:
            continue

        locusseqid = re.match('([^\t]+)', line).group(1)
        if seqid is None:
            seqid = locusseqid
        elif locusseqid != seqid:
            if len(prev_loci) > 0:
                yield merge_iloci(prev_loci)
                prev_loci = []
            seqid = locusseqid

        if ';gene=' in line:
            prev_loci.append(line)
            continue
        else:
            if len(prev_loci) > 0:
                yield merge_iloci(prev_loci)
                prev_loci = []
            line = re.sub('ID=[^;\n]+;*', 'geneless=true;', line)
            line = re.sub('Name=[^;\n]+;*', '', line)
            yield line
    if len(prev_loci) > 0:
        yield merge_iloci(prev_loci)


if __name__ == '__main__':
    for locus in parse_iloci(sys.stdin):
        print(locus)
