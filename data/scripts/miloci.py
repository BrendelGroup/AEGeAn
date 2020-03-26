#!/usr/bin/env python

# Copyright (c) 2010-2016, Daniel S. Standage and CONTRIBUTORS
#
# The AEGeAn Toolkit is distributed under the ISC License. See
# the 'LICENSE' file in the AEGeAn source code distribution or
# online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

from __future__ import print_function
import re
import sys


class Locus(object):
    def __init__(self, line):
        self._rawdata = line
        self.fields = line.strip().split('\t')
        assert len(self.fields) == 9

    @property
    def seqid(self):
        return self.fields[0]

    @property
    def start(self):
        return int(self.fields[3])

    @property
    def end(self):
        return int(self.fields[4])

    @property
    def ilocus_class(self):
        typematch = re.search('iLocus_type=([^;\n]+)', self.fields[8])
        assert typematch, 'could not determine iLocus type: ' + self._rawdata
        return typematch.group(1)

    @property
    def mergeable(self):
        if self.ilocus_class not in ['siLocus', 'niLocus']:
            return False
        if 'iiLocus_exception=intron-gene' in self.fields[8]:
            return False
        return True

    def __len__(self):
        return self.end - self.start + 1

    def __str__(self):
        return '\t'.join(self.fields)

    def strip(self):
        self.fields[8] = re.sub('ID=[^;\n]+;*', '', self.fields[8])
        self.fields[8] = re.sub('Name=[^;\n]+;*', '', self.fields[8])


def merge_iloci(loci, parts):
    """Merge adjacent or overlapping gene-containing iLoci."""
    assert len(loci) > 0
    if len(loci) == 1:
        loci[0].strip()
        return loci[0]

    seqid = None
    start, end = -1, -1
    attrs = {}
    for locus in loci:
        if seqid:
            assert locus.seqid == seqid
        seqid = locus.seqid
        if start == -1 or locus.start < start:
            start = locus.start
        end = max(end, locus.end)
        numeric_attrs = re.findall(r'([^;=]+=\d+)', locus.fields[8])
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

    attrstring = 'iLocus_type=miLocus'
    annotation = ''
    for key in sorted(attrs):
        attrstring += ';%s=%d' % (key, attrs[key])
    attrstring += ';'
    i = 0
    for gene in parts:
        i += 1
        geneL = Locus(gene)
        geneN = re.findall(r'(Name=[^;]+);', geneL.fields[8])
        if len(geneN) > 0:
            mgname = "miLocusGene%d=" % i
            geneN[0] = geneN[0].replace("Name=", mgname, 1)
            annotation += "%s on %s strand from %s to %s;" % \
                (geneN[0], geneL.fields[6], geneL.fields[3], geneL.fields[4])
        else:
            mgname = "miLocusGene%d=unnamed" % i
            annotation += "%s on %s strand from %s to %s;" % \
                (mgname, geneL.fields[6], geneL.fields[3], geneL.fields[4])
    attrstring += annotation
    gff3 = [seqid, 'AEGeAn::miloci.py', 'locus', str(start), str(end),
            str(len(loci)), '.', '.', attrstring]
    line = '\t'.join(gff3)
    return Locus(line)


def parse_iloci(fp):
    """
    Input: a GFF3 file containing iLoci (LocusPocus output)
    Output: merged iLoci; gene-containing iLoci that are adjacent or
            overlapping are combined
    """
    locus_buffer = []
    parts_buffer = []
    for line in fp:
        if '\tlocus\t' not in line:
            if '\tgene\t' in line:
                parts_buffer.append(line)
            continue
        else:
            locus = Locus(line)

        if len(locus_buffer) > 0 and locus.seqid != locus_buffer[0].seqid:
            yield merge_iloci(locus_buffer, parts_buffer)
            locus_buffer = []
            parts_buffer = []

        if locus.mergeable:
            locus_buffer.append(locus)
            continue
        else:
            if len(locus_buffer) > 0:
                yield merge_iloci(locus_buffer, parts_buffer)
                locus_buffer = []
            parts_buffer = []
            locus.strip()
            yield locus

    if len(locus_buffer) > 0:
        yield merge_iloci(locus_buffer, parts_buffer)


if __name__ == '__main__':
    for locus in parse_iloci(sys.stdin):
        print(locus)
