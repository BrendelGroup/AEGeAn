#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015-2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import re
import sys
import fidibus


class FeatureFormatter(object):
    """Load features from GFF3, parse and (re-)attach accession numbers."""

    def __init__(self, instream, source):
        self.instream = instream
        self.source = source

        self.id2type = dict()
        self.id2acc = dict()
        self.filters = ['\tregion\t', '\tmatch\t', '\tcDNA_match\t',
                        '##species']

    def __iter__(self):
        for line in self.instream:
            line = line.rstrip()
            if self.match_filter(line):
                continue
            if self.pseudogenic_cds(line):
                continue

            self.parse_type(line)
            line = self.parse_gene(line)
            line = self.parse_transcript(line)
            line = self.parse_vdj(line)
            line = self.parse_feature(line)

            yield line

    def parse_type(self, line):
        fields = line.split('\t')
        if len(fields) != 9:
            return False

        ftype = fields[2]
        attributes = fields[8]
        idmatch = re.search(r'ID=([^;\n]+)', attributes)
        if idmatch:
            featureid = idmatch.group(1)
            self.id2type[featureid] = ftype

    def match_filter(self, line):
        for filt in self.filters:
            if filt in line:
                return True
        return False

    def pseudogenic_cds(self, line):
        """
        Test whether the given entry is a pseudogene-associated CDS.

        We want to ignore these!
        """
        if self.source == 'tair':
            return False

        fields = line.split('\t')
        if len(fields) != 9:
            return False

        ftype = fields[2]
        attributes = fields[8]
        if ftype != 'CDS':
            return False

        parentmatch = re.search(r'Parent=([^;\n]+)', attributes)
        assert parentmatch, attributes
        parentid = parentmatch.group(1)
        if self.id2type[parentid] == 'pseudogene':
            return True
        return False

    def parse_gene(self, line):
        """Parse accession for gene features."""
        fields = line.split('\t')
        if len(fields) != 9:
            return line

        ftype = fields[2]
        attributes = fields[8]
        if ftype != 'gene' or '\tAEGeAn::tidygff3\t' in line:
            return line

        accmatch = None
        if self.source == 'refseq':
            accmatch = re.search(r'GeneID:([^;,\n]+)', line)
        elif self.source == 'crg':
            accmatch = re.search(r'ID=([^;\n]+)', line)
        elif self.source in ['genbank', 'pdom', 'tair', 'beebase']:
            accmatch = re.search(r'Name=([^;\n]+)', line)
        elif self.source == 'local':
            accmatch = re.search(r'accession=([^;\n]+)', attributes)
            if not accmatch:
                accmatch = re.search(r'Name=([^;\n]+)', attributes)
        else:
            pass
        assert accmatch, 'unable to parse gene accession: %s' % line
        accession = accmatch.group(1)

        idmatch = re.search(r'ID=([^;\n]+)', attributes)
        if idmatch:
            geneid = idmatch.group(1)
            self.id2acc[geneid] = accession
        else:
            print('Warning: gene has no ID: %s' % attributes, file=sys.stderr)

        if 'accession=' in line:
            return line
        return line + ';accession=' + accession

    def parse_transcript(self, line):
        """Parse accession for transcript features."""
        fields = line.split('\t')
        if len(fields) != 9:
            return line

        ftype = fields[2]
        attributes = fields[8]
        ttypes = [
            'mRNA', 'tRNA', 'rRNA', 'transcript', 'primary_transcript',
            'ncRNA', 'miRNA', 'snRNA', 'snoRNA', 'lnc_RNA', 'scRNA', 'SRP_RNA',
            'antisense_RNA', 'RNase_P_RNA', 'telomerase_RNA', 'piRNA',
            'RNase_MRP_RNA', 'guide_RNA',
        ]
        if ftype not in ttypes:
            return line

        accmatch = None
        idmatch = None
        parentaccession = None
        if self.source == 'refseq':
            accmatch = re.search(r'transcript_id=([^;\n]+)', attributes)
            idmatch = re.search(r'GeneID:([^;,\n]+)', attributes)
        elif self.source == 'genbank':
            parentid = re.search(r'Parent=([^;\n]+)', attributes).group(1)
            parentaccession = self.id2acc[parentid]
        elif self.source in ['crg', 'pdom']:
            accmatch = re.search(r'ID=([^;\n]+)', attributes)
        elif self.source in ['beebase', 'tair', 'am10']:
            accmatch = re.search(r'Name=([^;\n]+)', attributes)
        elif self.source == 'local':
            accmatch = re.search(r'protein_id=([^;\n]+)', attributes)
            if not accmatch:
                accmatch = re.search(r'accession=([^;\n]+)', attributes)
            if not accmatch:
                accmatch = re.search(r'Name=([^;\n]+)', attributes)
        else:
            pass
        assert accmatch or idmatch or parentaccession, \
            'unable to parse transcript accession: %s' % line
        if accmatch:
            accession = accmatch.group(1)
        elif parentaccession:
            accession = '{}.{}'.format(parentaccession, ftype)
        else:
            accession = '%s:%s' % (idmatch.group(1), ftype)

        rnaidmatch = re.search(r'ID=([^;\n]+)', attributes)
        if rnaidmatch:
            rnaid = rnaidmatch.group(1)
            self.id2acc[rnaid] = accession
        else:
            print('Warning: RNA has no ID: %s' % attributes, file=sys.stderr)

        if 'accession=' in line:
            return line
        return line + ';accession=' + accession

    def parse_vdj(self, line):
        """Parse accessions for features of V(D)J genes."""
        fields = line.split('\t')
        if len(fields) != 9:
            return line

        ftype = fields[2]
        if ftype not in ['V_gene_segment', 'D_gene_segment', 'J_gene_segment',
                         'C_gene_segment']:
            return line

        vdjid = re.search(r'ID=([^;\n]+)', line).group(1)
        accession = re.search(r'GeneID:([^;,\n]+)', line).group(1)
        self.id2acc[vdjid] = accession

        return line + ';accession=' + accession

    def parse_feature(self, line):
        """Parse accession for exons, introns, and coding sequences"""
        if 'accession=' in line:
            return line

        fields = line.split('\t')
        if len(fields) != 9:
            return line

        ftype = fields[2]
        if ftype not in ['exon', 'intron', 'CDS']:
            return line

        parentid = re.search(r'Parent=([^;\n]+)', line).group(1)
        if self.source == 'tair':
            for pid in parentid.split(','):
                if 'RNA' in pid:
                    parentid = pid
        assert ',' not in parentid, parentid
        assert parentid in self.id2acc, parentid
        accession = self.id2acc[parentid]
        return line + ';accession=' + accession


def parse_args():
    """Define the command-line interface."""
    desc = 'Filter features and parse accession values'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-v', '--version', action='version',
                        version='fidibus v%s' % fidibus.__version__)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-p', '--prefix', default=None, metavar='PFX',
                        help='attach the given prefix to each sequence ID')
    parser.add_argument('--source', default='refseq', choices=fidibus.sources,
                        help='data source; default is "refseq"')
    parser.add_argument('gff3', type=argparse.FileType('r'))
    return parser.parse_args()


def format_prefix(line, prefix):
    """Apply a prefix to each sequence ID."""
    if len(line.split('\t')) == 9:
        return prefix + line
    elif line.startswith('##sequence-region'):
        return re.sub(r'##sequence-region(\s+)(\S+)',
                      r'##sequence-region\g<1>%s\g<2>' % args.prefix, line)


def main():
    args = parse_args()
    formatter = FeatureFormatter(args.gff3, args.source)
    for line in formatter:
        if args.prefix:
            line = format_prefix(line, args.prefix)
        print(line, file=args.outfile)


if __name__ == '__main__':
    main()
