#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015-2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""
Retrieve and format data from the *Polistes dominula* genome project.

This is for the original annotation (version r1.2) published by the Toth Lab
(http://pdomgenomeproject.github.io/). The RefSeq annotation for *P. dominula*
is handled via the refseq.py module.
"""

from __future__ import print_function
import filecmp
import gzip
import re
import subprocess
import sys
import LocusPocus 


class PdomDB(LocusPocus.genomedb.GenomeDB):

    def __init__(self, label, conf, workdir='.'):
        super(PdomDB, self).__init__(label, conf, workdir)
        assert self.config['source'] == 'pdom'

    def __repr__(self):
        return 'figshare'

    @property
    def gdnaurl(self):
        return 'https://ndownloader.figshare.com/files/3557633'

    @property
    def gff3url(self):
        return 'https://ndownloader.figshare.com/files/3558071'

    @property
    def proturl(self):
        return 'https://ndownloader.figshare.com/files/3558059'

    def format_gdna(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            # No processing required currently.
            # If any is ever needed, do it here.
            print(line, end='', file=outstream)

    def format_prot(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            # No processing required currently.
            # If any is ever needed, do it here.
            print(line, end='', file=outstream)

    def format_gff3(self, logstream=sys.stderr, debug=False):
        command = ['fidibus-format-gff3.py', '--source', 'pdom', '--outfile',
                   self.gff3file, self.gff3path]
        subprocess.check_call(command)

    def gff3_protids(self, instream):
        for line in instream:
            if '\tmRNA\t' not in line:
                continue
            namematch = re.search('Name=([^;\n]+)', line)
            assert namematch, 'cannot parse mRNA name: ' + line
            yield namematch.group(1)

    def protein_mapping(self, instream):
        locusid2name = dict()
        gene2loci = dict()
        for line in instream:
            fields = line.split('\t')
            if len(fields) != 9:
                continue
            feattype = fields[2]
            attrs = fields[8]

            if feattype == 'locus':
                idmatch = re.search('ID=([^;\n]+);.*Name=([^;\n]+)', attrs)
                if idmatch:
                    locusid = idmatch.group(1)
                    locusname = idmatch.group(2)
                    locusid2name[locusid] = locusname
            elif feattype == 'gene':
                idmatch = re.search('ID=([^;\n]+);Parent=([^;\n]+)', attrs)
                assert idmatch, \
                    'Unable to parse gene and iLocus IDs: %s' % attrs
                geneid = idmatch.group(1)
                ilocusid = idmatch.group(2)
                gene2loci[geneid] = ilocusid
            elif feattype == 'mRNA':
                pattern = 'Parent=([^;\n]+);.*Name=([^;\n]+)'
                idmatch = re.search(pattern, attrs)
                assert idmatch, \
                    'Unable to parse mRNA and gene IDs: %s' % attrs
                protid = idmatch.group(2)
                geneid = idmatch.group(1)
                locusid = gene2loci[geneid]
                locusname = locusid2name[locusid]
                yield protid, locusname


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def test_download():
    """Pdom r1.2: figshare download"""
    pdom_db = LocusPocus.test_registry.genome('Pdtl')
    assert pdom_db.gdnaurl == 'https://ndownloader.figshare.com/files/3557633'
    assert pdom_db.gff3url == 'https://ndownloader.figshare.com/files/3558071'
    assert pdom_db.proturl == 'https://ndownloader.figshare.com/files/3558059'
    assert '%r' % pdom_db == 'figshare'


def test_format():
    """Pdom r1.2: formatting task"""
    pdom_db = LocusPocus.test_registry_supp.genome('Pdtl',
                                               workdir='testdata/demo-workdir')
    pdom_db.preprocess_gdna(logstream=None)
    pdom_db.preprocess_gff3(logstream=None)
    pdom_db.preprocess_prot(logstream=None)


def test_protids():
    """Pdom r1.2: extract protein IDs from GFF3"""
    db = LocusPocus.test_registry.genome('Pdtl')
    protids = ['PdomMRNAr1.2-08518.1', 'PdomMRNAr1.2-11420.1',
               'PdomMRNAr1.2-08519.1']
    infile = 'testdata/gff3/pdom-266.gff3'
    testids = list()
    with open(infile, 'r') as instream:
        for protid in db.gff3_protids(instream):
            testids.append(protid)
    assert sorted(protids) == sorted(testids), \
        'protein ID mismatch: %r %r' % (protids, testids)


def test_protmap():
    """Pdom r1.2: extract protein-->iLocus mapping from GFF3"""
    db = LocusPocus.test_registry.genome('Pdtl')
    mapping = {'PdomMRNAr1.2-08518.1': 'PdomILC-18235',
               'PdomMRNAr1.2-11420.1': 'PdomILC-18237',
               'PdomMRNAr1.2-08519.1': 'PdomILC-18238'}
    infile = 'testdata/gff3/pdom-266-loci.gff3'
    testmap = dict()
    with open(infile, 'r') as instream:
        for protid, locid in db.protein_mapping(instream):
            testmap[protid] = locid
    assert mapping == testmap, \
        'protein mapping mismatch: %r %r' % (mapping, testmap)
