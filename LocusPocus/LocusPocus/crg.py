#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015-2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""
Retrieve and format data from CRG.

GenomeDB implementation for two data sets residing in an unofficial
distribution site at the Centre de Regulacio Genomica.
"""

from __future__ import print_function
import re
import subprocess
import sys
import LocusPocus


class CrgDB(LocusPocus.genomedb.GenomeDB):

    def __init__(self, label, conf, workdir='.'):
        super(CrgDB, self).__init__(label, conf, workdir)
        assert self.config['source'] == 'crg'
        self.specbase = 'http://wasp.crg.eu'

    def __repr__(self):
        return 'CRG'

    @property
    def gdnaurl(self):
        return '%s/%s' % (self.specbase, self.gdnafilename)

    @property
    def gff3url(self):
        return '%s/%s' % (self.specbase, self.config['annotation'])

    @property
    def proturl(self):
        return '%s/%s' % (self.specbase, self.protfilename)

    @property
    def gff3filename(self):
        return '%s.gz' % self.config['annotation']

    def format_gdna(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            if line.startswith('>'):
                line = line.replace('scaffold_', '%sScf_' % self.label)
                line = line.replace('scaffold', '%sScf_' % self.label)
            print(line, end='', file=outstream)

    def format_prot(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            # No processing required currently.
            # If any is ever needed, do it here.
            print(line, end='', file=outstream)

    def format_gff3(self, logstream=sys.stderr, debug=False):
        cmds = list()
        cmds.append('gunzip -c %s' % self.gff3path)
        cmds.append("sed 's/	transcript	/	mRNA	/'")
        cmds.append("sed 's/scaffold_/%sScf_/'" % self.label)
        cmds.append("sed 's/scaffold/%sScf_/'" % self.label)
        cmds.append('fidibus-format-gff3.py --source crg -')
        cmds.append('seq-reg.py - %s' % self.gdnafile)
        cmds.append('gt gff3 -sort -tidy -o %s -force' % self.gff3file)

        commands = ' | '.join(cmds)
        if debug:  # pragma: no cover
            print('DEBUG: running command: %s' % commands, file=logstream)
        proc = subprocess.Popen(commands, shell=True, universal_newlines=True,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        for line in stderr.split('\n'):  # pragma: no cover
            if 'has not been previously introduced' not in line and \
               'does not begin with "##gff-version"' not in line and \
               line != '':
                print(line, file=logstream)
        assert proc.returncode == 0, \
            'annot cleanup command failed: %s' % commands

    def gff3_protids(self, instream):
        protids = dict()
        for line in instream:
            if '\tCDS\t' not in line:
                continue
            idmatch = re.search(r'Target=(\S+)', line)
            assert idmatch, 'cannot parse protein_id: ' + line
            protid = idmatch.group(1)
            if protid not in protids:
                protids[protid] = True
                yield protid

    def protein_mapping(self, instream):
        locusid2name = dict()
        gene2loci = dict()
        mrna2gene = dict()
        proteins = dict()
        for line in instream:
            fields = line.split('\t')
            if len(fields) != 9:
                continue
            feattype = fields[2]
            attrs = fields[8]

            if feattype == 'locus':
                idmatch = re.search(r'ID=([^;\n]+);.*Name=([^;\n]+)', attrs)
                if idmatch:
                    locusid = idmatch.group(1)
                    locusname = idmatch.group(2)
                    locusid2name[locusid] = locusname
            elif feattype == 'gene':
                idmatch = re.search(r'ID=([^;\n]+);Parent=([^;\n]+)', attrs)
                assert idmatch, \
                    'Unable to parse gene and iLocus IDs: %s' % attrs
                geneid = idmatch.group(1)
                ilocusid = idmatch.group(2)
                gene2loci[geneid] = ilocusid
            elif feattype == 'mRNA':
                idmatch = re.search(r'ID=([^;\n]+);Parent=([^;\n]+)', attrs)
                assert idmatch, \
                    'Unable to parse mRNA and gene IDs: %s' % attrs
                mrnaid = idmatch.group(1)
                geneid = idmatch.group(2)
                mrna2gene[mrnaid] = geneid
            elif feattype == 'CDS':
                idmatch = re.search(r'Parent=([^;\n]+).*Target=(\S+)', attrs)
                assert idmatch, \
                    'Unable to parse protein and mRNA IDs: %s' % attrs
                mrnaid = idmatch.group(1)
                proteinid = idmatch.group(2)
                if proteinid not in proteins:
                    geneid = mrna2gene[mrnaid]
                    ilocusid = gene2loci[geneid]
                    ilocusname = locusid2name[ilocusid]
                    proteins[proteinid] = mrnaid
                    yield proteinid, ilocusname


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------


def test_scaffolds():
    """CRG: scaffolds download"""
    dqua_db = fidibus.test_registry.genome('Dqcr')
    testurl = 'http://wasp.crg.eu/DQUA.v01.fa.gz'
    testpath = './Dqcr/DQUA.v01.fa.gz'
    assert dqua_db.gdnaurl == testurl, \
        'scaffold URL mismatch\n%s\n%s' % (dqua_db.gdnaurl, testurl)
    assert dqua_db.gdnapath == testpath, \
        'scaffold path mismatch\n%s\n%s' % (dqua_db.gdnapath, testpath)
    assert '%r' % dqua_db == 'CRG'


def test_annot():
    """CRG: annotation download"""
    dqua_db = fidibus.test_registry.genome('Dqcr', workdir='CRG')
    testurl = 'http://wasp.crg.eu/DQUA.v01.gff3'
    testpath = 'CRG/Dqcr/DQUA.v01.gff3.gz'
    assert dqua_db.gff3url == testurl, \
        'annotation URL mismatch\n%s\n%s' % (dqua_db.gff3url, testurl)
    assert dqua_db.gff3path == testpath, \
        'annotation path mismatch\n%s\n%s' % (dqua_db.gff3path, testpath)


def test_proteins():
    """CRG: protein download"""
    dqua_db = fidibus.test_registry.genome('Dqcr', workdir='/opt/db/fidibus')
    testurl = 'http://wasp.crg.eu/DQUA.v01.pep.fa.gz'
    testpath = '/opt/db/fidibus/Dqcr/DQUA.v01.pep.fa.gz'
    assert dqua_db.proturl == testurl, \
        'protein URL mismatch\n%s\n%s' % (dqua_db.proturl, testurl)
    assert dqua_db.protpath == testpath, \
        'protein path mismatch\n%s\n%s' % (dqua_db.protpath, testpath)


def test_protids():
    """CRG: extract protein IDs from GFF3"""
    db = fidibus.test_registry.genome('Dqcr')
    protids = ['DQUA011a006022P1', 'DQUA011a006023P1', 'DQUA011a006024P1']
    infile = 'testdata/gff3/dqua-275.gff3'
    testids = list()
    with open(infile, 'r') as instream:
        for protid in db.gff3_protids(instream):
            testids.append(protid)
    assert sorted(protids) == sorted(testids), \
        'protein ID mismatch: %r %r' % (protids, testids)


def test_protmap():
    """CRG: extract protein-->iLocus mapping from GFF3"""
    db = fidibus.test_registry.genome('Dqcr')
    mapping = {'DQUA011a006022P1': 'DquaILC-14465',
               'DQUA011a006023P1': 'DquaILC-14466',
               'DQUA011a006024P1': 'DquaILC-14467'}
    infile = 'testdata/gff3/dqua-275-loci.gff3'
    testmap = dict()
    with open(infile, 'r') as instream:
        for protid, locid in db.protein_mapping(instream):
            testmap[protid] = locid
    assert mapping == testmap, \
        'protein mapping mismatch: %r %r' % (mapping, testmap)


def test_format():
    """GenomeDB task drivers"""
    db = fidibus.test_registry_supp.genome('Pccr',
                                          workdir='testdata/demo-workdir')
    db.prep(logstream=None)
