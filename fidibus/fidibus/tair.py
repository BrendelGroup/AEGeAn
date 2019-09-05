#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""
Retrieve and format data from The Arabidopsis Information Resource.

GenomeDB implementation for data residing in TAIR. Currently tested only on
TAIR6, but presumably trivial to extend to other versions using `Att6.yml` as a
template.
"""

from __future__ import print_function
import filecmp
import gzip
import os
import re
import subprocess
import sys
import fidibus


class TairDB(fidibus.genomedb.GenomeDB):

    def __init__(self, label, conf, workdir='.'):
        super(TairDB, self).__init__(label, conf, workdir)

        assert self.config['source'] == 'tair'
        assert 'version' in self.config
        assert str(self.config['version']) == '6'
        assert 'species' in self.config
        assert self.config['species'] == 'Arabidopsis thaliana'

        self.specbase = 'ftp://ftp.arabidopsis.org/home/tair'
        self.format_gdna = self.format_fasta
        self.format_prot = self.format_fasta

    def __repr__(self):
        return 'TAIR' + self.version

    @property
    def version(self):
        return str(self.config['version'])

    @property
    def gdnafilename(self):
        return self.config['gdna_filename']

    @property
    def gff3filename(self):
        return self.config['gff3_filename']

    @property
    def protfilename(self):
        return self.config['prot_filename']

    @property
    def gdnaurl(self):
        return self.config['gdna_url']

    @property
    def gff3url(self):
        return self.config['gff3_url']

    @property
    def proturl(self):
        return self.config['prot_url']

    def format_fasta(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            # No processing required.
            # If any is ever needed again, do it here.
            print(line, end='', file=outstream)

    def format_gff3(self, logstream=sys.stderr, debug=False):
        cmds = list()
        cmds.append('gunzip -c %s' % self.gff3path)
        if 'annotfilter' in self.config:  # pragma: no cover
            excludefile = self.filter_file()
            cmds.append('grep -vf %s' % excludefile.name)
        cmds.append('sed "s/Index=/index=/"')
        cmds.append('tidygff3')
        cmds.append('fidibus-format-gff3.py --source tair -')
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
               'more than one pseudogene attribute' not in line and \
               line != '':
                print(line, file=logstream)
        assert proc.returncode == 0, \
            'annot cleanup command failed: %s' % commands
        if 'annotfilter' in self.config:  # pragma: no cover
            os.unlink(excludefile.name)

    def gff3_protids(self, instream):
        protids = dict()
        for line in instream:
            if '\tmRNA\t' not in line:
                continue
            idmatch = re.search('Name=([^;\n]+)', line)
            assert idmatch, 'cannot parse protein ID: ' + line
            protid = idmatch.group(1)
            assert protid not in protids, protid
            protids[protid] = True
            yield protid

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
                if idmatch:
                    geneid = idmatch.group(1)
                    ilocusid = idmatch.group(2)
                    gene2loci[geneid] = ilocusid
                else:  # pragma: no cover
                    print('Unable to parse gene and iLocus IDs: %s' % attrs,
                          file=sys.stderr)
            elif feattype == 'mRNA':
                idmatch = re.search('Parent=([^;\n]+);Name=([^;\n]+)', attrs)
                assert idmatch, \
                    'Unable to parse mRNA and gene IDs: %s' % attrs.rstrip()
                geneid = idmatch.group(1)
                mrnaid = idmatch.group(2)
                ilocusid = gene2loci[geneid]
                ilocusname = locusid2name[ilocusid]
                yield mrnaid, ilocusname


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------


def test_gdna_format():
    """TAIR6: gDNA pre-processing"""
    db = fidibus.test_registry.genome('Att6', workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Att6/Att6.gdna.fa'
    testoutfile = 'testdata/fasta/tair6-gdna-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'TAIR6 gDNA formatting failed'

    assert repr(db) == 'TAIR6'
    ftpbase = ('ftp://ftp.arabidopsis.org/home/tair/Sequences/'
               'whole_chromosomes/OLD/chr%d.fas')
    urls = [ftpbase % x for x in (1, 2, 3, 4, 5)]
    assert db.gdnaurl == urls, db.gdnaurl
    url = ('ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR6_genome_release/'
           'TAIR6_GFF3_genes.gff')
    assert db.gff3url == url, db.gff3url
    url = ('ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR6_genome_release/'
           'TAIR6_pep_20060907')
    assert db.proturl == url, db.proturl


def test_annot_format():
    """TAIR6: annotation pre-processing"""
    db = fidibus.test_registry.genome('Att6', workdir='testdata/demo-workdir')
    db.preprocess_gff3(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Att6/Att6.gff3'
    testfile = 'testdata/gff3/tair6-format.gff3'
    assert filecmp.cmp(outfile, testfile), 'TAIR6 annotation formatting failed'


def test_protids():
    """TAIR6: extract protein IDs from GFF3"""
    db = fidibus.test_registry.genome('Att6')
    protids = ['AT5G01010.1', 'AT5G01015.1', 'AT5G01020.1', 'AT5G01030.1',
               'AT5G01030.2', 'AT5G01040.1', 'AT5G01050.1', 'AT5G01060.1',
               'AT5G01070.1', 'AT5G01075.1']
    infile = 'testdata/gff3/tair6-format.gff3'
    testids = list()
    with open(infile, 'r') as instream:
        for protid in db.gff3_protids(instream):
            testids.append(protid)
    assert sorted(protids) == sorted(testids), \
        'protein ID mismatch: %r %r' % (protids, testids)


def test_protmap():
    """TAIR6: extract protein-->iLocus mapping from GFF3"""
    db = fidibus.test_registry.genome('Att6')
    mapping = {'AT2G01800.1': 'Att6ILC-09797',
               'AT2G01810.1': 'Att6ILC-09799',
               'AT2G01820.1': 'Att6ILC-09801',
               'AT2G01830.1': 'Att6ILC-09803',
               'AT2G01830.2': 'Att6ILC-09803',
               'AT2G01830.3': 'Att6ILC-09803'}
    infile = 'testdata/gff3/tair6-loci.gff3'
    testmap = dict()
    with open(infile, 'r') as instream:
        for protid, locid in db.protein_mapping(instream):
            testmap[protid] = locid
    assert mapping == testmap, \
        'protein mapping mismatch: %r %r' % (mapping, testmap)
