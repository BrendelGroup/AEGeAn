#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""Custom Genome database implementation for *Apis mellifera* OGS 1.0."""

from __future__ import print_function
import filecmp
import gzip
import re
import subprocess
import sys
import LocusPocus


class Am10DB(LocusPocus.genomedb.GenomeDB):

    def __init__(self, label, conf, workdir='.'):
        super(Am10DB, self).__init__(label, conf, workdir)
        assert self.config['source'] == 'am10'
        self.specbase = ('http://hymenopteragenome.org/drupal/sites/'
                         'hymenopteragenome.org.beebase/files/data')

    def __repr__(self):
        return 'OGS1.0'

    @property
    def gdnaurl(self):
        return '%s/%s' % (self.specbase, self.gdnafilename)

    @property
    def gff3url(self):
        return '%s/%s' % (self.specbase, self.gff3filename)

    @property
    def proturl(self):
        return '%s/%s' % (self.specbase, self.protfilename)

    def format_gdna(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            deflinematch = re.search(r'>gnl\|[^\|]+\|(\S+)', line)
            if deflinematch:
                seqid = deflinematch.group(1)
                line = line.replace('>', '>%s ' % seqid)
            print(line, end='', file=outstream)

    def format_prot(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            # No processing required currently.
            # If any is ever needed, do it here.
            print(line, end='', file=outstream)

    def format_gff3(self, logstream=sys.stderr, debug=False):
        cmds = list()
        cmds.append('gunzip -c %s' % self.gff3path)
        cmds.append('fidibus-glean-to-gff3.py')
        cmds.append('tidygff3')
        cmds.append('fidibus-format-gff3.py --source am10 -')
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
               'illegal uppercase attribute "Shift"' not in line and \
               'has the wrong phase' not in line and \
               line != '':
                print(line, file=logstream)
        assert proc.returncode == 0, \
            'annot cleanup command failed: %s' % commands

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


def test_gdna_format():
    """Amel OGSv1.0: gDNA formatting"""
    db = LocusPocus.test_registry.genome('Am10',
                                         workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Am10/Am10.gdna.fa'
    testoutfile = 'testdata/fasta/am10-gdna-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Am10 gDNA formatting failed'

    assert repr(db) == 'OGS1.0'
    url = ('http://hymenopteragenome.org/drupal/sites/'
           'hymenopteragenome.org.beebase/files/data/Amel_2.0_scaffolds.fa.gz')
    assert db.gdnaurl == url, db.gdnaurl
    url = ('http://hymenopteragenome.org/drupal/sites/'
           'hymenopteragenome.org.beebase/files/data/amel_OGSv1.0.gff.gz')
    assert db.gff3url == url, db.gff3url
    url = ('http://hymenopteragenome.org/drupal/sites/'
           'hymenopteragenome.org.beebase/files/data/amel_OGSv1.0_pep.fa.gz')
    assert db.proturl == url, db.proturl


def test_annotation_am10():
    """Amel OGSv1.0: annotation pre-processing"""
    db = LocusPocus.test_registry.genome('Am10',
                                         workdir='testdata/demo-workdir')
    db.preprocess_gff3(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Am10/Am10.gff3'
    testfile = 'testdata/gff3/am10-format.gff3'
    assert filecmp.cmp(outfile, testfile), 'Am10 annotation formatting failed'


def test_proteins_am10():
    """Amel OGSv1.0: protein formatting"""
    db = LocusPocus.test_registry.genome('Am10',
                                         workdir='testdata/demo-workdir')
    db.preprocess_prot(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Am10/Am10.all.prot.fa'
    testoutfile = 'testdata/fasta/am10-prot-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Am10 protein formatting failed'


def test_protids():
    """Amel OGSv1.0: extract protein IDs from GFF3"""
    db = LocusPocus.test_registry.genome('Am10')
    protids = ['GB18127-PA', 'GB12533-PA', 'GB12334-PA', 'GB13832-PA',
               'GB10374-PA', 'GB11673-PA', 'GB17405-PA', 'GB15709-PA',
               'GB13454-PA', 'GB16202-PA', 'GB15597-PA', 'GB13411-PA',
               'GB18700-PA', 'GB10338-PA', 'GB12987-PA', 'GB19188-PA',
               'GB18764-PA', 'GB10796-PA', 'GB12389-PA', 'GB17508-PA',
               'GB10290-PA', 'GB11017-PA', 'GB15828-PA', 'GB12008-PA',
               'GB15518-PA', 'GB14106-PA', 'GB10878-PA', 'GB17802-PA',
               'GB17122-PA', 'GB18635-PA']
    infile = 'testdata/gff3/am10-format.gff3'
    testids = list()
    with open(infile, 'r') as instream:
        for protid in db.gff3_protids(instream):
            testids.append(protid)
    assert sorted(protids) == sorted(testids), \
        'protein ID mismatch: %r %r' % (protids, testids)


def test_protmap():
    """Amel OGSv1.0: extract protein-->iLocus mapping from GFF3"""
    db = LocusPocus.test_registry.genome('Am10')
    mapping = {'GB18571-PA': 'Am10ILC-00345', 'GB13022-PA': 'Am10ILC-00347'}
    infile = 'testdata/gff3/am10-loci.gff3'
    testmap = dict()
    with open(infile, 'r') as instream:
        for protid, locid in db.protein_mapping(instream):
            testmap[protid] = locid
    assert mapping == testmap, \
        'protein mapping mismatch: %r %r' % (mapping, testmap)
