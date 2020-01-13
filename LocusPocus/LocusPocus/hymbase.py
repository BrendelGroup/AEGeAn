#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017   Regents of the University of California.
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""
Retrieve and format data from BeeBase.

GenomeDB implementation for BeeBase consortium data provisioned by
HymenopteraBase.
"""

from __future__ import print_function
import filecmp
import gzip
import re
import subprocess
import sys
import LocusPocus


class HymBaseDB(LocusPocus.genomedb.GenomeDB):

    def specbase(self):
        genus = self.config['species'].split()[0].lower()
        base = 'http://hymenopteragenome.org/drupal/sites/'
        base += 'hymenopteragenome.org.%s/files/data' % genus
        return base

    def __init__(self, label, conf, workdir='.'):
        super(HymBaseDB, self).__init__(label, conf, workdir)
        assert self.config['source'] in ['hymbase', 'beebase']

    def __repr__(self):
        return 'HymenopteraBase'

    @property
    def gdnaurl(self):
        return '%s/%s' % (self.specbase(), self.gdnafilename)

    @property
    def gff3url(self):
        return '%s/%s' % (self.specbase(), self.gff3filename)

    @property
    def proturl(self):
        return '%s/%s' % (self.specbase(), self.protfilename)

    def format_gdna(self, instream, outstream, logstream=sys.stderr):
        self.format_fasta(instream, outstream, logstream)

    def format_prot(self, instream, outstream, logstream=sys.stderr):
        self.format_fasta(instream, outstream, logstream)

    def format_fasta(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            if line.startswith('>gnl|'):
                deflinematch = re.search(r'>gnl\|[^\|]+\|(\S+)', line)
                assert deflinematch, line
                protid = deflinematch.group(1)
                line = line.replace('>', '>%s ' % protid)
            print(line, end='', file=outstream)

    def format_gff3(self, logstream=sys.stderr, debug=False):
        cmds = list()
        cmds.append('gunzip -c %s' % self.gff3path)
        if 'annotfilter' in self.config:
            excludefile = self.filter_file()
            cmds.append('grep -vf %s' % excludefile.name)
        cmds.append("grep -v '\tregion\t'")
        cmds.append('fidibus-uniq.py')
        cmds.append('fidibus-namedup.py')
        cmds.append('tidygff3')
        cmds.append('fidibus-format-gff3.py --source beebase -')
        cmds.append('seq-reg.py - %s' % self.gdnafile)
        cmds.append('gt gff3 -sort -tidy -o %s -force' % self.gff3file)

        commands = 'bash -o pipefail -c "%s"' % ' | '.join(cmds)
        if debug:  # pragma: no cover
            print('DEBUG: running command: %s' % commands, file=logstream)
        proc = subprocess.Popen(commands, shell=True, universal_newlines=True,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        for line in stderr.split('\n'):  # pragma: no cover
            if 'has not been previously introduced' not in line and \
               'does not begin with "##gff-version"' not in line and \
               'has the wrong phase' not in line and \
               line != '':
                print(line, file=logstream)
        assert proc.returncode == 0, \
            'annot cleanup command failed: %s' % commands

    def gff3_protids(self, instream):
        for line in instream:
            if '\tmRNA\t' not in line:
                continue
            namematch = re.search(r'Name=([^;\n]+)', line)
            assert namematch, 'cannot parse mRNA name: ' + line
            protid = namematch.group(1).replace('-RA', '-PA')
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
                pattern = 'Parent=([^;\n]+);.*Name=([^;\n]+)'
                idmatch = re.search(pattern, attrs)
                assert idmatch, \
                    'Unable to parse mRNA and gene IDs: %s' % attrs
                protid = idmatch.group(2).replace('-RA', '-PA')
                geneid = idmatch.group(1)
                locusid = gene2loci[geneid]
                locusname = locusid2name[locusid]
                yield protid, locusname


class BeeBaseDB(HymBaseDB):

    def specbase(self):
        return ('http://hymenopteragenome.org/beebase/sites/'
                'hymenopteragenome.org.beebase/files/data/consortium_data')

    def __init__(self, label, conf, workdir='.'):
        super(BeeBaseDB, self).__init__(label, conf, workdir)
        assert self.config['source'] == 'beebase'

    def __repr__(self):
        return 'BeeBase'

    def format_gdna(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            if line.startswith('>scaffold'):
                line = line.replace('scaffold', '%sScf_' % self.label)
            elif line.startswith('>Group'):
                line = line.replace('Group', '%sGroup' % self.label)
            print(line, end='', file=outstream)

    def format_prot(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            if line.startswith('>gnl|'):
                deflinematch = re.search(r'>gnl\|[^\|]+\|(\S+)', line)
                assert deflinematch, line
                protid = deflinematch.group(1)
                line = line.replace('>', '>%s ' % protid)
            print(line, end='', file=outstream)

    def format_gff3(self, logstream=sys.stderr, debug=False):
        cmds = list()
        cmds.append('gunzip -c %s' % self.gff3path)
        cmds.append('fidibus-namedup.py')
        cmds.append("sed 's/scaffold/%sScf_/'" % self.label)
        cmds.append("sed 's/Group/%sGroup/'" % self.label)
        cmds.append('tidygff3')
        cmds.append('fidibus-format-gff3.py --source beebase -')
        cmds.append('seq-reg.py - %s' % self.gdnafile)
        cmds.append('gt gff3 -sort -tidy -o %s -force' % self.gff3file)

        commands = 'bash -o pipefail -c "%s"' % ' | '.join(cmds)
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


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------


def test_scaffolds_download():
    """HymBase: scaffolds download"""
    cfdb = fidibus.test_registry.genome('Cfhb')
    testurl = ('http://hymenopteragenome.org/drupal/sites/hymenopteragenome.'
               'org.camponotus/files/data/Cflo_3.3_scaffolds.fa.gz')
    testpath = './Cfhb/Cflo_3.3_scaffolds.fa.gz'
    assert cfdb.gdnaurl == testurl, \
        'scaffold URL mismatch\n%s\n%s' % (cfdb.gdnaurl, testurl)
    assert cfdb.gdnapath == testpath, \
        'scaffold path mismatch\n%s\n%s' % (cfdb.gdnapath, testpath)
    assert '%r' % cfdb == 'HymenopteraBase'


def test_annotation_hymbase():
    """HymBase: annotation pre-processing"""
    db = fidibus.test_registry.genome('Cfhb', workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    db.preprocess_gff3(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Cfhb/Cfhb.gff3'
    testfile = 'testdata/gff3/hymbase-format-cfhb.gff3'
    assert filecmp.cmp(outfile, testfile), 'Cfhb annotation formatting failed'

    db = fidibus.test_registry.genome('Sihb', workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    db.preprocess_gff3(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Sihb/Sihb.gff3'
    testfile = 'testdata/gff3/hymbase-format-sihb.gff3'
    assert filecmp.cmp(outfile, testfile), 'Sihb annotation formatting failed'


def test_gdna_format_hymbase():
    """HymBase: gDNA pre-processing"""
    db = fidibus.test_registry.genome('Sihb', workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Sihb/Sihb.gdna.fa'
    testoutfile = 'testdata/fasta/hymbase-format-sihb.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Sihb gDNA formatting failed'


def test_proteins_hymbase():
    """HymBase: protein pre-processing"""
    db = fidibus.test_registry.genome('Sihb', workdir='testdata/demo-workdir')
    db.preprocess_prot(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Sihb/Sihb.all.prot.fa'
    testoutfile = 'testdata/fasta/hymbase-format-Sihb-prot.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Sihb protein formatting failed'


def test_scaffolds_download_beebase():
    """BeeBase: scaffolds download"""
    emex_db = fidibus.test_registry.genome('Emex')
    testurl = ('http://hymenopteragenome.org/beebase/sites/'
               'hymenopteragenome.org.beebase/files/data/consortium_data/'
               'Eufriesea_mexicana.v1.0.fa.gz')
    testpath = './Emex/Eufriesea_mexicana.v1.0.fa.gz'
    assert emex_db.gdnaurl == testurl, \
        'scaffold URL mismatch\n%s\n%s' % (emex_db.gdnaurl, testurl)
    assert emex_db.gdnapath == testpath, \
        'scaffold path mismatch\n%s\n%s' % (emex_db.gdnapath, testpath)
    assert '%r' % emex_db == 'BeeBase'


def test_annot_download():
    """BeeBase: annotation download"""
    dnov_db = fidibus.test_registry.genome('Dnov', workdir='BeeBase')
    testurl = ('http://hymenopteragenome.org/beebase/sites/'
               'hymenopteragenome.org.beebase/files/data/consortium_data/'
               'Dufourea_novaeangliae_v1.1.gff.gz')
    testpath = 'BeeBase/Dnov/Dufourea_novaeangliae_v1.1.gff.gz'
    assert dnov_db.gff3url == testurl, \
        'annotation URL mismatch\n%s\n%s' % (dnov_db.gff3url, testurl)
    assert dnov_db.gff3path == testpath, \
        'annotation path mismatch\n%s\n%s' % (dnov_db.gff3path, testpath)


def test_proteins_download():
    """BeeBase: protein download"""
    hlab_db = fidibus.test_registry.genome('Hlab', workdir='/opt/db/fidibus')
    testurl = ('http://hymenopteragenome.org/beebase/sites/'
               'hymenopteragenome.org.beebase/files/data/consortium_data/'
               'Habropoda_laboriosa_v1.2.pep.fa.gz')
    testpath = '/opt/db/fidibus/Hlab/Habropoda_laboriosa_v1.2.pep.fa.gz'
    assert hlab_db.proturl == testurl, \
        'protein URL mismatch\n%s\n%s' % (hlab_db.proturl, testurl)
    assert hlab_db.protpath == testpath, \
        'protein path mismatch\n%s\n%s' % (hlab_db.protpath, testpath)


def test_gdna_format():
    """BeeBase: gDNA pre-processing"""
    db = fidibus.test_registry.genome('Hlab', workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Hlab/Hlab.gdna.fa'
    testoutfile = 'testdata/fasta/hlab-first-6-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Hlab gDNA formatting failed'

    db = fidibus.test_registry.genome('Am32', workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Am32/Am32.gdna.fa'
    testoutfile = 'testdata/fasta/am32-gdna-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Am32 gDNA formatting failed'


def test_annotation_beebase():
    """BeeBase: annotation pre-processing"""
    db = fidibus.test_registry.genome('Hlab', workdir='testdata/demo-workdir')
    db.preprocess_gff3(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Hlab/Hlab.gff3'
    testfile = 'testdata/gff3/beebase-format-hlab.gff3'
    assert filecmp.cmp(outfile, testfile), 'Hlab annotation formatting failed'


def test_proteins_beebase():
    """BeeBase: protein pre-processing"""
    db = fidibus.test_registry.genome('Hlab', workdir='testdata/demo-workdir')
    db.preprocess_prot(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Hlab/Hlab.all.prot.fa'
    testoutfile = 'testdata/fasta/hlab-first-20-prot-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Hlab protein formatting failed'

    db = fidibus.test_registry.genome('Am32', workdir='testdata/demo-workdir')
    db.preprocess_prot(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Am32/Am32.all.prot.fa'
    testoutfile = 'testdata/fasta/am32-prot-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Am32 protein formatting failed'


def test_protids():
    """BeeBase: extract protein IDs from GFF3"""
    db = fidibus.test_registry.genome('Hlab')
    protids = ['Hlab050%d' % x for x in range(62, 75)]
    infile = 'testdata/gff3/hlab-238.gff3'
    testids = list()
    with open(infile, 'r') as instream:
        for protid in db.gff3_protids(instream):
            testids.append(protid)
    assert sorted(protids) == sorted(testids), \
        'protein ID mismatch: %r %r' % (protids, testids)


def test_protmap():
    """BeeBase: extract protein-->iLocus mapping from GFF3"""
    db = fidibus.test_registry.genome('Hlab')
    mapping = {'Hlab05074': 'HlabILC-653102', 'Hlab05071': 'HlabILC-653096',
               'Hlab05062': 'HlabILC-653079', 'Hlab05068': 'HlabILC-653090',
               'Hlab05072': 'HlabILC-653098', 'Hlab05070': 'HlabILC-653094',
               'Hlab05066': 'HlabILC-653086', 'Hlab05067': 'HlabILC-653088',
               'Hlab05069': 'HlabILC-653092', 'Hlab05073': 'HlabILC-653100',
               'Hlab05064': 'HlabILC-653083', 'Hlab05065': 'HlabILC-653085',
               'Hlab05063': 'HlabILC-653081'}
    infile = 'testdata/gff3/hlab-238-loci.gff3'
    testmap = dict()
    with open(infile, 'r') as instream:
        for protid, locid in db.protein_mapping(instream):
            testmap[protid] = locid
    assert mapping == testmap, \
        'protein mapping mismatch: %r %r' % (mapping, testmap)
