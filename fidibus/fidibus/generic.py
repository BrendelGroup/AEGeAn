#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2016   Indiana University
#
# This file is part of fidibus (http://github.com/standage/fidibus) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

"""
GenomeDB implementation for a generic genome data set.
"""

from __future__ import print_function
import os
import re
import subprocess
import sys
import fidibus


class GenericDB(fidibus.genomedb.GenomeDB):

    def __init__(self, label, conf, workdir='.'):
        super(GenericDB, self).__init__(label, conf, workdir)
        assert 'gdna' in self.config
        assert 'gff3' in self.config
        assert 'prot' in self.config
        assert self.config['source'] == 'local'

    @property
    def gdnapath(self):
        return self.config['gdna']

    @property
    def gff3path(self):
        return self.config['gff3']

    @property
    def protpath(self):
        return self.config['prot']

    def download(self, logstream=sys.stderr):
        subprocess.call(['mkdir', '-p', self.dbdir])
        if logstream is not None:  # pragma: no cover
            msg = '[GenHub: %s] checking input files' % self.config['species']
            print(msg, file=logstream)
        assert os.path.isfile(self.gdnapath), \
            'gDNA file {} does not exist'.format(self.gdnapath)
        assert os.path.isfile(self.gff3path), \
            'GFF3 file {} does not exist'.format(self.gff3apath)
        assert os.path.isfile(self.protpath), \
            'proetin file {} does not exist'.format(self.protpath)

    def format_gdna(self, instream, outstream, logstream=sys.stderr):
        subprocess.call(['mkdir', '-p', self.dbdir])
        for line in instream:
            if line.strip() == '':
                continue
            print(line, end='', file=outstream)

    def format_prot(self, instream, outstream, logstream=sys.stderr):
        for line in instream:
            if line.strip() == '':
                continue
            print(line, end='', file=outstream)

    def format_gff3(self, logstream=sys.stderr, debug=False):
        cmds = list()
        if self.gff3path.endswith('.gz'):  # pragma: no cover
            cmds.append('gunzip -c %s' % self.gff3path)
        else:
            cmds.append('cat %s' % self.gff3path)
        cmds.append('seq-reg.py - %s' % self.gdnafile)
        cmds.append('fidibus-format-gff3.py --source local -')
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
            namematch = re.search('protein_id=([^;\n]+)', line)
            if not namematch:  # pragma: no cover
                namematch = re.search('Name=([^;\n]+)', line)
            assert namematch, 'cannot parse protein ID/name/accession: ' + line
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
                idmatch = re.search('Parent=([^;\n]+)', attrs)
                protmatch = re.search('protein_id=([^;\n]+)', attrs)
                if not protmatch:  # pragma: no cover
                    protmatch = re.search('Name=([^;\n]+)', attrs)
                assert idmatch and protmatch, \
                    'Unable to parse protein and gene IDs: %s' % attrs
                protid = protmatch.group(1)
                geneid = idmatch.group(1)
                locusid = gene2loci[geneid]
                locusname = locusid2name[locusid]
                yield protid, locusname


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def test_all():
    """GenericDB: the whole enchilada"""
    config = {
        'gdna': 'testdata/fasta/generic.gdna.fa.gz',
        'gff3': 'testdata/gff3/generic.gff3',
        'prot': 'testdata/fasta/generic.prot.fa',
        'source': 'local',
        'species': 'Gnrc',
    }
    db = GenericDB('Gnrc', config, workdir='testdata/demo-workdir')
    db.download(logstream=None)
    db.prep(logstream=None)
    fidibus.iloci.prepare(db, ilcformat='{}ILC-%05lu', logstream=None)
    fidibus.proteins.prepare(db, logstream=None)
    fidibus.mrnas.prepare(db, logstream=None)
    fidibus.exons.prepare(db, logstream=None)
    fidibus.stats.compute(db, logstream=None)

    sha1 = 'f3629aedbd683dd4dcf158ac12d3549e5c9081a0'
    testsha1 = db.file_sha1('testdata/demo-workdir/Gnrc/Gnrc.iloci.tsv')
    assert testsha1 == sha1, ('generic iLocus stats checksum failed')
