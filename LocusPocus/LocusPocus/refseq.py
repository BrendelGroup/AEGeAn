#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015-2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""
Retrieve and format data from RefSeq.

GenomeDB implementation for data residing in NCBI's RefSeq database.
"""

from __future__ import print_function
import filecmp
import gzip
import os
import re
import subprocess
import sys
import LocusPocus
from LocusPocus import test_registry


class RefSeqDB(LocusPocus.genomedb.GenomeDB):

    @classmethod
    def base(self):
        return 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq'

    def __init__(self, label, conf, workdir='.'):
        super(RefSeqDB, self).__init__(label, conf, workdir)

        assert self.config['source'] in ['refseq', 'genbank']
        assert 'branch' in self.config
        assert 'species' in self.config
        assert 'accession' in self.config
        assert 'build' in self.config

        species = self.config['species'].replace(' ', '_')
        self.acc = self.config['accession'] + '_' + self.config['build']

        pathbase = [self.config['branch'], species, 'all_assembly_versions']
        if 'suppressed' in self.config and self.config['suppressed'] is True:
            pathbase.append('suppressed')
        pathbase.extend([self.acc, self.acc])
        self.specbase = '/'.join(pathbase)
        self.urlbase = '/'.join([self.base()] + pathbase)
        self.format_gdna = self.format_fasta
        self.format_prot = self.format_fasta

    def __repr__(self):
        return 'RefSeq'

    @property
    def gdnafilename(self):
        return '%s_genomic.fna.gz' % self.acc

    @property
    def gff3filename(self):
        return '%s_genomic.gff.gz' % self.acc

    @property
    def protfilename(self):
        return '%s_protein.faa.gz' % self.acc

    @property
    def gdnaurl(self):
        return '%s_genomic.fna.gz' % self.urlbase

    @property
    def gff3url(self):
        return '%s_genomic.gff.gz' % self.urlbase

    @property
    def proturl(self):
        return '%s_protein.faa.gz' % self.urlbase

    def download(self, localpath=None, logstream=sys.stderr):  # pragma: no cover  # noqa
        """Override the download task to enable connection to local mirrors."""
        if localpath is None:
            super(RefSeqDB, self).download(logstream=logstream)
        else:
            subprocess.call(['mkdir', '-p', self.dbdir])
            asmblpath = '/'.join(localpath, specbase, '_genomic.fna.gz')
            annotpath = '/'.join(localpath, specbase, '_genomic.gff.gz')
            protpath = '/'.join(localpath, specbase, '_protein.faa.gz')
            copy_file(asmblpath, self.gdnafilename)
            copy_file(annotpath, self.gff3filename)
            copy_file(protpath, self.protfilename)

    def format_fasta(self, instream, outstrm, logstream=sys.stderr):
        for defline, sequence in LocusPocus.fasta.parse(instream):
            if 'seqfilter' in self.config:
                discard = False
                for pattern in self.config['seqfilter']:
                    if pattern in defline:
                        discard = True
                        break
                if discard:
                    continue
            print(defline, file=outstrm)
            LocusPocus.fasta.format(sequence, linewidth=80, outstream=outstrm)

    def format_gff3(self, logstream=sys.stderr, debug=False):
        cmds = list()
        cmds.append('gunzip -c %s' % self.gff3path)
        if 'annotfilter' in self.config:
            excludefile = self.filter_file()
            cmds.append('grep -vf %s' % excludefile.name)
        cmds.append('tidygff3')
        cmds.append('fidibus-format-gff3.py --source %s -' % str(self).lower())
        if 'fixseqreg' in self.config and self.config['fixseqreg'] is True:
            cmds.append('seq-reg.py - %s' % self.gdnafile)  # pragma: no cover
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
               'more than one pseudogene attribute' not in line and \
               line != '':
                print(line, file=logstream)
        assert proc.returncode == 0, \
            'annot cleanup command failed: %s' % commands
        if 'annotfilter' in self.config:
            os.unlink(excludefile.name)

    def gff3_protids(self, instream):
        protids = dict()
        for line in instream:
            if '\tCDS\t' not in line:
                continue
            idmatch = re.search('protein_id=([^;\n]+)', line)
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
        protdups = dict()
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
                idmatch = re.search('ID=([^;\n]+);Parent=([^;\n]+)', attrs)
                assert idmatch, \
                    'Unable to parse mRNA and gene IDs: %s' % attrs.rstrip()
                mrnaid = idmatch.group(1)
                geneid = idmatch.group(2)
                mrna2gene[mrnaid] = geneid
            elif feattype == 'CDS':
                if 'exception=rearrangement required for product' in attrs:
                    continue
                idmatch = re.search('Parent=([^;\n]+).*protein_id=([^;\n]+)',
                                    attrs)
                assert idmatch, \
                    'Unable to parse protein and mRNA IDs: %s' % attrs
                mrnaid = idmatch.group(1)
                proteinid = idmatch.group(2)
                if proteinid not in proteins:
                    geneid = mrna2gene[mrnaid]
                    proteins[proteinid] = mrnaid
                    ilocusid = gene2loci[geneid]
                    ilocusname = locusid2name[ilocusid]
                    yield proteinid, ilocusname


class GenbankDB(RefSeqDB):

    @classmethod
    def base(self):
        return 'https://ftp.ncbi.nlm.nih.gov/genomes/genbank'

    def __repr__(self):
        return 'Genbank'


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------


def test_genome_download():
    """RefSeq: gDNA download"""
    ador_db = test_registry.genome('Ador')
    testurl = ('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/'
               'Apis_dorsata/all_assembly_versions/'
               'GCF_000469605.1_Apis_dorsata_1.3/'
               'GCF_000469605.1_Apis_dorsata_1.3_genomic.fna.gz')
    testpath = './Ador/GCF_000469605.1_Apis_dorsata_1.3_genomic.fna.gz'
    assert '%r' % ador_db == 'RefSeq'
    assert ador_db.gdnaurl == testurl, \
        'scaffold URL mismatch\n%s\n%s' % (ador_db.gdnaurl, testurl)
    assert ador_db.gdnapath == testpath, \
        'scaffold path mismatch\n%s\n%s' % (ador_db.gdnapath, testpath)
    assert ador_db.compress_gdna is False

    amel_db = test_registry.genome('Amel')
    testurl = ('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/'
               'Apis_mellifera/all_assembly_versions/suppressed/'
               'GCF_000002195.4_Amel_4.5/'
               'GCF_000002195.4_Amel_4.5_genomic.fna.gz')
    testpath = './Amel/GCF_000002195.4_Amel_4.5_genomic.fna.gz'
    assert amel_db.gdnaurl == testurl, \
        'chromosome URL mismatch\n%s\n%s' % (amel_db.gdnaurl, testurl)
    assert amel_db.gdnapath == testpath, \
        'chromosome path mismatch\n%s\n%s' % (amel_db.gdnapath, testpath)
    assert ador_db.compress_gdna is False


def test_annot_download():
    """RefSeq: annotation download"""
    ador_db = test_registry.genome('Ador')
    testurl = ('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/'
               'Apis_dorsata/all_assembly_versions/'
               'GCF_000469605.1_Apis_dorsata_1.3/'
               'GCF_000469605.1_Apis_dorsata_1.3_genomic.gff.gz')
    testpath = './Ador/GCF_000469605.1_Apis_dorsata_1.3_genomic.gff.gz'
    assert ador_db.gff3url == testurl, \
        'annotation URL mismatch\n%s\n%s' % (ador_db.gff3url, testurl)
    assert ador_db.gff3path == testpath, \
        'annotation path mismatch\n%s\n%s' % (ador_db.gff3path, testpath)
    assert ador_db.compress_gff3 is False


def test_proteins_download():
    """RefSeq: protein download"""
    db = test_registry.genome('Ador', workdir='/home/gandalf/HymHub')
    testurl = ('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/'
               'Apis_dorsata/all_assembly_versions/'
               'GCF_000469605.1_Apis_dorsata_1.3/'
               'GCF_000469605.1_Apis_dorsata_1.3_protein.faa.gz')
    testpath = ('/home/gandalf/HymHub/Ador/'
                'GCF_000469605.1_Apis_dorsata_1.3_protein.faa.gz')
    assert db.proturl == testurl, \
        'protein URL mismatch\n%s\n%s' % (db.proturl, testurl)
    assert db.protpath == testpath, \
        'protein path mismatch\n%s\n%s' % (db.protpath, testpath)
    assert db.compress_prot is False


def test_gdna_format():
    """RefSeq: gDNA pre-processing"""
    db = test_registry.genome('Hsal', workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Hsal/Hsal.gdna.fa'
    testoutfile = 'testdata/fasta/hsal-first-7-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Hsal gDNA formatting failed'

    conf = test_registry.genome('Tcas')
    db = test_registry.genome('Tcas', workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Tcas/Tcas.gdna.fa'
    testoutfile = 'testdata/fasta/tcas-first-33-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Tcas gDNA formatting failed'

    conf = test_registry.genome('Mmus')
    db = test_registry.genome('Mmus', workdir='testdata/demo-workdir')
    db.preprocess_gdna(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Mmus/Mmus.gdna.fa'
    testoutfile = 'testdata/fasta/mmus-gdna.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Mmus gDNA formatting failed'


def test_annot_format():
    """RefSeq: annotation pre-processing"""
    db = test_registry.genome('Aech', workdir='testdata/demo-workdir')
    db.preprocess_gff3(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Aech/Aech.gff3'
    testfile = 'testdata/gff3/ncbi-format-aech.gff3'
    assert filecmp.cmp(outfile, testfile), 'Aech annotation formatting failed'

    db = test_registry.genome('Pbar', workdir='testdata/demo-workdir')
    db.config['annotfilter'] = 'NW_011933506.1'
    db.preprocess_gff3(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Pbar/Pbar.gff3'
    testfile = 'testdata/gff3/ncbi-format-pbar.gff3'
    assert filecmp.cmp(outfile, testfile), 'Pbar annotation formatting failed'

    db = test_registry.genome('Ador', workdir='testdata/demo-workdir')
    db.config['annotfilter'] = ['NW_006264094.1', 'NW_006263516.1']
    db.preprocess_gff3(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Ador/Ador.gff3'
    testfile = 'testdata/gff3/ncbi-format-ador.gff3'
    assert filecmp.cmp(outfile, testfile), 'Ador annotation formatting failed'

    db.config['checksums']['gff3'] = 'b0gU$h@sH'
    passed = False
    try:
        db.preprocess_gff3(logstream=None, verify=True)
    except Exception as e:
        if 'integrity check failed' in str(e):
            passed = True
        else:  # pragma: no cover
            raise e
    assert passed is True, 'Ador bogus hash validated'

    db.preprocess_gff3(logstream=None, verify=True, strict=False)


def test_prot_ncbi():
    """RefSeq: protein pre-processing"""
    db = test_registry.genome('Hsal', workdir='testdata/demo-workdir')
    db.preprocess_prot(logstream=None, verify=False)
    outfile = 'testdata/demo-workdir/Hsal/Hsal.all.prot.fa'
    testoutfile = 'testdata/fasta/hsal-13-prot-out.fa'
    assert filecmp.cmp(testoutfile, outfile), 'Hsal protein formatting failed'


def test_protids():
    """RefSeq: extract protein IDs from GFF3"""
    db = test_registry.genome('Xtro')
    protids = ['XP_012809995.1', 'XP_012809996.1', 'XP_012809997.1',
               'XP_012809998.1']
    infile = 'testdata/gff3/xtro-3genes.gff3'
    testids = list()
    with open(infile, 'r') as instream:
        for protid in db.gff3_protids(instream):
            testids.append(protid)
    assert sorted(protids) == sorted(testids), \
        'protein ID mismatch: %r %r' % (protids, testids)


def test_protmap():
    """RefSeq: extract protein-->iLocus mapping from GFF3"""
    db = test_registry.genome('Xtro')
    mapping = {'XP_012809997.1': 'XtroILC-43374',
               'XP_012809996.1': 'XtroILC-43374',
               'XP_012809995.1': 'XtroILC-43373',
               'XP_012809998.1': 'XtroILC-43374'}
    infile = 'testdata/gff3/xtro-3genes-loci.gff3'
    testmap = dict()
    with open(infile, 'r') as instream:
        for protid, locid in db.protein_mapping(instream):
            testmap[protid] = locid
    assert mapping == testmap, \
        'protein mapping mismatch: %r %r' % (mapping, testmap)


def test_cleanup():
    """RefSeq: cleanup task"""
    db = test_registry.genome('Aech', workdir='testdata/demo-workdir')
    delfiles = ['testdata/demo-workdir/Aech/Aech.gff3']
    testfiles = db.cleanup(None, False, True)
    assert testfiles == delfiles, '%r %r' % (testfiles, delfiles)
    delfiles = ['testdata/demo-workdir/Aech/Aech.gff3',
                'testdata/demo-workdir/Aech/GCF_000204515.1_Aech_3.9_genomic.'
                'gff.gz']
    testfiles = db.cleanup(None, True, True)
    assert set(testfiles) == set(delfiles), '%r %r' % (testfiles, delfiles)

    db = test_registry.genome('Bdis', workdir='testdata/demo-workdir')
    delfiles = ['testdata/demo-workdir/Bdis/Bdis.gdna.fa',
                'testdata/demo-workdir/Bdis/Bdis.gff3',
                'testdata/demo-workdir/Bdis/Bdis.ilocus.mrnas.gff3',
                'testdata/demo-workdir/Bdis/Bdis.miloci.fa',
                'testdata/demo-workdir/Bdis/Bdis.mrnas.txt',
                'testdata/demo-workdir/Bdis/Bdis.simple-iloci.txt',
                'testdata/demo-workdir/Bdis/ilens.temp']
    testfiles = db.cleanup(None, False, True)
    assert set(testfiles) == set(delfiles), '%r %r' % (testfiles, delfiles)
    delfiles = ['testdata/demo-workdir/Bdis/Bdis.gdna.fa',
                'testdata/demo-workdir/Bdis/Bdis.gff3',
                'testdata/demo-workdir/Bdis/Bdis.ilocus.mrnas.gff3',
                'testdata/demo-workdir/Bdis/Bdis.mrnas.txt',
                'testdata/demo-workdir/Bdis/Bdis.simple-iloci.txt',
                'testdata/demo-workdir/Bdis/ilens.temp']
    testfiles = db.cleanup(['.miloci.'], False, True)
    assert set(testfiles) == set(delfiles), '%r %r' % (testfiles, delfiles)
    delfiles = ['testdata/demo-workdir/Bdis/Bdis.gdna.fa',
                'testdata/demo-workdir/Bdis/Bdis.gff3',
                'testdata/demo-workdir/Bdis/Bdis.ilocus.mrnas.gff3',
                'testdata/demo-workdir/Bdis/Bdis.mrnas.txt',
                'testdata/demo-workdir/Bdis/ilens.temp']
    testfiles = db.cleanup(['.miloci.', 'simple'], False, True)
    assert set(testfiles) == set(delfiles), '%r %r' % (testfiles, delfiles)

    db = test_registry.genome('Vcar', workdir='testdata/demo-workdir')
    nodelfile = 'testdata/demo-workdir/Vcar/Vcar.protein2ilocus.tsv'
    testfiles = db.cleanup(None, False, True)
    assert nodelfile not in testfiles, 'incorrectly deleted file'


def test_get_map():
    """RefSeq: get prot map"""
    db = test_registry.genome('Vcar', workdir='testdata/demo-workdir')
    mapping = dict()
    for protid, locid in db.get_prot_map():
        mapping[protid] = locid
    testmap = {'XP_002945607.1': 'VcarILC-00002',
               'XP_002945976.1': 'VcarILC-00003',
               'XP_002945608.1': 'VcarILC-00004',
               'XP_002945609.1': 'VcarILC-00005'}
    assert mapping == testmap, 'get prot map fail: %r %r' % (testmap, mapping)


def test_genbank():
    """Genbank: smoke test"""
    db = test_registry.genome('Znev')
    assert str(db) == 'Genbank'
    assert db.base() == 'https://ftp.ncbi.nlm.nih.gov/genomes/genbank'
