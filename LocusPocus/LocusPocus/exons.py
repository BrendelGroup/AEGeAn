#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import filecmp
import re
import subprocess
import sys
import LocusPocus


def cds_sequences(db, logstream=sys.stderr):
    if logstream is not None:  # pragma: no cover
        logmsg = '[Genome: %s] ' % db.config['species']
        logmsg += 'extracting coding sequences'
        print(logmsg, file=logstream)
    specdir = '%s/%s' % (db.workdir, db.label)

    gff3infile = '%s/%s.gff3' % (specdir, db.label)
    fastainfile = '%s/%s.gdna.fa' % (specdir, db.label)
    outfile = '%s/%s.all.cds.fa' % (specdir, db.label)
    command = 'xtractore --type=CDS --outfile=%s ' % outfile
    command += '%s %s' % (gff3infile, fastainfile)
    cmd = command.split(' ')
    subprocess.check_call(cmd)

    gff3infile = '%s/%s.ilocus.mrnas.gff3' % (specdir, db.label)
    fastainfile = '%s/%s.gdna.fa' % (specdir, db.label)
    outfile = '%s/%s.cds.fa' % (specdir, db.label)
    command = 'xtractore --type=CDS --outfile=%s ' % outfile
    command += '%s %s' % (gff3infile, fastainfile)
    cmd = command.split(' ')
    subprocess.check_call(cmd)


def exon_sequences(db, logstream=sys.stderr):
    if logstream is not None:  # pragma: no cover
        logmsg = '[Genome: %s] ' % db.config['species']
        logmsg += 'extracting exon sequences'
        print(logmsg, file=logstream)
    specdir = '%s/%s' % (db.workdir, db.label)

    gff3infile = '%s/%s.ilocus.mrnas.gff3' % (specdir, db.label)
    fastainfile = '%s/%s.gdna.fa' % (specdir, db.label)
    outfile = '%s/%s.exons.fa' % (specdir, db.label)
    command = 'xtractore --type=exon --outfile=%s ' % outfile
    command += '%s %s' % (gff3infile, fastainfile)
    cmd = command.split(' ')
    subprocess.check_call(cmd)


def parse_intron_accessions(instream):
    moltypes = [
        'mRNA', 'tRNA', 'ncRNA', 'transcript', 'primary_transcript',
        'guide_RNA', 'V_gene_segment', 'D_gene_segment', 'J_gene_segment',
        'C_gene_segment'
    ]
    id_to_accession = dict()
    for line in instream:
        line = line.rstrip()
        idmatch = re.search('ID=([^;\n]+)', line)
        accmatch = re.search('accession=([^;\n]+)', line)
        if idmatch and accmatch:
            molid = idmatch.group(1)
            accession = accmatch.group(1)
            id_to_accession[molid] = accession

        if '\tintron\t' in line:
            parentid = re.search('Parent=([^;\n]+)', line).group(1)
            assert ',' not in parentid, parentid
            accession = id_to_accession[parentid]
            line += ';accession=%s' % accession

        yield line


def intron_sequences(db, logstream=sys.stderr):
    if logstream is not None:  # pragma: no cover
        logmsg = '[Genome: %s] ' % db.config['species']
        logmsg += 'extracting intron sequences'
        print(logmsg, file=logstream)
    specdir = '%s/%s' % (db.workdir, db.label)

    infile = '%s/%s.ilocus.mrnas.gff3' % (specdir, db.label)
    outfile = '%s/%s.with-introns.gff3' % (specdir, db.label)
    command = 'canon-gff3 --outfile=%s %s' % (outfile, infile)
    cmd = command.split(' ')
    subprocess.check_call(cmd)

    infile = '%s/%s.ilocus.mrnas.gff3' % (specdir, db.label)
    outfile = '%s/%s.with-introns.gff3' % (specdir, db.label)
    with open(infile, 'r') as instream, open(outfile, 'w') as outstream:
        for line in parse_intron_accessions(instream):
            print(line, file=outstream)

    gff3infile = '%s/%s.with-introns.gff3' % (specdir, db.label)
    fastainfile = '%s/%s.gdna.fa' % (specdir, db.label)
    outfile = '%s/%s.introns.fa' % (specdir, db.label)
    command = 'xtractore --type=intron --outfile=%s ' % outfile
    command += '%s %s' % (gff3infile, fastainfile)
    cmd = command.split(' ')
    subprocess.check_call(cmd)


# -----------------------------------------------------------------------------
# Driver function
# -----------------------------------------------------------------------------


def prepare(db, logstream=sys.stderr):  # pragma: no cover
    exon_sequences(db, logstream=logstream)
    intron_sequences(db, logstream=logstream)
    cds_sequences(db, logstream=logstream)


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------


def test_coding_sequences():
    """Breakdown: extract coding sequences"""
    db = LocusPocus.test_registry.genome('Atha',
                                         workdir='testdata/demo-workdir')
    cds_sequences(db, logstream=None)
    outfile = 'testdata/demo-workdir/Atha/Atha.cds.fa'
    testfile = 'testdata/fasta/atha-cds.fa'
    assert filecmp.cmp(outfile, testfile), 'coding sequence extraction failed'


def test_exon_sequences():
    """Breakdown: extract exon sequences"""
    db = LocusPocus.test_registry.genome('Atha',
                                         workdir='testdata/demo-workdir')
    exon_sequences(db, logstream=None)
    outfile = 'testdata/demo-workdir/Atha/Atha.exons.fa'
    testfile = 'testdata/fasta/atha-exons.fa'
    assert filecmp.cmp(outfile, testfile), 'exon sequence extraction failed'


def test_intron_sequences():
    """Breakdown: extract intron sequences"""
    db = LocusPocus.test_registry.genome('Atha',
                                         workdir='testdata/demo-workdir')
    intron_sequences(db, logstream=None)
    outfile = 'testdata/demo-workdir/Atha/Atha.introns.fa'
    testfile = 'testdata/fasta/atha-introns.fa'
    assert filecmp.cmp(outfile, testfile), 'intron sequence extraction failed'
