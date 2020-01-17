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


def mrna_exons(instream, convert=False, keepMrnas=False, usecds=False):
    """
    Parse exons from a data stream in GFF3 format.

    - If `convert` is true, convert the feature type of the exons to `mRNA`,
      creating a mature mRNA multi-feature.
    - If `keepMrnas` is true, return mRNA features in addition to exon
      features.
    - If `usecds` is true, parse exon structure from `CDS` features instead of
      `exon` features.
    """
    mrnaids = {}
    for line in instream:
        line = line.rstrip()
        fields = line.split('\t')
        if len(fields) != 9:
            continue
        exontype = 'exon'
        if usecds:
            exontype = 'CDS'

        if fields[2] == 'mRNA':
            mrnaid = re.search('ID=([^;\n]+)', fields[8]).group(1)
            accmatch = re.search('accession=([^;\n]+)', fields[8])
            assert accmatch, 'Unable to parse mRNA accession: %s' % fields[8]
            mrnaacc = accmatch.group(1)
            mrnaids[mrnaid] = 1
            if not convert and keepMrnas:  # pragma: no cover
                fields[8] = re.sub('Parent=[^;\n]+;*', '', fields[8])
                yield '\t'.join(fields)

        elif fields[2] == exontype:
            parentid = re.search('Parent=([^;\n]+)', fields[8]).group(1)
            fields[7] = '.'
            if parentid in mrnaids:
                if convert:
                    fields[2] = 'mRNA'
                    fields[8] = re.sub('ID=[^;\n]+;*', '', fields[8])
                    fields[8] = fields[8].replace('Parent=', 'ID=')
                    if 'accession=' not in fields[8]:  # pragma: no cover
                        fields[8] += ';accession=' + mrnaacc
                else:
                    if not keepMrnas:  # pragma: no cover
                        fields[8] = re.sub('Parent=[^;\n]+;*', '', fields[8])
                yield '\t'.join(fields)


def mature_mrna_intervals(db, logstream=sys.stderr):
    """
    Parse gene model structures and create mRNA mutli-features.

    Extracting the sequence of a pre-mRNA is trivial, but extracting the
    sequence of a mature mRNA (sans introns) requires some additional work.
    This function creates a new GFF3 file containing mRNA multi-features,
    enabling sequence extraction via xtractore.
    """
    if logstream is not None:  # pragma: no cover
        logmsg = '[Genome: %s] ' % db.config['species']
        logmsg += 'calculating mature mRNA intervals'
        print(logmsg, file=logstream)
    specdir = '%s/%s' % (db.workdir, db.label)

    infile = '%s/%s.gff3' % (specdir, db.label)
    outfile = '%s/%s.mrnas.temp' % (specdir, db.label)
    usecds = False
    if repr(db) in ['BeeBase', 'OGS1.0']:
        usecds = True
    with open(infile, 'r') as instream, open(outfile, 'w') as outstream:
        for exon in mrna_exons(instream, convert=True, usecds=usecds):
            print(exon, file=outstream)

    infile = '%s/%s.ilocus.mrnas.gff3' % (specdir, db.label)
    outfile = '%s/%s.ilocus.mrnas.temp' % (specdir, db.label)
    with open(infile, 'r') as instream, open(outfile, 'w') as outstream:
        for exon in mrna_exons(instream, convert=True, usecds=usecds):
            print(exon, file=outstream)

    inpatterns = ['%s/%s.mrnas.temp', '%s/%s.ilocus.mrnas.temp']
    outpatterns = ['%s/%s.all.mrnas.gff3', '%s/%s.mrnas.gff3']
    for inpattern, outpattern in zip(inpatterns, outpatterns):
        inf = inpattern % (specdir, db.label)
        otf = outpattern % (specdir, db.label)
        command = 'gt gff3 -retainids -sort -tidy -force -o %s %s' % (otf, inf)
        cmd = command.split(' ')
        proc = subprocess.Popen(cmd, stderr=subprocess.PIPE,
                                universal_newlines=True)
        _, stderr = proc.communicate()
        for line in stderr.split('\n'):  # pragma: no cover
            if 'has not been previously introduced' not in line and \
               'does not begin with "##gff-version"' not in line and \
               line != '':
                print(line, file=logstream)


def sequences(db, logstream=sys.stderr):
    if logstream is not None:  # pragma: no cover
        logmsg = '[Genome: %s] ' % db.config['species']
        logmsg += 'extracting pre-mRNA and mature mRNA sequences'
        print(logmsg, file=logstream)
    specdir = '%s/%s' % (db.workdir, db.label)

    # All pre-mRNA sequences
    gff3infile = '%s/%s.gff3' % (specdir, db.label)
    fastainfile = '%s/%s.gdna.fa' % (specdir, db.label)
    outfile = '%s/%s.all.pre-mrnas.fa' % (specdir, db.label)
    command = 'xtractore --type=mRNA --outfile=%s ' % outfile
    command += '%s %s' % (gff3infile, fastainfile)
    cmd = command.split(' ')
    subprocess.check_call(cmd)

    # All mature mRNA sequences
    gff3infile = '%s/%s.all.mrnas.gff3' % (specdir, db.label)
    fastainfile = '%s/%s.gdna.fa' % (specdir, db.label)
    outfile = '%s/%s.all.mrnas.fa' % (specdir, db.label)
    command = 'xtractore --type=mRNA --outfile=%s ' % outfile
    command += '%s %s' % (gff3infile, fastainfile)
    cmd = command.split(' ')
    subprocess.check_call(cmd)

    # Representative pre-mRNA sequences
    idfile = '%s/%s.mrnas.txt' % (specdir, db.label)
    seqfile = '%s/%s.all.pre-mrnas.fa' % (specdir, db.label)
    outfile = '%s/%s.pre-mrnas.fa' % (specdir, db.label)
    with open(idfile, 'r') as idstream, \
            open(seqfile, 'r') as seqstream, \
            open(outfile, 'w') as outstream:
        for defline, seq in LocusPocus.fasta.select(idstream, seqstream):
            print(defline, file=outstream)
            LocusPocus.fasta.format(seq, outstream=outstream)

    # Representative mature mRNA sequences
    idfile = '%s/%s.mrnas.txt' % (specdir, db.label)
    seqfile = '%s/%s.all.mrnas.fa' % (specdir, db.label)
    outfile = '%s/%s.mrnas.fa' % (specdir, db.label)
    with open(idfile, 'r') as idstream, \
            open(seqfile, 'r') as seqstream, \
            open(outfile, 'w') as outstream:
        for defline, seq in LocusPocus.fasta.select(idstream, seqstream):
            print(defline, file=outstream)
            LocusPocus.fasta.format(seq, outstream=outstream)


# -----------------------------------------------------------------------------
# Driver function
# -----------------------------------------------------------------------------

def prepare(db, logstream=sys.stderr):  # pragma: no cover
    mature_mrna_intervals(db, logstream=logstream)
    sequences(db, logstream=logstream)


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------


def test_mature_mrna_intervals():
    """Breakdown: define mature mRNA intervals"""
    db = LocusPocus.test_registry.genome('Atha', workdir='testdata/demo-workdir')
    mature_mrna_intervals(db, logstream=None)

    outfile = 'testdata/demo-workdir/Atha/Atha.all.mrnas.gff3'
    testfile = 'testdata/gff3/atha-all-mrnas.gff3'
    assert filecmp.cmp(outfile, testfile), 'mature mRNA interval ID failed'

    outfile = 'testdata/demo-workdir/Atha/Atha.mrnas.gff3'
    testfile = 'testdata/gff3/atha-mrnas.gff3'
    assert filecmp.cmp(outfile, testfile), 'mature mRNA interval ID failed'

    db = LocusPocus.test_registry.genome('Dnov', workdir='testdata/demo-workdir')
    mature_mrna_intervals(db, logstream=None)

    outfile = 'testdata/demo-workdir/Dnov/Dnov.all.mrnas.gff3'
    testfile = 'testdata/gff3/dnov-all-mrnas.gff3'
    assert filecmp.cmp(outfile, testfile), 'mature mRNA interval ID failed'

    outfile = 'testdata/demo-workdir/Dnov/Dnov.mrnas.gff3'
    testfile = 'testdata/gff3/dnov-mrnas.gff3'
    assert filecmp.cmp(outfile, testfile), 'mature mRNA interval ID failed'


def test_mrna_sequences():
    """Breakdown: extract mRNA sequences"""
    db = LocusPocus.test_registry.genome('Atha', workdir='testdata/demo-workdir')
    sequences(db, logstream=None)

    outfile = 'testdata/demo-workdir/Atha/Atha.all.pre-mrnas.fa'
    testfile = 'testdata/fasta/atha-all-pre-mrnas.fa'
    with open(outfile, 'r') as out, open(testfile, 'r') as test:
        assert LocusPocus.fasta.compare(out, test) is True, \
            'all pre-mRNA seq extraction failed'

    outfile = 'testdata/demo-workdir/Atha/Atha.pre-mrnas.fa'
    testfile = 'testdata/fasta/atha-pre-mrnas.fa'
    with open(outfile, 'r') as out, open(testfile, 'r') as test:
        assert LocusPocus.fasta.compare(out, test) is True, \
            'pre-mRNA seq extraction failed'

    outfile = 'testdata/demo-workdir/Atha/Atha.all.pre-mrnas.fa'
    testfile = 'testdata/fasta/atha-all-pre-mrnas.fa'
    with open(outfile, 'r') as out, open(testfile, 'r') as test:
        assert LocusPocus.fasta.compare(out, test) is True, \
            'all mRNA seq extraction failed'

    outfile = 'testdata/demo-workdir/Atha/Atha.mrnas.fa'
    testfile = 'testdata/fasta/atha-mrnas.fa'
    with open(outfile, 'r') as out, open(testfile, 'r') as test:
        assert LocusPocus.fasta.compare(out, test) is True, \
            'mature mRNA seq extraction failed'
