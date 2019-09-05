#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Indiana University
# Copyright (c) 2016   The Regents of the University of California
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import filecmp
import sys
import fidibus


def ids(db, logstream=sys.stderr):  # pragma: no cover
    """
    Retrieve protein IDs/accessions from the genome annotation.

    The `db` variable, a `GenomeDB` object, must implement a `gff3_protids`
    method for this retrieval.
    """
    if logstream is not None:
        logmsg = '[Genome: %s] retrieving protein IDs' % db.config['species']
        print(logmsg, file=logstream)

    specdir = '%s/%s' % (db.workdir, db.label)
    infile = '%s/%s.ilocus.mrnas.gff3' % (specdir, db.label)
    outfile = '%s/%s.protids.txt' % (specdir, db.label)
    with open(infile, 'r') as instream, open(outfile, 'w') as outstream:
        for protid in db.gff3_protids(instream):
            print(protid, file=outstream)


def sequences(db, logstream=sys.stderr):
    """Extract protein sequences."""
    if logstream is not None:  # pragma: no cover
        logmsg = '[Genome: %s] ' % db.config['species']
        logmsg += 'extracting protein sequences'
        print(logmsg, file=logstream)

    specdir = '%s/%s' % (db.workdir, db.label)
    idfile = '%s/%s.protids.txt' % (specdir, db.label)
    seqfile = '%s/%s.all.prot.fa' % (specdir, db.label)
    outfile = '%s/%s.prot.fa' % (specdir, db.label)
    with open(idfile, 'r') as idstream, \
            open(seqfile, 'r') as seqstream, \
            open(outfile, 'w') as outstream:
        for defline, seq in fidibus.fasta.select(idstream, seqstream):
            defline = '>gnl|%s|%s' % (db.label, defline[1:])
            print(defline, file=outstream)
            fidibus.fasta.format(seq, outstream=outstream)


def mapping(db, only_reps=False, logstream=sys.stderr):
    """
    Retrieve mapping of protein IDs to iLocus IDs.

    The `db` variable, a `GenomeDB` object, must implement a `protein_mapping`
    method for this retrieval.
    """

    specdir = '%s/%s' % (db.workdir, db.label)
    infile = '%s/%s.iloci.gff3' % (specdir, db.label)
    if only_reps:
        if logstream is not None:  # pragma: no cover
            logmsg = '[Genome: %s] ' % db.config['species']
            logmsg += 'parsing protein->iLocus reps mapping'
            print(logmsg, file=logstream)
        outfile = '%s/%s.protein2ilocus.repr.tsv' % (specdir, db.label)
        repfile = '%s/%s.protids.txt' % (specdir, db.label)
        protreps = dict()
        with open(repfile, 'r') as repstream:
            for line in repstream:
                protid = line.strip()
                protreps[protid] = True
    else:
        if logstream is not None:  # pragma: no cover
            logmsg = '[Genome: %s] ' % db.config['species']
            logmsg += 'parsing protein->iLocus alls mapping'
            print(logmsg, file=logstream)
        outfile = '%s/%s.protein2ilocus.tsv' % (specdir, db.label)
    with open(infile, 'r') as instream, open(outfile, 'w') as outstream:
        print('ProteinID', 'piLocusID', sep='\t', file=outstream)
        for protid, ilocusid in db.protein_mapping(instream):
            if not only_reps or protid in protreps:
                print(protid, ilocusid, sep='\t', file=outstream)


# -----------------------------------------------------------------------------
# Driver function
# -----------------------------------------------------------------------------

def prepare(db, logstream=sys.stderr):  # pragma: no cover
    ids(db, logstream=logstream)
    sequences(db, logstream=logstream)
    mapping(db, only_reps=False, logstream=logstream)
    mapping(db, only_reps=True, logstream=logstream)


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------


def test_protein_sequence():
    """Breakdown: select protein sequences"""
    db = fidibus.test_registry.genome('Scer', workdir='testdata/demo-workdir')
    sequences(db, logstream=None)
    outfile = 'testdata/demo-workdir/Scer/Scer.prot.fa'
    testfile = 'testdata/fasta/scer-few-prots.fa'
    assert filecmp.cmp(outfile, testfile), 'Protein sequence selection failed'
