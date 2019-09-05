#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""Simple module for reading, writing, subsetting, and comparing sequences."""

from __future__ import print_function
import sys
try:
    from StringIO import StringIO
except ImportError:  # pragma: no cover
    from io import StringIO


def parse(data):
    """
    Load sequences in Fasta format.

    This generator function yields a tuple containing a defline and a sequence
    for each record in the Fasta data. Stolen shamelessly from
    http://stackoverflow.com/a/7655072/459780.
    """
    name, seq = None, []
    for line in data:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


def format(seq, linewidth=70, outstream=sys.stdout):
    """Print a sequence in a readable format."""
    if linewidth == 0 or len(seq) <= linewidth:
        print(seq, file=outstream)
        return

    i = 0
    while i < len(seq):
        print(seq[i:i+linewidth], file=outstream)
        i += linewidth


def select(idstream, seqstream):
    ids = dict()
    for line in idstream:
        seqid = line.rstrip()
        ids[seqid] = True

    for defline, seq in parse(seqstream):
        seqid = defline[1:].split()[0]
        if seqid in ids:
            yield defline, seq


def compare(stream1, stream2):
    seqs1 = dict()
    seqs2 = dict()
    for defline, seq in parse(stream1):
        seqs1[defline] = seq
    for defline, seq in parse(stream2):
        seqs2[defline] = seq
    return seqs1 == seqs2


def test_parse():
    """Fasta: parsing"""
    data = ('>seq1\n'
            'ACGT\n'
            '>seq2\n'
            'ACGTACGTACGTACGT\n'
            'ACGTACGTACGTACGT\n')
    seqs = {'>seq1': 'ACGT',
            '>seq2': 'ACGTACGTACGTACGTACGTACGTACGTACGT'}
    testseqs = dict()
    for defline, seq in parse(data.split('\n')):
        testseqs[defline] = seq
    assert seqs == testseqs, \
        'parsed sequence mismatch: %r %r' % (seqs, testseqs)
    seqs = []
    testseqs = list(parse(['']))
    assert seqs == testseqs, 'empty seqfile fail %r %r' % (seqs, testseqs)


def test_format_seq():
    """Fasta: sequence formatting"""
    seq = ('TCTCCCTCCA'
           'ACGCCCGAAC'
           'GTGTCTGCTC'
           'ATTTCAAGCA'
           'CACGCATGAA'
           'CGGCATCGCG'
           'CAGACGTGCG'
           'AGCGAGCGCA')

    sio = StringIO()
    format(seq, linewidth=0, outstream=sio)
    assert sio.getvalue() == seq + '\n'
    sio.close()

    sio = StringIO()
    format(seq, linewidth=40, outstream=sio)
    assert sio.getvalue() == ('TCTCCCTCCA'
                              'ACGCCCGAAC'
                              'GTGTCTGCTC'
                              'ATTTCAAGCA\n'
                              'CACGCATGAA'
                              'CGGCATCGCG'
                              'CAGACGTGCG'
                              'AGCGAGCGCA\n')
    sio.close()

    sio = StringIO()
    format(seq, linewidth=20, outstream=sio)
    assert sio.getvalue() == ('TCTCCCTCCA'
                              'ACGCCCGAAC\n'
                              'GTGTCTGCTC'
                              'ATTTCAAGCA\n'
                              'CACGCATGAA'
                              'CGGCATCGCG\n'
                              'CAGACGTGCG'
                              'AGCGAGCGCA\n')
    sio.close()


def test_select():
    """Fasta: sequence extraction"""
    data = ('>seq1\n'
            'ACGT\n'
            '>seq2\n'
            'ACGTACGTACGTACGT\n'
            'ACGTACGTACGTACGT')
    seqs = {'>seq1': 'ACGT'}
    testseqs = dict()
    for defline, seq in select(['seq1'], data.split('\n')):
        testseqs[defline] = seq
    assert seqs == testseqs, \
        'extracted sequence mismatch: %r %r' % (seqs, testseqs)


def test_compare():
    """Fasta: order-independent sequence comparison"""
    data1 = ('>seq1\n'
             'ACGT\n'
             '>seq2\n'
             'ACGTACGTACGTACGT\n'
             'ACGTACGTACGTACGT')
    data2 = ('>seq2\n'
             'ACGTACGTACGTACGT\n'
             'ACGTACGTACGTACGT\n'
             '>seq1\n'
             'ACGT')

    assert compare(data1.split('\n'), data2.split('\n')), \
        'sequence comparison failed'
