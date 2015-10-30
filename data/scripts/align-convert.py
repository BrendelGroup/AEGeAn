#!/usr/bin/env python

# Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS
#
# The AEGeAn Toolkit is distributed under the ISC License. See
# the 'LICENSE' file in the AEGeAn source code distribution or
# online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

from __future__ import print_function
import re
import sys


class Peeker:
    """
    Minimal wrapper around an iterator that buffers lines for peeking ahead.
    Stolen shamelessly from http://stackoverflow.com/a/1517965/459780.
    """

    def __init__(self, iter):
        self.iter = iter
        self.buffer = []

    def __iter__(self):
        return self

    def next(self):
        if self.buffer:
            return self.buffer.pop(0)
        else:
            return next(self.iter)

    def __next__(self):
        return self.next()

    def peek(self, n=0):
        """Return an item n entries ahead in the iteration."""
        while n >= len(self.buffer):
            try:
                self.buffer.append(next(self.iter))
            except StopIteration:
                return None
        return self.buffer[n]


def align_convert(fp):
    """
    If an alignment feature is encountered, convert it from the 2-tiered match/
    match_part encoding to a 1-tiered match multifeature encoding.
    """
    moltypes = {'cDNA_match': 1, 'EST_match': 1, 'nucleotide_match': 1,
                'protein_match': 1}
    buffer = None
    for line in fp:
        line = line.rstrip()
        fields = line.split("\t")
        if len(fields) != 9:
            yield line
            continue
        ftype = fields[2]
        if ftype in moltypes:
            matchid = re.search('ID=([^;\n]+)', fields[8]).group(1)
            peekline = fp.peek()
            while not (peekline is None) and '\tmatch_part\t' in peekline:
                line = next(fp).rstrip()
                fields = line.split('\t')
                fields[8] = re.sub('ID=[^;\n]+;*', '', fields[8])
                fields[8] = re.sub('Parent=[^;\n]+;*', '', fields[8])
                fields[8] = 'ID=' + matchid + ';' + fields[8]
                fields[2] = ftype
                yield '\t'.join(fields)
                peekline = fp.peek()
        else:
            yield line

if __name__ == '__main__':
    fqiter = Peeker(sys.stdin)
    for entry in align_convert(fqiter):
        print(entry)
