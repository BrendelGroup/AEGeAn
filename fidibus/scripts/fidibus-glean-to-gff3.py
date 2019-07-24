#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2016   Indiana University
#
# This file is part of fidibus (http://github.com/standage/fidibus) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from __future__ import print_function
import re
import sys

for line in sys.stdin:
    fields = line.split('\t')
    assert len(fields) == 9
    ftype = fields[2]
    assert ftype in ['mRNA', 'CDS']

    attrstring = fields[8]
    protmatch = re.search(r'GenePrediction (\S+)', attrstring)
    assert protmatch, line
    protid = protmatch.group(1)
    featid = protid
    if ftype == 'CDS':
        featid += '-CDS'
    attrs = list()
    attrs.append('ID=%s' % featid)
    if ftype == 'CDS':
        attrs.append('Parent=%s' % protid)
    attrs.append('Name=%s' % protid)
    fields[8] = ';'.join(attrs)

    newentry = '\t'.join(fields)
    print(newentry)
