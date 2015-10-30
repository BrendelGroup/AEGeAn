#!/usr/bin/env python

# Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS
#
# The AEGeAn Toolkit is distributed under the ISC License. See
# the 'LICENSE' file in the AEGeAn source code distribution or
# online at https://github.com/standage/AEGeAn/blob/master/LICENSE.
#
# Usage: check-gff3.py < annot.gff3

from __future__ import print_function
import re
import sys

types = {}
counts = {}
topleveltypes = {}
relationships = {}

for line in sys.stdin:
    line = line.rstrip()
    if line == '' or line.startswith('#'):
        continue

    values = line.split('\t')
    assert len(values) == 9

    if not values[2] in counts:
        counts[values[2]] = 0
    counts[values[2]] += 1

    m = re.search('ID=([^;\n]+)', values[8])
    if m:
        featureid = m.group(1)
        types[featureid] = values[2]

    m = re.search('Parent=([^;\n]+)', values[8])
    if m:
        parent = m.group(1)
        parenttype = types[parent]
        assert parenttype
        if not values[2] in relationships:
            relationships[values[2]] = {}
        if parenttype not in relationships[values[2]]:
            relationships[values[2]][parenttype] = 0
        relationships[values[2]][parenttype] += 1
    else:
        if values[2] not in topleveltypes:
            topleveltypes[values[2]] = 0
        topleveltypes[values[2]] += 1

print('\n===== Types =====')
featuretypes = counts.keys()
featuretypes.sort()
for featuretype in featuretypes:
    count = counts[featuretype]
    topcount = 0
    if featuretype in topleveltypes:
        topcount = topleveltypes[featuretype]
    print('%s (%lu total, %lu top-level)' % (featuretype, count, topcount))

print('\n===== Relationships (child --> parent) =====')
childtypes = relationships.keys()
childtypes.sort()
for childtype in childtypes:
    parenttypes = list(relationships[childtype].keys())
    parenttypes.sort()
    for parenttype in parenttypes:
        count = relationships[childtype][parenttype]
        print('%s --> %s (%lu occurrences)' % (childtype, parenttype, count))
