#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument('species', help='4-letter species label')
parser.add_argument('gff3', type=argparse.FileType('r'), default=sys.stdin)
args = parser.parse_args()

for line in args.gff3:
    liilmatch = re.search(r'liil=(\d+)', line)
    riilmatch = re.search(r'riil=(\d+)', line)
    namematch = re.search(r'Name=([^;\n]+)', line)
    if not liilmatch or not riilmatch:
        continue

    lname = namematch.group(1)
    liil = liilmatch.group(1)
    riil = riilmatch.group(1)
    fields = '\t'.join([args.species, lname, liil, riil])
    print(fields)
