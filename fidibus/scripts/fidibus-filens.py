#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015   Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2015   Indiana University
#
# This file is part of fidibus (http://github.com/standage/fidibus) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
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
