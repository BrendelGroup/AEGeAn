#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------
"""
Process a GFF3 file, and for any feature that has an ID but lacks a `Name`
attribute, copy the ID attribute to the Name attribute.
"""

from __future__ import print_function
import re
import sys

for line in sys.stdin:
    line = line.rstrip()
    idmatch = re.search("ID=([^;\n]+)", line)
    namematch = re.search("Name=([^;\n]+)", line)
    if idmatch and not namematch:
        line += ";Name=%s" % idmatch.group(1)
    print(line)
