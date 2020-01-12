#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import sys
lines = set()
for line in sys.stdin:
    if line in lines:
        continue
    lines.add(line)
    print(line, end='')
