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
import sys
lines = set()
for line in sys.stdin:
    if line in lines:
        continue
    lines.add(line)
    print(line, end='')
