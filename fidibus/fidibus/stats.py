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
import subprocess
import sys


def compute(db, logstream=sys.stderr):  # pragma: no cover
    if logstream is not None:
        logmsg = '[Genome: %s] ' % db.config['species']
        logmsg += 'calculating feature statistics'
        print(logmsg, file=logstream)

    prefix = '%s/%s/%s' % (db.workdir, db.label, db.label)
    prefix3 = (prefix, prefix, prefix)

    command = 'fidibus-stats.py --species ' + db.label
    command += ' --iloci %s.iloci.gff3 %s.iloci.fa %s.iloci.tsv' % prefix3
    command += ' --miloci %s.miloci.gff3 %s.miloci.fa %s.miloci.tsv' % prefix3
    command += (' --prnas %s.ilocus.mrnas.gff3 %s.pre-mrnas.fa '
                '%s.pre-mrnas.tsv' % prefix3)
    command += ' --mrnas %s.mrnas.gff3 %s.mrnas.fa %s.mrnas.tsv' % prefix3
    command += ' --cds %s.ilocus.mrnas.gff3 %s.cds.fa %s.cds.tsv' % prefix3
    command += (' --exons %s.ilocus.mrnas.gff3 %s.exons.fa '
                '%s.exons.tsv' % prefix3)
    command += (' --introns %s.ilocus.mrnas.gff3 %s.introns.fa '
                '%s.introns.tsv' % prefix3)

    cmd = command.split(' ')
    subprocess.check_call(cmd)
