#!/usr/bin/env python

# Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS
#
# The AEGeAn Toolkit is distributed under the ISC License. See
# the 'LICENSE' file in the AEGeAn source code distribution or
# online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

from __future__ import print_function
import argparse
import os
import re
import subprocess
import sys


def run_locuspocus(infile, outfile, delta, ilenfile=None, debug=False):
    command = 'locuspocus --verbose --namefmt=prelocus%lu'
    command += ' --delta %d' % delta
    command += ' --outfile %s.lp' % outfile
    command += ' --cds'
    if ilenfile:
        command += ' --ilens %s' % ilenfile
    command += ' %s' % infile
    if debug:
        print('command: %s' % command, file=sys.stderr)

    cmd = command.split(' ')
    subprocess.check_call(cmd)

    numloci = 0
    with open('%s.lp' % outfile, 'r') as fp:
        for line in fp:
            if '\tlocus\t' in line:
                numloci += 1
    if debug:
        print('numloci: %s' % numloci, file=sys.stderr)
    return numloci


def run_uloci(infile, outfile, counter, debug=False):
    command = 'uloci.py --counter %d' % counter
    command += ' --src AEGeAn::uloci.py'
    command += ' %s' % infile
    if debug:
        print('command: %s' % command, file=sys.stderr)

    cmd = command.split(' ')
    with open('%s.ul' % outfile, 'w') as fp:
        result = subprocess.call(cmd, stdout=fp)
    if result:
        print('command resulted in return status %d: %s' % (result, command),
              file=sys.stderr)
        exit(result)


def combine_output(outfile, namefmt):
    locusids = {}
    with open(outfile, 'w') as fp:
        command = 'gt gff3 -retainids -sort -tidy '
        command += '%s.lp %s.ul' % (outfile, outfile)
        cmd = command.split(' ')
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             universal_newlines=True)
        counter = 0
        while True:
            line = p.stdout.readline().rstrip()
            if line:
                if '\tlocus\t' in line:
                    counter += 1
                    locusname = namefmt % counter
                    line = re.sub('Name=[^;\n]+', 'Name=%s' % locusname, line)
                print(line, file=fp)
            else:
                break


if __name__ == '__main__':
    desc = 'Report all iLoci for an annotation file'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Debug mode')
    parser.add_argument('--ilenfile', type=str, default=None,
                        help='File to which iiLocus lengths will be written')
    parser.add_argument('--namefmt', type=str, default='locus%d',
                        help='An Name with a serial number is assigned to each'
                        ' locus; default format is "locus%%d"')
    parser.add_argument('--delta', type=int, default=500,
                        help='Delta for extending iLoci; default is 500')
    parser.add_argument('--out', type=str, default=None,
                        help='Output file; default is $infile.loci')
    parser.add_argument('infile', help='Input data in GFF3 format')
    args = parser.parse_args()
    if not args.out:
        args.out = '%s.loci' % args.infile

    numloci = run_locuspocus(args.infile, args.out, args.delta, args.ilenfile,
                             args.debug)
    run_uloci(args.infile, args.out, numloci + 1, args.debug)
    combine_output(args.out, args.namefmt)

    os.unlink('%s.lp' % args.out)
    os.unlink('%s.ul' % args.out)
