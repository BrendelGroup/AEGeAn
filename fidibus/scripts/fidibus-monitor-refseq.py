#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2016   Indiana University
# Copyright (c) 2016   The Regents of the University of California
#
# This file is part of fidibus (http://github.com/standage/fidibus) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import datetime
import glob
import os
import shutil
import subprocess
import sys
import fidibus


class GenomeDBCache(object):
    def __init__(self, db, cachedir='cache', tempdir='temp'):
        self.db = db
        self.cachedir = cachedir
        self.tempdir = tempdir
        self.quiet = False

    @property
    def cacheroot(self):
        return self.cachedir + '/' + self.db.label

    @property
    def newprefix(self):
        date = datetime.datetime.now().strftime('%Y-%m-%d')
        return '{c}/{d}-{a}'.format(c=self.cacheroot, d=date, a=self.db.acc)

    @property
    def oldprefix(self):
        return '{c}/*-{a}'.format(c=self.cacheroot, a=self.db.acc)

    def printlog(self, message):
        if self.quiet:
            return
        print(message, file=sys.stderr)

    def copy2cache(self, existingfile, newfile):
        newdir = os.path.dirname(newfile)
        subprocess.call(['mkdir', '-p', newdir])
        shutil.copy2(existingfile, newfile)

    def file_test(self, cachefile, testfile, newfile):
        """If test file is different from cache, copy to new file."""
        assert os.path.isfile(testfile)
        testsha1 = self.db.file_sha1(testfile)
        assert os.path.isfile(cachefile)
        cachesha1 = self.db.file_sha1(cachefile)
        if testsha1 == cachesha1:
            message = (
                'Testfile "{tf}" and cachefile "{cf}" match ({sha}); cache is '
                'up-to-date!'.format(tf=testfile, cf=cachefile, sha=testsha1)
            )
            self.printlog(message)
        else:
            message = (
                'Testfile "{tf}" (sha1={ts}) and cachefile "{cf}" (sha1={cs}) '
                'do not match match; cache is outdated; copying test file to '
                '"{nf}"!'.format(tf=testfile, ts=testsha1, cf=cachefile,
                                 cs=cachesha1, nf=newfile)
            )
            self.printlog(message)
            self.copy2cache(testfile, newfile)

    def run_tests(self):
        suffixes = ['genomic.fna.gz', 'genomic.gff.gz', 'protein.faa.gz']
        cachefilepatterns = [self.oldprefix + '_' + suf for suf in suffixes]
        newfiles = [self.newprefix + '_' + suf for suf in suffixes]
        tempfiles = [self.db.gdnapath, self.db.gff3path, self.db.protpath]
        for cfp, nf, tf in zip(cachefilepatterns, newfiles, tempfiles):
            cachefiles = sorted(glob.glob(cfp))
            if len(cachefiles) == 0:
                message = (
                    'No cachefile for comparison with {tf}; creating new cache'
                    ' and copying to "{nf}"'.format(tf=tf, nf=nf)
                )
                self.printlog(message)
                self.copy2cache(tf, nf)
            else:
                cachefile = cachefiles[-1]
                self.file_test(cachefile, tf, nf)


def get_parser():
    desc = 'Script to monitor RefSeq genomes and keep a local cache'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-v', '--version', action='version',
                        version='fidibus v%s' % fidibus.__version__)
    parser.add_argument('-w', '--work', default='DELETEME', metavar='DIR',
                        help='temporary working directory')
    parser.add_argument('-c', '--cache', default='cache', metavar='DIR',
                        help='cache directory; default is "cache/"')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='do not print debugging output')
    return parser


def main(args):
    registry = fidibus.registry.Registry()
    for label, config in registry.list_genomes():
        db = registry.genome(label, workdir=args.work)
        if 'source' not in db.config or db.config['source'] != 'refseq':
            continue
        if 'known_failing' in db.config:
            continue

        db.download()
        cache = GenomeDBCache(db, tempdir=args.work, cachedir=args.cache)
        cache.run_tests()

    shutil.rmtree(args.work)


if __name__ == '__main__':
    main(get_parser().parse_args())
