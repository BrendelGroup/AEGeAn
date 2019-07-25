#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016   Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2016   Indiana University
#
# This file is part of fidibus (http://github.com/standage/fidibus) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

"""Module for handling cd-hit output."""

from __future__ import print_function
import re


class ClusterSeq(object):
    """Object representing a sequence entry in a cd-hit cluster."""

    def __init__(self, line):
        self.rawdata = line.rstrip()
        values = re.compile(r'\s+').split(self.rawdata)
        self.index = int(values[0])
        self.length = int(values[1][:-3])
        self.defline = values[2]

    @property
    def accession(self):
        """
        Parse accession number from the specified format.

            >gnl|Tcas|XP_008191512.1
        """
        assert self.defline.startswith('>gnl|')
        acc = re.match(r'>gnl\|[^\|]+\|([^\|\n]+)', self.defline).group(1)
        assert acc.endswith('...')
        acc = acc[:-3]
        return acc

    @property
    def species(self):
        """
        Parse species label from the specified format.

            >gnl|Tcas|XP_008191512.1
        """
        return self.defline.split('|')[1]

    def __len__(self):
        return self.length


def parse_clusters(filehandle):
    """
    Iterate over clusters from a CD-HIT output file.

    Yields the cluster ID (a numeric string) and a list of sequence objects.
    """
    clusterid = None
    clusterseqs = list()
    for line in filehandle:
        if line.startswith('>'):
            if clusterid is not None:
                yield clusterid, clusterseqs
            clusterid = line.rstrip()[9:]  # Strip '>Cluster ' from front
            clusterseqs = list()
        else:
            seqinfo = ClusterSeq(line)
            clusterseqs.append(seqinfo)

    yield clusterid, clusterseqs


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def test_parse_clusters():
    """CD-HIT: parse clusters"""
    clusters = list()
    with open('testdata/misc/hymhub-head.clstr', 'r') as infile:
        for clusterid, clusterseqs in parse_clusters(infile):
            clusters.append(clusterseqs)

    assert len(clusters) == 3

    assert len(clusters[0]) == 1
    assert len(clusters[0][0]) == 25481
    assert clusters[0][0].species == 'Tcas'
    assert clusters[0][0].accession == 'XP_008191512.1'

    assert len(clusters[1]) == 1
    assert len(clusters[1][0]) == 22949
    assert clusters[1][0].species == 'Dmel'
    assert clusters[1][0].accession == 'NP_001260032.1'

    assert len(clusters[2]) == 3
    assert clusters[2][0].species == 'Amel'
    assert clusters[2][1].species == 'Bimp'
    assert clusters[2][2].species == 'Bter'
