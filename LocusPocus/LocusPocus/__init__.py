#!/usr/bin/env python
# coding: utf-8
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015-2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""Package-wide configuration"""

# Package modules
from __future__ import print_function
from . import registry
from . import download
from . import fasta
from . import cdhit
from . import genomedb
from . import refseq
from . import crg
from . import hymbase
from . import tair
from . import generic
from . import iloci
from . import proteins
from . import mrnas
from . import exons
from . import stats
try:
    FileNotFoundError
except NameError:  # pragma: no cover
    FileNotFoundError = IOError

# Custom modules
from . import am10
from . import pdom

# Versioneer
from . import _version
__version__ = _version.get_versions()['version']


# Unit test fixtures (can't figure out how to do package-scope global fixtures
# with nose's setup and teardown mechanism).
test_registry = registry.Registry()
test_registry_supp = registry.Registry()
try:
    # This will only work when the current working directory is the LocusPocus
    # root directory. Fine since it's only for development.
    test_registry_supp.update('testdata/conf')
except FileNotFoundError:  # pragma: no cover
    pass


sources = {
    'refseq': 'NCBI RefSeq',
    'genbank': 'NCBI Genbank',
    'beebase': 'BeeBase Consortium',
    'crg': 'Wasp/ant genome project (Centro de Regulación Genómica)',
    'hymbase': 'Hymenoptera Genome Database',
    'pdom': 'Paper wasp genome project (Toth Lab)',
    'tair': 'TAIR6 (The Arabidopsis Information Resource)',
    'am10': 'Amel OGSv1.0 (Honeybee Genome Sequencing Consortium)',
    'local': 'user-supplied genome (local file system)'
}

dbtype = {
    'refseq': refseq.RefSeqDB,
    'genbank': refseq.GenbankDB,
    'beebase': hymbase.BeeBaseDB,
    'crg': crg.CrgDB,
    'hymbase': hymbase.HymBaseDB,
    'pdom': pdom.PdomDB,
    'tair': tair.TairDB,
    'am10': am10.Am10DB,
    'local': generic.GenericDB,
}
