#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015-2016   Indiana University
# Copyright (c) 2016        The Regents of the University of California
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""Setup configuration for LocusPocus"""

import setuptools
import versioneer


setuptools.setup(name='LocusPocus',
                 version=versioneer.get_version(),
                 cmdclass=versioneer.get_cmdclass(),
                 description=('Explore eukaryotic genome composition and '
                              'organization with iLoci'),
                 url='http://github.com/BrendeGroup/AEGeAn',
                 author='Daniel Standage and collaborators',
                 author_email='daniel.standage@gmail.com',
                 license='BSD-3',
                 packages=['LocusPocus'],
                 scripts=['scripts/fidibus',
                          'scripts/fidibus-filens.py',
                          'scripts/fidibus-format-gff3.py',
                          'scripts/fidibus-glean-to-gff3.py',
                          'scripts/fidibus-namedup.py',
                          'scripts/fidibus-ilocus-summary.py',
                          'scripts/fidibus-pilocus-summary.py',
                          'scripts/fidibus-milocus-summary.py',
                          'scripts/fidibus-stats.py',
                          'scripts/fidibus-compact.py',
                          'scripts/fidibus-monitor-refseq.py',
                          'scripts/fidibus-uniq.py'],
                 install_requires=['pyyaml', 'pycurl'],
                 package_data={'LocusPocus': ['genome_configs/*.yml', 'genome_configs/*.txt']},
                 classifiers=[
                    'Development Status :: 4 - Beta',
                    'Environment :: Console',
                    'License :: OSI Approved :: BSD License',
                    'Programming Language :: Python :: 2.7',
                    'Programming Language :: Python :: 3.3',
                    'Programming Language :: Python :: 3.4',
                    'Programming Language :: Python :: 3.5',
                    'Topic :: Scientific/Engineering :: Bio-Informatics'
                 ],
                 zip_safe=False)
