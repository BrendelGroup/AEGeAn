#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015-2016   Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2015-2016   Indiana University
# Copyright (c) 2016        The Regents of the University of California
#
# This file is part of fidibus (http://github.com/standage/fidibus) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

"""Simple module for downloading data with PycURL"""

from __future__ import print_function
import gzip
import pycurl
import sys


def url_download(urldata, localpath, compress=False, follow=True):
    """
    Helper function for downloading remote data files with PycURL.

    - urldata: string(s), URL or list of URLs
    - localpath: path of the filename to which output will be written
    - compress: output compression
    """
    urls = urldata
    if isinstance(urldata, str):
        urls = [urldata]

    openfunc = open
    if compress is True:
        openfunc = gzip.open

    with openfunc(localpath, 'wb') as out:
        for url in urls:
            try:
                c = pycurl.Curl()
                c.setopt(c.URL, url)
                c.setopt(c.WRITEDATA, out)
                if follow:
                    c.setopt(c.FOLLOWLOCATION, True)
                c.perform()
                c.close()
            except pycurl.error as e:
                print('Error: unable to download URL::', url, file=sys.stderr)
                raise e
