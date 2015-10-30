#!/usr/bin/env python

# Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS
#
# The AEGeAn Toolkit is distributed under the ISC License. See
# the 'LICENSE' file in the AEGeAn source code distribution or
# online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

from __future__ import print_function
import re
import subprocess
import sys

with open("VERSION", "r") as vfile:
  semverstr = vfile.read().replace('\n', '')
  semver, stability = semverstr.split(' ')

try:
  logproc = subprocess.Popen(["git", "log"], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True)
  logout, logerr = logproc.communicate()
except:
  logerr = True

if logerr:
  sha1 = ""
  sha1slug = ""
  link = "https://github.com/standage/AEGeAn/releases/tag/"+ semver
  year = "2015"
else:
  sha1match = re.search("commit (\S+)", logout)
  assert sha1match, "could not find latest commit SHA1 hash"
  sha1 = sha1match.group(1)
  sha1slug = sha1[:10]
  link = "https://github.com/standage/AEGeAn/tree/"+ sha1

  yearmatch = re.search("Date:\s+.+(\d{4}) ", logout)
  assert yearmatch, "could not find year of latest commit"
  year = yearmatch.group(1)

print('#ifndef AEGEAN_VERSION_H')
print('#define AEGEAN_VERSION_H')
print('#define AGN_SEMANTIC_VERSION  "%s"' % semver)
print('#define AGN_VERSION_STABILITY "%s"' % stability)
print('#define AGN_VERSION_HASH      "%s"' % sha1)
print('#define AGN_VERSION_HASH_SLUG "%s"' % sha1slug)
print('#define AGN_VERSION_LINK      "%s"' % link)
print('#define AGN_COPY_YEAR         "%s"' % year)
print('#endif')
