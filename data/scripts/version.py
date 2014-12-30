#!/usr/bin/env python
import re
import subprocess
import sys

with open("VERSION", "r") as vfile:
  semverstr = vfile.read().replace('\n', '')
  semver, stability = semverstr.split(' ')

logproc = subprocess.Popen(["git", "log"], stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
logout, logerr = logproc.communicate()
if logerr:
  sha1 = ""
  sha1slug = ""
  link = "https://github.com/standage/AEGeAn/releases/tag/"+ semver
  year = "2014"
else:
  sha1match = re.search("commit (\S+)", logout)
  assert sha1match, "could not find latest commit SHA1 hash"
  sha1 = sha1match.group(1)
  sha1slug = sha1[:10]
  link = "https://github.com/standage/AEGeAn/tree/"+ sha1
  
  yearmatch = re.search("Date:\s+.+(\d{4}) ", logout)
  assert yearmatch, "could not find year of latest commit"
  year = yearmatch.group(1)

print '#ifndef AEGEAN_VERSION_H'
print '#define AEGEAN_VERSION_H'
print '#define AGN_SEMANTIC_VERSION  "%s"' % semver
print '#define AGN_VERSION_STABILITY "%s"' % stability
print '#define AGN_VERSION_HASH      "%s"' % sha1
print '#define AGN_VERSION_HASH_SLUG "%s"' % sha1slug
print '#define AGN_VERSION_LINK      "%s"' % link
print '#define AGN_COPY_YEAR         "%s"' % year
print '#endif'
