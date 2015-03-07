#!/usr/bin/env python

# Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS
#
# The AEGeAn Toolkit is distributed under the ISC License. See
# the 'LICENSE' file in the AEGeAn source code distribution or
# online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

import argparse
import os
import re
import subprocess
import sys

def parse_fasta(fp):
  """
  Stolen shamelessly from http://stackoverflow.com/a/7655072/459780.
  """
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield (name, ''.join(seq))
      name, seq = line, []
    else:
      seq.append(line)
  if name: yield (name, ''.join(seq))

def seq_len(fp):
  for defline, seq in parse_fasta(fp):
    seqid = defline[1:].split(" ")[0]
    yield seqid, len(seq)

def parse_gff3(fp):
  for line in fp:
    if line.startswith("##gff-version") or line.startswith("##sequence-region"):
      continue
    yield line.rstrip()

if __name__ == "__main__":
  desc = "Correct a GFF3's sequence-region declarations with Fasta sequences"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument("gff3", help="Annotation in GFF3 format; use - for stdin")
  parser.add_argument("fasta", help="Corresponding sequences in Fasta format")
  args = parser.parse_args()

  print "##gff-version   3"
  with open(args.fasta, "r") as fp:
    for seqid, seqlen in seq_len(fp):
      print "##sequence-region   %s 1 %d" % (seqid, seqlen)

  gff3in = sys.stdin
  if args.gff3 != "-":
    gff3in = open(args.gff3, "r")
  for entry in parse_gff3(gff3in):
    print entry
