#!/usr/bin/env python

# Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS
#
# The AEGeAn Toolkit is distributed under the ISC License. See
# the 'LICENSE' file in the AEGeAn source code distribution or
# online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

import argparse
import re
import sys

def parse_gff3(fp):
  """
  Process entire (sorted) GFF3 file. Store lengths of sequences as reported by
  `##sequence-region` pragmas, and then discard lengths of sequences with
  annotated features. Sequence IDs and lengths yielded by the function
  correspond to unannotated sequences.
  """
  seqlens = {}
  for line in fp:
    if line.startswith("##sequence-region"):
      pattern = "##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)"
      seqmatch = re.search(pattern, line)
      assert seqmatch, "unable to parse sequence region pragma: %s" % line
      seqid = seqmatch.group(1)
      start = int(seqmatch.group(2))
      end   = int(seqmatch.group(3))
      assert start == 1, "assumption: start == 1: %s" % line
      seqlens[seqid] = end - start + 1
    elif len(line.split("\t")) == 9:
      seqid = line.split("\t")[0]
      seqlens.pop(seqid, None)

  for seqid in sorted(seqlens.keys()):
    yield seqid, seqlens[seqid]

def run_parse_gff3(gff3string, src=".", idformat="locus%d", counter=1):
  """
  Driver function for small examples (i.e. testing).
  """
  output = ""
  for seqid, seqlen in parse_gff3(gff3string.splitlines()):
    locusid = idformat % counter
    counter += 1
    attrs = "ID=%s;fragment=true;unannot=true;" % locusid
    attrs += "iLocus_type=iiLocus;effective_length=%d\n" % seqlen
    fields = [seqid, src, "locus", "1", str(seqlen),
              ".", ".", ".", attrs]
    output += "\t".join(fields)
  return output

def test_parse_gff3():
  """
  Unit tests for `parse_gff3` function. Throws exception in case of error, so no
  output is good output!
  """
  gff3 = ("##gff-version   3\n"
          "##sequence-region   chr1 1 1000000\n"
          "##sequence-region   chr2 1 100000\n"
          "##sequence-region   chr3 1 10000\n"
          "##sequence-region   chr4 1 1000\n"
          "chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1\n"
          "chr2\t.\tgene\t100\t200\t.\t+\t.\tID=gene2\n"
          "chr3\t.\tgene\t100\t200\t.\t+\t.\tID=gene3\n")
  out  = "\t".join(["chr4", ".", "locus", "1", "1000", ".", ".", ".",
                    "ID=locus1;fragment=true;unannot=true;"
                    "iLocus_type=iiLocus;effective_length=1000\n"])
  testout = run_parse_gff3(gff3)
  assert out == testout, "test 1 failed"

  gff3 = ("##gff-version\t3\n"
          "##sequence-region\tscaf00001   1   150000\n"
          "##sequence-region\tscaf00002   1   100000\n"
          "##sequence-region\tscaf00003   1   50000\n"
          "##sequence-region\tscaf00004   1   10000\n"
          "scaf00003\t.\tgene\t10000\t12000\t.\t+\t.\tID=gene1\n")
  out  = "\t".join(["scaf00001", ".", "locus", "1", "150000", ".", ".", ".",
                    "ID=locus001;fragment=true;unannot=true;iLocus_type=iiLocus;"
                    "effective_length=150000\n"])
  out += "\t".join(["scaf00002", ".", "locus", "1", "100000", ".", ".", ".",
                    "ID=locus002;fragment=true;unannot=true;iLocus_type=iiLocus;"
                    "effective_length=100000\n"])
  out += "\t".join(["scaf00004", ".", "locus", "1", "10000", ".", ".", ".",
                    "ID=locus003;fragment=true;unannot=true;iLocus_type=iiLocus;"
                    "effective_length=10000\n"])
  testout = run_parse_gff3(gff3, idformat="locus%03d")
  assert out == testout, "test 2 failed"

test_parse_gff3()

if __name__ == "__main__":
  desc = "Report iLoci for unannotated sequences"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument("--idfmt", type=str, default="locus%d",
                      help="An ID with a serial number is assigned to each "+
                      "locus; default format is 'locus%%d'")
  parser.add_argument("--counter", type=int, default=1,
                      help="Serial number to assign to first iLocus; default"+
                      " is 1")
  parser.add_argument("--src", type=str, default=".", help="Source label to "+
                      "use for GFF3 output (column 2); default is '.'")
  parser.add_argument("infile", help="Input data; use - to read from stdin")
  args = parser.parse_args()

  if args.infile == "-":
    fp = sys.stdin
  else:
    fp = open(args.infile, "r")

  print "##gff-version   3"
  for seqid, seqlen in parse_gff3(fp):
    locusid = args.idfmt % args.counter
    args.counter += 1
    attrs = "ID=%s;fragment=true;unannot=true;iLocus_type=iiLocus;" % locusid
    attrs += "effective_length=%d" % seqlen
    fields = [seqid, args.src, "locus", "1",
              str(seqlen), ".", ".", ".", attrs]
    print "\t".join(fields)
