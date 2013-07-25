#!/usr/bin/env bash

testname="FBgn0035002"
echo "    GFF3 I/O: '$testname'"
for test in codons sansexons sansutrs
do
  tempfile="${testname}-${test}-temp.gff3"
  bin/canon-gff3 -o $tempfile -s FlyBase data/gff3/FBgn0035002-${test}.gff3
  diff $tempfile data/gff3/FBgn0035002.gff3 >/dev/null || \
       diff $tempfile data/gff3/FBgn0035002-altsort1.gff3 >/dev/null || \
       diff $tempfile data/gff3/FBgn0035002-altsort2.gff3 >/dev/null
  status=$?
  result="FAIL"
  if [ $status == 0 ]; then
    result="PASS"
  fi
  printf "        | %-30s| %s\n" $test $result
  rm $tempfile
done
