#!/usr/bin/env bash

testname="FBgn0035002"
echo "Test '$testname'"
for test in codons sansexons sansutrs
do
  tempfile="${testname}-${test}-temp.gff3"
  bin/canon-gff3 -o $tempfile -s FlyBase data/gff3/FBgn0035002-${test}.gff3
  diff $tempfile data/gff3/FBgn0035002.gff3 >/dev/null || \
       diff $tempfile data/gff3/FBgn0035002-altsort1.gff3 >/dev/null || \
       diff $tempfile data/gff3/FBgn0035002-altsort2.gff3 >/dev/null
  status=$?
  result="fail"
  if [ $status == 0 ]; then
    result="pass"
  fi
  printf "%12s: %s\n" $test $result
  rm $tempfile
done
