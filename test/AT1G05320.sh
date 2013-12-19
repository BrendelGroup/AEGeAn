#!/usr/bin/env bash

testname="AT1G05320"
echo "    GFF3 I/O: '$testname'"
for test in codons sansexons sansutrs withprot
do
  tempfile="${testname}-${test}-temp.gff3"
  bin/canon-gff3 -o $tempfile -s TAIR10 data/gff3/AT1G05320-${test}.gff3
  bin/parseval -s $tempfile data/gff3/AT1G05320.gff3 2> /dev/null | grep 'perfect matches\.' | grep '100.0%' > /dev/null
  status=$?
  result="FAIL"
  if [ $status == 0 ]; then
    result="PASS"
  fi
  printf "        | %-36s | %s\n" $test $result
  rm $tempfile
done
