#!/usr/bin/env bash

if [ "$1" == "memcheck" ]; then
  memcheckcmd="valgrind --leak-check=full --show-reachable=yes --suppressions=data/share/libpixman.supp --error-exitcode=1"
fi

testname="FBgn0035002"
echo "    GFF3 I/O: '$testname'"
for test in codons sansexons sansutrs
do
  tempfile="${testname}-${test}-temp.gff3"
  $memcheckcmd bin/canon-gff3 -o $tempfile -s FlyBase data/gff3/FBgn0035002-${test}.gff3
  if [ $? != 0 ]; then
    exit 1
  fi
  bin/parseval -s $tempfile data/gff3/FBgn0035002.gff3 2> /dev/null | grep 'perfect matches\.' | grep '100.0%' > /dev/null
  status=$?
  result="FAIL"
  if [ $status == 0 ]; then
    result="PASS"
  fi
  printf "        | %-36s | %s\n" $test $result
  rm $tempfile
done
