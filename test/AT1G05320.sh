#!/usr/bin/env bash

if [ "$1" == "memcheck" ]; then
  memcheckcmd="valgrind --leak-check=full --show-reachable=yes --suppressions=data/share/libpixman.supp --error-exitcode=1"
fi

FAILURES=0
testname="AT1G05320"
echo "    GFF3 I/O: '$testname'"
for test in codons sansexons sansutrs withprot
do
  tempfile="${testname}-${test}-temp.gff3"
  $memcheckcmd bin/canon-gff3 -o $tempfile -s TAIR10 data/gff3/AT1G05320-${test}.gff3
  if [ $? != 0 ]; then
    $((FAILURES++))
  fi
  bin/parseval -s $tempfile data/gff3/AT1G05320.gff3 2> /dev/null | grep 'perfect matches\.' | grep '100.0%' > /dev/null
  status=$?
  result="PASS"
  if [ $status != 0 ]; then
    result="FAIL"
    $((FAILURES++))
  fi
  printf "        | %-36s | %s\n" $test $result
  rm $tempfile
done

exit $FAILURES
