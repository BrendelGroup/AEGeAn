#!/usr/bin/env bash
set -eo pipefail

if [[ $1 == "memcheck" ]]; then
  memcheckcmd="valgrind --leak-check=full --show-reachable=yes --suppressions=data/misc/libpixman.supp --suppressions=data/misc/libpango.supp --error-exitcode=1"
fi
echo "    AEGeAn::GAEVAL"
tempfile="gaeval.temp"


$memcheckcmd \
bin/gaeval data/gff3/gaeval-stream-unit-test-1.gff3 \
           data/gff3/gaeval-stream-unit-test-1.gff3 \
    > $tempfile

diff $tempfile data/gff3/gaeval-stream-unit-test-1-out.gff3 > /dev/null
status=$?
result="FAIL"
if [[ $status == 0 ]]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "unit test 1" $result
rm $tempfile
