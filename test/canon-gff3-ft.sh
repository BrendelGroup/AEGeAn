#!/usr/bin/env bash
set -eo pipefail

if [[ $1 == "memcheck" ]]; then
  memcheckcmd="valgrind --leak-check=full --show-reachable=yes --suppressions=data/misc/libpixman.supp --suppressions=data/misc/libpango.supp --error-exitcode=1"
fi
echo "    AEGeAn::CanonGFF3"
tempfile="canon.temp"

$memcheckcmd \
bin/canon-gff3 --outfile $tempfile data/gff3/amel-gene-multitrans.gff3

diff $tempfile data/gff3/amel-gene-multitrans-canon.gff3 > /dev/null
status=$?
result="FAIL"
if [[ $status == 0 ]]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "A. mellifera gene multitrans" $result
rm $tempfile
