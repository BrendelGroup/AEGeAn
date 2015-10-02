#!/usr/bin/env bash
set -eo pipefail

if [[ $1 == "memcheck" ]]; then
  memcheckcmd="valgrind --leak-check=full --show-reachable=yes --suppressions=data/misc/libpixman.supp --suppressions=data/misc/libpango.supp --error-exitcode=1"
fi
echo "    AEGeAn::pmrna"
tempfile="misc.temp"


$memcheckcmd \
bin/pmrna < data/gff3/amel-shal.gff3 > $tempfile

diff $tempfile data/gff3/amel-shal-pmrna.gff3 > /dev/null
status=$?
result="FAIL"
if [[ $status == 0 ]]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "A. mellifera Shal" $result
rm $tempfile


$memcheckcmd \
bin/pmrna < data/gff3/amel-dock.gff3 > $tempfile

diff $tempfile data/gff3/amel-dock-pmrna.gff3 > /dev/null
status=$?
result="FAIL"
if [[ $status == 0 ]]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "A. mellifera dock" $result
rm $tempfile


$memcheckcmd \
bin/pmrna < data/gff3/amel-pseudo.gff3 > $tempfile

diff $tempfile data/gff3/amel-pseudo-pmrna.gff3 > /dev/null
status=$?
result="FAIL"
if [[ $status == 0 ]]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "A. mellifera pseudo" $result
rm $tempfile



echo "    AEGeAn::tidygff3"
$memcheckcmd \
bin/tidygff3 < data/gff3/ador-except-in.gff3 > $tempfile

diff $tempfile data/gff3/ador-except-out.gff3 > /dev/null
status=$?
result="FAIL"
if [[ $status == 0 ]]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "A. dorsata exception" $result
rm $tempfile
