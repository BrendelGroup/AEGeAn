#!/usr/bin/env bash
set -eo pipefail

echo "    AEGeAn::align-convert.py"
tempfile="align-convert.temp"

data/scripts/align-convert.py \
    < data/gff3/gaeval-example1-in.gff3 \
    > $tempfile

diff $tempfile data/gff3/gaeval-example1-fixed.gff3 > /dev/null
status=$?
result="FAIL"
if [[ $status == 0 ]]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "GAEVAL example 1" $result
rm $tempfile
