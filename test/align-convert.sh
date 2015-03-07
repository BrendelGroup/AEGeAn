#!/usr/bin/env bash
set -eo pipefail

echo "    AEGeAn::align-convert.py"
tempfile="align-convert.temp"

for label in example1 example5
do
  data/scripts/align-convert.py \
      < data/gff3/gaeval-${label}-in.gff3 \
      > $tempfile

  diff $tempfile data/gff3/gaeval-${label}-fixed.gff3 > /dev/null
  status=$?
  result="FAIL"
  if [[ $status == 0 ]]; then
    result="PASS"
  fi
  printf "        | %-36s | %s\n" "GAEVAL $label" $result
  rm $tempfile
done
