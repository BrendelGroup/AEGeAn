#!/usr/bin/env bash

if [[ $1 == "memcheck" ]]; then
  memcheckcmd="valgrind --leak-check=full --show-reachable=yes --suppressions=data/share/libpixman.supp --error-exitcode=1"
fi
export LD_LIBRARY_PATH=./src/genometools/lib
export DYLD_LIBRARY_PATH=./src/genometools/lib
echo "    AEGeAn::xtractore"
tempfile="xtract.temp"

$memcheckcmd \
bin/xtractore --type CDS \
              --outfile $tempfile \
              --width 80 \
              data/gff3/mrj.gff3 data/fasta/mrj.gdna.fa

diff $tempfile data/fasta/mrj.cds.fa > /dev/null
status=$?
result="FAIL"
if [[ $status == 0 ]]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "major royal jelly" $result
rm $tempfile
