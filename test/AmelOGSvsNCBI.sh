#!/usr/bin/env bash

if [[ $1 == "memcheck" ]]; then
  memcheckcmd="valgrind --leak-check=full --show-reachable=yes --suppressions=data/misc/libpixman.supp --error-exitcode=1"
fi
echo "    AEGeAn::ParsEval"
tempfile="pe.temp"

$memcheckcmd \
bin/parseval --refrlabel=OGS \
             --predlabel=NCBI \
             data/gff3/amel-ogs-g716.gff3 \
             data/gff3/amel-ncbi-g716.gff3 \
  | grep -v '^Started' \
  > $tempfile

grep -v '^Started' data/misc/amel-ogs-vs-ncbi-parseval.txt \
  > ${tempfile}.orig

diff $tempfile ${tempfile}.orig > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "Amel Group7.16 (text)" $result
rm $tempfile ${tempfile}.orig



$memcheckcmd \
bin/parseval --datashare=data/share/ \
             --outformat=html \
             --png --outfile=$tempfile \
             --overwrite \
             --refrlabel=OGS \
             --predlabel=NCBI \
             data/gff3/amel-ogs-g716.gff3 \
             data/gff3/amel-ncbi-g716.gff3

cp -r data/misc/amel-ogs-vs-ncbi-parseval-html ${tempfile}.orig
grep -v '^Started' ${tempfile}/index.html > ${tempfile}.idx && mv ${tempfile}.idx ${tempfile}/index.html
grep -v '^Started' ${tempfile}.orig/index.html > ${tempfile}.idx && mv ${tempfile}.idx ${tempfile}.orig/index.html
rm ${tempfile}/*/*.png ${tempfile}.orig/*/*.png

diff -r $tempfile ${tempfile}.orig > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "Amel Group7.16 (HTML)" $result
rm -r $tempfile ${tempfile}.orig
