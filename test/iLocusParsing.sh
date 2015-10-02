#!/usr/bin/env bash
set -eo pipefail

if [ "$1" == "memcheck" ]; then
  memcheckcmd="valgrind --leak-check=full --show-reachable=yes --suppressions=data/misc/libpixman.supp --suppressions=data/misc/libpango.supp --error-exitcode=1"
fi
tempfile="ilocus-tempfilefile.gff3"
failures=0

run_func_test()
{
  testlabel="$1"
  shift
  testoutfile="$1"
  shift
  filelabel=$(echo "$testlabel" | tr -d '(' | tr -d ')' | tr ' ' '_')-test.gff3

  $memcheckcmd bin/locuspocus --retainids $@ > /dev/null 2>&1
  if [ $? != 0 ]; then
    echo "Error running functional test '$testlabel'"
    exit 1
  fi
  
  set +e
  diff $tempfile $testoutfile > /dev/null 2>&1
  status=$?
  set -e
  
  result="FAIL"
  if [ $status == 0 ]; then
    result="PASS"
    rm -f $tempfile $filelabel
  else
    mv $tempfile $filelabel
    failures=$((failures+1))
  fi
  printf "        | %-36s | %s\n" "$testlabel" $result
}

echo "    AEGeAn::LocusPocus"

run_func_test "default" data/gff3/ilocus.out.noskipends.gff3 --delta=200 --outfile=${tempfile} --parent mRNA:gene data/gff3/ilocus.in.gff3
run_func_test "end skip" data/gff3/ilocus.out.skipends.gff3 --delta=200 --outfile=${tempfile} --skipends --parent mRNA:gene data/gff3/ilocus.in.gff3
run_func_test "Apis mellifera plap" data/gff3/amel-plap-out.gff3 --outfile=${tempfile} --skipends data/gff3/amel-plap.gff3
run_func_test "Apis mellifera plap (CDS)" data/gff3/amel-plap-out-cds.gff3 --outfile=${tempfile} --skipends --cds data/gff3/amel-plap.gff3
run_func_test "Apis mellifera LSM" data/gff3/amel-lsm-out-cds.gff3 --outfile=${tempfile} --skipends --cds data/gff3/amel-lsm.gff3
run_func_test "Megachile rotundata CST (intron)" data/gff3/mrot-cst-out-cds.gff3 --outfile=${tempfile} --skipends --cds data/gff3/mrot-cst.gff3
run_func_test "iiLocus lengths (Amel OGS Group7.16)" data/misc/amel-ogs-ilens.txt --delta=300 --ilens=${tempfile} --cds data/gff3/amel-ogs-g716.gff3
run_func_test "Nasonia vitripennis (intron gene)" data/gff3/nvit-exospindle-out.gff3 --outfile=${tempfile} --cds data/gff3/nvit-exospindle.gff3
run_func_test "A. echinatior (intron gene + ncRNA)" data/gff3/aech-dachsous-out.gff3 --outfile=${tempfile} --cds data/gff3/aech-dachsous.gff3

exit $failures
