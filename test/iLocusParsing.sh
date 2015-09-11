#!/usr/bin/env bash
set -eo pipefail

if [ "$1" == "memcheck" ]; then
  memcheckcmd="valgrind --leak-check=full --show-reachable=yes --suppressions=data/misc/libpixman.supp --suppressions=data/misc/libpango.supp --error-exitcode=1"
fi

echo "    AEGeAn::LocusPocus"

temp="iLocusParseTest.gff3"
$memcheckcmd bin/locuspocus --delta=200 --outfile=${temp} --parent mRNA:gene data/gff3/ilocus.in.gff3 > /dev/null 2>&1
if [ $? != 0 ]; then
  exit 1
fi
diff ${temp} data/gff3/ilocus.out.noskipends.gff3 > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "default" $result
rm ${temp}

$memcheckcmd bin/locuspocus --delta=200 --outfile=${temp} --skipends --parent mRNA:gene data/gff3/ilocus.in.gff3 > /dev/null 2>&1
if [ $? != 0 ]; then
  exit 1
fi
diff ${temp} data/gff3/ilocus.out.skipends.gff3 > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "end skip" $result
rm ${temp}

$memcheckcmd bin/locuspocus --outfile=${temp} --skipends data/gff3/amel-plap.gff3 > /dev/null 2>&1
if [ $? != 0 ]; then
  exit 1
fi
diff ${temp} data/gff3/amel-plap-out.gff3 > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "Apis mellifera plap" $result
rm ${temp}

$memcheckcmd bin/locuspocus --outfile=${temp} --skipends --cds data/gff3/amel-plap.gff3 > /dev/null 2>&1
if [ $? != 0 ]; then
  exit 1
fi
diff ${temp} data/gff3/amel-plap-out-cds.gff3 > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "Apis mellifera plap (CDS)" $result
rm ${temp}

$memcheckcmd bin/locuspocus --outfile=${temp} --skipends --cds data/gff3/amel-lsm.gff3 > /dev/null 2>&1
if [ $? != 0 ]; then
  exit 1
fi
diff ${temp} data/gff3/amel-lsm-out-cds.gff3 > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "Apis mellifera LSM" $result
rm ${temp}

$memcheckcmd bin/locuspocus --outfile=${temp} --skipends --cds data/gff3/mrot-cst.gff3 > /dev/null 2>&1
if [ $? != 0 ]; then
  exit 1
fi
diff ${temp} data/gff3/mrot-cst-out-cds.gff3 > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "Megachile rotundata CST (intron)" $result
rm ${temp}

$memcheckcmd bin/locuspocus --delta=300 --ilens=${temp} --cds data/gff3/amel-ogs-g716.gff3 > /dev/null 2>&1
if [ $? != 0 ]; then
  exit 1
fi
diff ${temp} data/misc/amel-ogs-ilens.txt > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "iiLocus lengths (Amel OGS Group7.16)" $result
rm ${temp}
