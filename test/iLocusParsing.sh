#!/usr/bin/env bash

echo "    iLocus Parsing"

temp="iLocusParseTest.gff3"
bin/locuspocus --intloci --delta=200 --outfile=${temp} --parent mRNA:gene data/gff3/ilocus.in.gff3 > /dev/null 2>&1
diff ${temp} data/gff3/ilocus.out.noskipends.gff3 > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "default" $result
rm ${temp}

bin/locuspocus --intloci --delta=200 --outfile=${temp} --skipends --parent mRNA:gene data/gff3/ilocus.in.gff3 > /dev/null 2>&1
diff ${temp} data/gff3/ilocus.out.skipends.gff3 > /dev/null 2>&1
status=$?
result="FAIL"
if [ $status == 0 ]; then
  result="PASS"
fi
printf "        | %-36s | %s\n" "end skip" $result
rm ${temp}

