#!/bin/sh
# set -x
PYTHON="python3"
if [ -d ../src ]; then
  PROG="$PYTHON ../src/readCounts.py"
else
  PROG=../bin/readCounts
fi
rm -rf example-output.tsv
set -x
$PROG -r example-NRNA.bam -s example-NRNA.vcf -o example-output.tsv
