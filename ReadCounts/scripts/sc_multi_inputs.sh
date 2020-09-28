#!/bin/sh
# 
# sc_multi_inputs.sh "<variants-file>" "<inputs-file>"
# 
#
# File <variants-file> contains loci of interest in VCF format. 
# 
# File <inputs-file> contains each of the input filenames without 
# path and extension.
# 
# select read filtering using -f option, see filter.ini.
#

if [ "$2" = "" ]; then
    echo "Usage: sc_multi_inputs.sh <variants-file> <inputs-file>" 1>&2
    exit 1;
fi

if [ -d ../src ]; then
  PROG="${PYTHON:=python3} ../src/scReadCounts.py"
else
  PROG=../bin/scReadCounts
fi

VARIANTS="$1"
INPUTS="$2"

while read line
do
{
  $PROG \
    -s "$VARIANTS" \
    -r $line".bam" \
    -o $line".tsv" \
    -f "Basic" \
    -m 5
}
done < "$INPUTS"
