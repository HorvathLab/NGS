#!/bin/sh

# multi-inputs.sh "<file>"

#
# File <file> contains each of the input filenames without 
# and extension.
# 
# remove -F to get rid of extra diagnostic columns
# remove -f to use higher stringency read filtering
# remove -U to turn off uniqueness filtering
#

if [ "$1" = "" ]; then
    echo "Usage: multi-inputs.sh <file>" 1>&2
    exit 1;
fi

PYTHON="python3"
if [ -d ../src ]; then
  PROG="$PYTHON ../src/readCounts.py"
else
  PROG=../bin/readCounts
fi

LINE="$1"

while read line
do
{
  "$PROG" \
    -s "variants.txt" \
    -r $line".sorted.bam" \
    -o $line"_10rnaed.csv" \
    -m 10 \
    -F \
    -f \ 
    -U \
    -t 20
}
done < "$LINE"
