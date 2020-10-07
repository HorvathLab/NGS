#!/bin/sh
# 
# sc_multi_inputs.sh "<inputs-file>"
#
# File <inputs-file> contains each of the read alignments (BAM) files
# without the .bam extension
# 

if [ "$2" = "" ]; then
    echo "Usage: sc_multi_inputs.sh <inputs-file>" 1>&2
    exit 1;
fi

INPUTS="$1"

while read line
do
{
  scReadCounts \
    -s "variants.txt" \
    -r $line".bam" \
    -o $line".tsv" \
    -f "Basic" \
    -m 5
}
done < "$INPUTS"
