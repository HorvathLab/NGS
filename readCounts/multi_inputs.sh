#!/bin/sh

# remove -F to get rid of extra diagnostic columns
# remove -f to use higher stringency read filtering
# remove -U to turn off uniqueness filtering

while read line
do
{
  python src/readCounts.py \
    -s "variants.txt" \
    -r $line".sorted.bam" \
    -o $line"_10rnaed.csv" \
    -m 10 \
    -F \
    -f \ 
    -U \
    -t 20
}
done < list
