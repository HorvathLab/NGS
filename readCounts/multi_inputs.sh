#!/bin/sh
while read line
do
{
  python src/readCounts.py \
    -s "variants.txt" \
    -r $line".sorted.bam" \
    -o $line"_10rnaed.csv" \
    -m 10 \
    -F False \
    -f False \
    -U False \
    -t 20
}
done < list
