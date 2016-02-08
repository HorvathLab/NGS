#!/bin/sh
if [ -d ../src ]; then
  PROG=../src/RNA2DNA.py
else
  PROG=../src/RNA2DNA
fi
rm -rf example-output
DARNEDFILE=../nodist/hg19.txt
DARNED=
if [ -f $DARNEDFILE ]; then
  DARNED="-d $DARNEDFILE"
fi
COSMICFILE=../nodist/CosmicMutantExport.tsv.gz
COSMIC=
if [ -f $COSMICFILE ]; then
  COSMIC="-c $COSMICFILE"
fi
$PROG -r "example-*.bam" -s "example-*.vcf" -m 3 -e UCSC_Human_hg19_RefSeq_CDS_exon_coordinates.txt $DARNED $COSMIC -o example-output
