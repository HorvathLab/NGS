#!/bin/sh
if [ -d ../src ]; then
  PROG=../src/RNA2DNAlign.py
else
  PROG=../bin/RNA2DNAlign
fi
rm -rf example-output
DARNEDFILE="DARNED_hg19.txt"
DARNED=""
if [ -f $DARNEDFILE ]; then
  DARNED="-d $DARNEDFILE"
fi
COSMICFILE="CosmicMutantExport_hg19.tsv.gz"
COSMIC=""
if [ -f $COSMICFILE ]; then
  COSMIC="-c $COSMICFILE"
fi
$PROG -r "example-*.bam" -s "example-*.vcf" -m 3 -e UCSC_Human_hg19_RefSeq_CDS_exon_coordinates.txt $DARNED $COSMIC -o example-output
