#!/bin/sh
# set -x
PYTHON="python3"
if [ -d ../src ]; then
  PROG="$PYTHON ../src/RNA2DNAlign.py"
else
  PROG=../bin/RNA2DNAlign
fi
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
rm -rf example-output
set -x
$PROG -r "example-*.bam" -s "example-*.vcf" -m 3 -e UCSC_Human_hg19_RefSeq_CDS_exon_coordinates.txt $DARNED $COSMIC -o example-output
