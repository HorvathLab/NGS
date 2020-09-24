#!/bin/sh

# scReadCounts.sh <reads>.bam <readCounts options>
#
# Extract read counts for various alleles, and then form VAF and counts
# matrix for cell barcodes added to the BAM file by UMITools.
# 

# Minimum # of reads for reporting VAF in VAF matrix
MINREADS=10

if [ "$1" = "" ]; then
    echo "Usage: scReadCounts.sh <reads>.bam <read counts options>" 1>&2
    exit 1;
fi

# Figure out if we are using the Python or binary distribution
if [ -d ../src ]; then
  PROG="${PYTHON3:-python3} ../src/readCounts.py"
  PROG1="${PYTHON3:-python3} ../src/readCountsMatrix.py"
else
  PROG=../bin/readCounts
  PROG1=../bin/readCountsMatrix
fi

BAM="$1"
shift

BASE=`basename "$BAM" .bam`

set -x
$PROG -r "$BAM" $@ -G UMITools -o "${BASE}_counts.tsv" && \
$PROG1 -c "${BASE}_counts.tsv" -M "VAF" -m $MINREADS -o "${BASE}_VAF_matrix.tsv" && \
$PROG1 -c "${BASE}_counts.tsv" -M "Ref:Var" -o "${BASE}_counts_matrix.tsv"

