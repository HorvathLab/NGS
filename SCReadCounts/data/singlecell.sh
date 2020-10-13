#!/bin/sh

# Clean up previous runs...
rm -f singlecell-output.tsv singlecell-output.cnt.matrix.tsv singlecell-output.vaf.matrix.tsv

#
# Equivalent to running these three commands...
#
# readCounts -s "singlecell_222_5_chr17.txt" -r "singlecell_chr17.bam" -m 0 -G "UMITools" -o "singlecell-output.tsv"
# readCountsMatrix -c "singlecell-output.tsv" -M Ref:Var -o "singlecell-output.cnt.matrix.tsv"                   
# readCountsMatrix -c "singlecell-output.tsv" -M VAF -m 3 -o "singlecell-output.vaf.matrix.tsv"
#
# readCountsMatrix -c "singlecell-output.tsv" -M VAF -m 5 -o "singlecell-output.vaf.matrix.tsv"
#

set -x
../bin/scReadCounts -r singlecell_chr17.bam -s singlecell_222_5_chr17.txt -m 3 -o singlecell-output.tsv

# Regenerate only the VAF matrix
rm -f singlecell-output.vaf.matrix.tsv
../bin/scReadCounts -r singlecell_chr17.bam -s singlecell_222_5_chr17.txt -m 5 -o singlecell-output.tsv

