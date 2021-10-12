#!/bin/sh

# Clean up previous runs...
rm -f singlecell-output.{tsv,*.matrix.tsv}
rm -f singlecell2-output.{tsv,*.matrix.tsv}
rm -f singlecell2-all-output.{tsv,*.matrix.tsv}

set -x

#
# Equivalent to running these three commands...
#
# UMI-tools is the default cell-barcode strategy
#
# readCounts -s "singlecell_222_5_chr17.txt" -r "singlecell_chr17.bam" -m 0 -G "UMI-tools" -o "singlecell-output.tsv"
# readCountsMatrix -c "singlecell-output.tsv" -M Ref:Var -o "singlecell-output.cnt.matrix.tsv"                   
# readCountsMatrix -c "singlecell-output.tsv" -M VAF -m 3 -o "singlecell-output.vaf-m3.matrix.tsv"
#
# readCountsMatrix -c "singlecell-output.tsv" -M VAF -m 5 -o "singlecell-output.vaf-m5.matrix.tsv"
#
../bin/scReadCounts -r singlecell_chr17.bam -s singlecell_222_5_chr17.txt -m 3 -o singlecell-output.tsv

#
# Regenerate the VAF matrix (different output name: singlecell-output.vaf-m5.matrix.tsv)
#
../bin/scReadCounts -r singlecell_chr17.bam -s singlecell_222_5_chr17.txt -m 5 -o singlecell-output.tsv

#
# STARsolo example
#
../bin/scReadCounts -r singlecell2_117.bam -s singlecell2_117_snvs.txt -m 5 \
                    -G STARsolo -b singlecell2_117_barcodes.tsv -o singlecell2-output.tsv

#
# Override accept list for barcodes, accept all barcodes
#
../bin/scReadCounts -r singlecell2_117.bam -s singlecell2_117_snvs.txt -m 5 \
                    -G STARsolo -b None -o singlecell2-all-output.tsv

