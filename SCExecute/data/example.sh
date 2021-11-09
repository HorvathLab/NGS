#!/bin/sh

set -x
../bin/scExecute -r singlecell2_117.bam -G STARsolo -i -C "samtools view {} | tail > test_{BAMBASE}_{BARCODE}_{CBINDEX}.sam" -t 3
