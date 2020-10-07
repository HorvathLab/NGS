#!/bin/sh

# Clean up previous runs
rm -rf example-output.tsv

# Execute read counts
set -x
../bin/readCounts -r example-NRNA.bam -s example-NRNA.vcf -o example-output.tsv
