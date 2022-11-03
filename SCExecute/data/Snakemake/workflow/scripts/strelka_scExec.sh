#!/bin/bash
# How to run the script: strelka Germline
#
# Assumes bcftools, bgzip, tabix, minDepth, python2, configureStrelkaGermlineWorkflow.py (streka script) are on the path. 
#
# Usage: stelka_scExec.sh <scBAM>.bam <reference_genome>.fasta [ <depth> [ <threads> ] ]
#
# If <depth> is omitted, no minimum depth is required and no intervals applied.
# If <threads> is omitted, strelka is run single-threaded. 
# The <scBAM>.bam BAM file should be indexed.
#
# Example scExcute usage:
#      scExecute -r <pooled-scBAM>.bam \
#                -b <barcodes-acceptlist>.tsv \
#                -i \
#                -t <worker-threads> \
#                -B 200 \
#                -D <pooled-scBAM>_strelka \
#                -C 'sh strelka_scExec.sh {{}} <reference_genome>.fasta <depth>' \
#                -F '{{BAMBASE}}.{{BARCODE}}.bam' \
#                -O '{{BAMBASE}}.{{BARCODE}}.log'
#
#
#

fa=$2
DEPTH=${3:-0}
THREADS=${4:-1}

function intervals() {
  minDepth $1 $2
}

line=${1}
l=$(echo ${line} | gawk -F"/" '{print $NF}' | sed 's/.bam//g')
rundir=${l}.run

if [ "$DEPTH" -gt 0 ]; then
    intervals ${line} $DEPTH | awk '{print $1"\t"$2"\t"$3+1}' > ${l}.bed
    rm -f ${l}.bed.gz
    bgzip -f ${l}.bed
    tabix -f -p bed ${l}.bed.gz
    python2 configureStrelkaGermlineWorkflow.py --bam ${line} --referenceFasta ${fa} --rna --runDir ${rundir} --callRegion ${l}.bed.gz
else
    python2 configureStrelkaGermlineWorkflow.py --bam ${line} --referenceFasta ${fa} --rna --runDir ${rundir}
fi
python2 ${rundir}/runWorkflow.py -m local -j "$THREADS"
mv ${rundir}"/results/variants/variants.vcf.gz" ${rundir}"/results/variants/"${l}"_strelka.vcf.gz"
mv ${rundir}"/results/variants/"${l}"_strelka.vcf.gz" .
bcftools index ${l}"_strelka.vcf.gz"
bcftools filter \
    -i 'TYPE="snp" && FILTER="PASS" && MIN(FORMAT/DP)>=3' \
    -r1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y ${l}"_strelka.vcf.gz" -Ov -o ${l}"_strelka_filt.vcf"
rm -f ${l}.bam*
rm -rf ${rundir}
rm -f "${l}".bed*
rm -f "${l}"_strelka.vcf.gz*
