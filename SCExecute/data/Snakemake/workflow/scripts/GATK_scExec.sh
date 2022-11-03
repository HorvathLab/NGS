#!/bin/bash
#
# Assumes samtools, bcftools, gatk3, minDepth are on the path
# Java virtual machine options can be set using the JAVAOPTS environment varialble.
# If unset -Xmx2048M (2G memory) is used. 
# 
# Usage: GATK_scExec.sh <scBAM>.bam <reference_genome>.fasta [ <depth> ]
#
# If <depth> is omitted, no minimum depth is required and no intervals applied.
# The <scBAM>.bam BAM file should be indexed.
#
# Example scExecute usage: 
#      scExecute -r <pooled-scBAM>.bam \
#                -b <barcodes-acceptlist>.tsv \
#                -i \
#                -t <worker-threads> \
#                --cpuaffinity \
#                -B 200 \
#                -D <pooled-scBAM>_gatk \
#                -F '{{BAMBASE}}.{{BARCODE}}.bam' \
#                -C 'sh GATK_scExec.sh {{}} <reference_genome>.fasta <depth>' \
#                -O '{{BAMBASE}}.{{BARCODE}}.log'
#
#

function intervals() {
  minDepth $1 $2
}

DEPTH=${3:-0}
line=${1}
pref=$(echo ${line} | gawk -F"/" '{print $NF}' | sed 's/.bam//g')

gatk3 AddOrReplaceReadGroups \
       --java-options "${JAVAOPTS:--Xmx2048M}" \
       I=$line \
       O=$pref"_rg.bam" \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20

samtools index -@ 10 $pref"_rg.bam"

gatk3 SplitNCigarReads \
        --java-options "${JAVAOPTS:--Xmx2048M}" \
	-R "$2" \
	-I $pref"_rg.bam" \
	-O $pref"_split_rg.bam" \
	-RF MappingQualityReadFilter

samtools index -@ 10 $pref"_split_rg.bam"

if [ "$DEPTH" -gt 0 ]; then

  intervals $pref"_split_rg.bam" $DEPTH | \
        awk '{printf("%s:%s-%s\n",$1,$2+1,(($2+1)>$3)?($2+1):$3);}' \
        > $pref"_split_rg.d${DEPTH}.list"

  gatk3 HaplotypeCaller \
        --java-options "${JAVAOPTS:--Xmx2048M}" \
	-R "$2" \
	-I $pref"_split_rg.bam" \
	-O $pref"_gatk.vcf" \
        -ip 100 \
        -L $pref"_split_rg.d${DEPTH}.list"

else

  gatk3 HaplotypeCaller \
        --java-options "${JAVAOPTS:--Xmx2048M}" \
	-R "$2" \
	-I $pref"_split_rg.bam" \
	-O $pref"_gatk.vcf"

fi

bcftools view $pref"_gatk.vcf" -Oz -o $pref"_gatk.vcf.gz"

bcftools index $pref"_gatk.vcf.gz"

bcftools filter \
	-i 'TYPE="snp" && MIN(INFO/DP)>=3 && QUAL>=100 && MQ==60' \
	-r1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $pref"_gatk.vcf.gz" -Ov -o $pref"_gatk_filt.vcf"

rm $pref".bam"
rm $pref".bam.bai"
rm $pref"_rg.bam"
rm $pref"_rg.bam.bai"
rm $pref"_split_rg.bam"
rm $pref"_split_rg.bai"
rm $pref"_split_rg.bam.bai"
rm $pref"_split_rg.d"*".list"
rm $pref"_gatk.vcf.idx"
rm $pref"_gatk.vcf.gz"
rm $pref"_gatk.vcf.gz.csi"
echo "DONE"
