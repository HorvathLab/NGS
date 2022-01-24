#!/bin/bash
#For RNA-Seq, preprocess the RNA-Seq data 
#splits reads into exon segments (getting rid of Ns but maintaining grouping information) 
#hard-clip any sequences overhanging into the intronic regions

ml gatk
ml samtools
ml bcftools

line=${1}
pref=$(echo ${line} | gawk -F"/" '{print $NF}' | sed 's/.bam//g')
#echo $pref
gatk AddOrReplaceReadGroups \
       I=$line \
       O=$pref"_rg.bam" \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20

samtools index -@ 10 $pref"_rg.bam"

gatk SplitNCigarReads \
	-R /data/lab/Homo_sapiens/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	-I $pref"_rg.bam" \
	-O $pref"_split_rg.bam" \
	-RF MappingQualityReadFilter

samtools index -@ 10 $pref"_split_rg.bam"

gatk --java-options "-Xmx4g" HaplotypeCaller \
	-R /data/lab/Homo_sapiens/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	-I $pref"_split_rg.bam" \
	-O $pref"_gatk.vcf"

bcftools view $pref"_gatk.vcf" -Oz -o $pref"_gatk.vcf.gz"

bcftools index $pref"_gatk.vcf.gz"

bcftools filter \
	-i 'TYPE="snp" && MIN(INFO/DP)>=3 && QUAL>=100 && MQ==60' \
	-r1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $pref"_gatk.vcf.gz" -Ov -o $pref"_gatk_filt.vcf"

rm $pref"_rg.bam"
rm $pref"_rg.bam.bai"
rm $line"_split_rg.bam"
rm $pref"_split_rg.bai"
rm $pref"_split_rg.bam.bai"
rm $pref"_gatk.vcf.idx"
rm $pref"_gatk.vcf.gz"
rm $pref"_gatk.vcf.gz.csi"
echo "DONE"
