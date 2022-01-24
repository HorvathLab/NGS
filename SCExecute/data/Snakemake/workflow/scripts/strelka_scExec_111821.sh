#!/bin/bash
# How to run the script:

# Germline: sh strelka.sh Germline bam_file 

fa='/data/lab/Homo_sapiens/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa'

ml strelka
ml bcftools
ml python/2.7.6
if [[ "$1" == "Germline" ]]
then
	line=${2}
	l=$(echo ${line} | gawk -F"/" '{print $NF}' | sed 's/.bam//g')
	configureStrelkaGermlineWorkflow.py --bam ${line} --referenceFasta ${fa} --rna --runDir "_"${l}
	python "_"${l}/runWorkflow.py -m local &> log_${l}
	mv "_"${l}"/results/variants/variants.vcf.gz" "_"${l}"/results/variants/"${l}"_strelka.vcf.gz"
	mv "_"${l}"/results/variants/"${l}"_strelka.vcf.gz" .
	bcftools index ${l}"_strelka.vcf.gz"
	bcftools filter \
           -i 'TYPE="snp" && FILTER="PASS"' \
           -r1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y ${l}"_strelka.vcf.gz" -Ov -o ${l}"_strelka_filt.vcf"
	
else
	echo "Strelka does not have this option"
fi
