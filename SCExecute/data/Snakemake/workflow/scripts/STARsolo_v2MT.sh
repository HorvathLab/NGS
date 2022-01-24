#!/bin/bash
#100 hours time limit
#SBATCH --time 100:00:00
#short,32 processors (two nodes worth)
#SBATCH -p defq -n 16


# File prefixed with "_wasp_Aligned.sortedByCoord.out.bam" is the original unfiltered bam file
# File called "*_wasp_Aligned.sortedByCoord_vW_filt.bam" is filtered, sorted and indexed bam file

GENOME_DIR="/data/lab/Homo_sapiens/genome"
GTF_FILE="/data/lab/Homo_sapiens/genome/Homo_sapiens.GRCh38.79.gtf"

module load star
module load samtools
ml cellranger

list=${1}
if [ ! -d STARsolo_wasp ]; then mkdir STARsolo_wasp; fi

while read line
do
{
STAR --genomeDir ${GENOME_DIR} --sjdbGTFfile ${GTF_FILE} --sjdbOverhang 150 --twopassMode Basic --readFilesIn ${line}"_2.fastq" ${line}"_1.fastq" --soloType CB_UMI_Simple --soloCBwhitelist /data/lab/737K-august-2016.txt --outSAMattributes NH HI AS nM NM MD CB UB CR UR GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --runThreadN 10 --outFileNamePrefix STAR_wasp_bams_filt2/${line}"_wasp_" --soloFeatures Gene GeneFull SJ Velocyto
# Store header file temporarily
samtools view -H STAR_wasp_bams_filt2/${line}"_wasp_Aligned.sortedByCoord.out.bam" > header_${line}".sam"
# extract only reads that have vW:i:1 or no vW tag
samtools view STAR_wasp_bams_filt2/${line}"_wasp_Aligned.sortedByCoord.out.bam" | awk -F "\t" '$0 !~ /vW:i:[2-4]/{print $0}' > temp_aa.sam
# sort and index to get a bam and bai file
cat header_${line}".sam" temp_aa.sam > ${line}"_wasp_vW_filt.sam"
samtools sort -T temp_aa -o STAR_wasp_bams_filt2/${line}"_wasp_Aligned.sortedByCoord_vW_filt.bam" ${line}"_wasp_vW_filt.sam"
samtools index STAR_wasp_bams_filt2/${line}"_wasp_Aligned.sortedByCoord_vW_filt.bam"
# remove the temporary files
rm temp_aa.sam
rm header_${line}".sam"
rm ${line}"_wasp_vW_filt.sam"
rm STAR_wasp_bams_filt2/${line}"_wasp_SJ.out.tab"
rm –r STAR_wasp_bams_filt2/${line}"_wasp__STARgenome"
rm –r STAR_wasp_bams_filt2/${line}"_wasp__STARpass1"
}
done < ${1}
echo "DONE"

