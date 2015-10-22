#!/bin/sh

# set -x

#
# Assumes a directory structure like:
# 
#  sample1/junctions1.bed
#          reads.bam
#          reads.bam.bai
#          snps.vcf
#  sample2/junctions1.bed
#          reads.bam
#          reads.bam.bai
#          snps.vcf
#  sample3/junctions1.bed
#          reads.bam
#          reads.bam.bai
#          snps.vcf
#  sample4/junctions1.bed
#          reads.bam
#          reads.bam.bai
#          snps.vcf
 
#
# Make a directory containing only the sample directories and cd to it. 
# Only the file extensions matter, vcf files may be as output from
# Samtools or Samtools + SeattleSeq
#

COMMON="-R 5 -U -d 300"
SNPLICE=`dirname $0`
EXT=".sh"
if [ \! -f "$SNPLICE/SNPlice$EXT" ]; then
  EXT=".py"
fi
for b in */*.bam; do
  s=`dirname "$b"`

  # Potentially, it would be appropriate to use all junction and snp
  # files for sample specific read analysis...
  # $SNPLICE/SNPlice -s "*/*.vcf" -r "$s/*.bam" -j "*/*.bed" -o "$a/snplice.tsv" $COMMON

  # However, for the time being, assume we only consider the SNPs and
  # junctions in each sample directory
  echo $SNPLICE/SNPlice$EXT -s "$s/*.vcf" -r "$s/*.bam" -j "$s/*.bed" -o "$s/snplice.tsv" $COMMON
  $SNPLICE/SNPlice$EXT -s "$s/*.vcf" -r "$s/*.bam" -j "$s/*.bed" -o "$s/snplice.tsv" $COMMON

done

# If the first execution style is used, then we can just use
# "SNPlice-Combine" instead of SNPlice all over again
# $SNPLICE/SNPlice-Combine$EXT -c "*/snplice.tsv" -o snplice.tsv

# However, for the time being, we cannot assume that all SNP loci and
# junctions have been counted in the sample directories
echo $SNPLICE/SNPlice$EXT -s "*/*.vcf" -r "*/*.bam" -j "*/*.bed" -o snplice.tsv $COMMON
$SNPLICE/SNPlice$EXT -s "*/*.vcf" -r "*/*.bam" -j "*/*.bed" -o snplice.tsv $COMMON
