# SCReadCount Output 

SCReadCounts output files are created based on the output filename provided. 

## SCReadCounts output file

A tab-separated values file consisting of the computed read-counts. This file contains the read counts
for each SNV locus in each BAM file.

CHROM
> Chromosome identifier

POS
> Chromosome position of the variant

REF
> Reference allele nucleotide

ALT
> Variant allele nucleotide

AlignedReads
> Name of the aligned reads file for the following read counts.

SNPCountForward
> Number of forward oriented variant reads in the paired end alignment.

SNPCountReverse
> Number of reverse oriented variant reads in the paired end alignment.

RefCountForward
> Number of forward oriented reference reads in the paired end alignment.

RefCountReverse
> Number of reverse oriented reference reads in the paired end alignment.

SNPCount
> Total number of variant reads.

RefCount
> Total number of reference reads.

R
> Proportion of variant reads.

