# readCount Output 

readCounts output file are created inside of the directory specified by the user.

## readCounts output file

A tab-separated values file consisting of the computed read-counts. This file contains the read counts
for each SNV locus in each BAM file and computes the various
statistical tests below above:

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

HomoVarSc
> Score of locus as homozygous variant.

HetSc
> Score of locus as heterozygous reference and variant.

HomoRefSc
> Score of locus as homozygous reference.

VarDomSc
> Score of locus as dominant for the variant allele.

RefDomSc
> Score of locus as dominant for the reference allele.

NotHomoVarpV
> p-Value of read counts with respect to homozygous variant null model.

NotHomoRefpV
> p-Value of read counts with respect to homozygous reference null model.

NotHetpV
> p-Value of read counts with respect to heterozygous reference and variant null model.

VarDompV
> p-Value of increased variant read counts with respect to heterozygous reference and variant null model.

RefDompV
> p-Value of increased reference read counts with respect to heterozygous reference and variant null model.

NotHomoVarFDR
> Multiple-test corrected FDR significance of read counts with respect to homozygous variant null model.

NotHomoRefFDR
> Multiple-test corrected FDR significance of read counts with respect to homozygous reference null model.

NotHetFDR
> Multiple-test corrected FDR significance of read counts with respect to heterozygous reference and variant null model.

VarDomFDR
> Multiple-test corrected FDR significance of increased variant read counts with respect to heterozygous reference and variant null model.

RefDomFDR
> Multiple-test corrected FDR significance of increased reference read counts with respect to heterozygous reference and variant null model.
