# ReadCounts Output 

A tab-separated values file consisting of the computed read-counts. This file contains the read counts
for each SNV locus in each read group (by default, each BAM file):

CHROM
> Chromosome identifier

POS
> Chromosome position of the variant

REF
> Reference allele nucleotide

ALT
> Variant allele nucleotide

ReadGroup
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

GoodReads
> Total number of reads considered.

%BadReads
> Proportion of reads not counted.

VAF
> Proportion of variant reads.

## Extended Output Columns

### Genotype Likelihood

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

### Read filtering statistics

OtherCountForward
> Reads with unexpected nucleotides in the forward oriented aligned read.
OtherCountReverse    
> Reads with unexpected nucleotides in the reverse oriented aligned read.
OtherCount 
> Reads with unexpected nucleotides.
FilteredSNVLociReads
> Total number of reads that passed the read-filters.
SNVLociReads    
> Total number of reads covering the site.
Alignment:IsDuplicate
> Alignment has the duplicate flag. 
Alignment:IsQCFail      
> Alignment has the QCFail flag.
Alignment:IsSecondary   
> Alignment has the secondary flag.
Alignment:IsUnmapped    
> Alignment has the unmapped flag.
BadCIGAROperation
> Unexpected CIGAR string operation.
BaseQualityTooLow
> Base quality at locus is too low.
DuplicateRead
> Duplicate read
GapInQueryAtSNVLocus
> Gap in the alignment at the site of interest.
MappingQualityTooLow
> Mapping quality too low.
MultipleAlignments
> Read aligns to more than one region.
OrphanRead
> Read's mate pair is unmapped.
OverlapRead
> Paired-end read covers the site of iterest redundantly.
QueryIndelAtSNVLocus    
> Indel in the read at the site of interest.
SNVLocusAtEndOfRead
> Site of interest is too close to the start or end of the read.
SubstitutionNearSNVLocus        
> Site of interest is too close to another substitution. 
TooManyEdits
> Read alignment requires too many edits.
TooManyEditsOtherThanSNV
> Read alignment requires too many edits other than at the site of interest.
TooManyQueryGaps
> Read is aligned using too many gaps. 
TooShort
> Read is too short.
