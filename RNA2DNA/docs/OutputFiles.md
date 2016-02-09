# RNA2DNA Output Files

RNA2DNA output files are created in the directory specified. The
folder will be created if necessary. Existing files will be
overwritten.

## Summary

The `summary_result.txt` file summarizes the count of each type of event observed.

See example output files in the `RNA2DNA/data` directory.

## Event Files

Each execution of RNA2DNA will create (up to) eight tab-separated value event files
representing the following events: RNA editing (`Events_RNAed.tsv`), tumor-specific
RNA editing (`Events_T-RNAed.tsv`), variant-specific expression (`Events_VSE.tsv`) or loss
(`Events_VSL.tsv`), tumor-specific variant expression (`Events_T-VSE.tsv`) or loss (`Events_T-VSL.tsv`),
somatic mutagenesis (`Events_SOM.tsv`), and loss of heterozygosity (`Events_LOH.tsv`).

See example output files in the `RNA2DNA/data` directory.

### Event File Fields

AlignedReads
> Name of the aligned reads file for the following read counts.

CHROM  
> Chromosome identifier

POS  
> Chromosome position of the variant

REF  
> Reference allele nucleotide

ALT  
> Variant allele nucleotide

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

## Read Counts

A tab-separated values file consisting of the computed read-counts is
also provided (`readCounts.tsv`). This file contains the read counts
for each SNV locus in each BAM file and computes the various
statistical tests described above, in "Event File Fields". The read
counts file can be used to investigate the computed values for
expected events that didn't pass filtering, significance, or scoring
thresholds.

## See Also

[RNA2DNA Home](..), [Usage](Usage.md), [Input Files](InputFiles.md), [Annotation Files](AnnotationFiles.md)

