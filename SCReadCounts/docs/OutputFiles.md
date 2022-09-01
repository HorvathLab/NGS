# SCReadCount Output 

SCReadCounts output files are created based on the output filename provided. 

## SCReadCounts output file

A tab-separated values file consisting of the computed counts for each SNV locus and cell-bardcode. 

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

SNVCountForward
> Number of forward oriented variant reads in the paired end alignment.

SNVCountReverse
> Number of reverse oriented variant reads in the paired end alignment.

RefCountForward
> Number of forward oriented reference reads in the paired end alignment.

RefCountReverse
> Number of reverse oriented reference reads in the paired end alignment.

SNVCount
> Total number of variant reads.

RefCount
> Total number of reference reads.

GoodReads
> Total number of good reads.

%BadRead
> Percentage of bad reads.

VAF
> Variant allele fraction

## SCReadCounts Counts matrix

File with extention `*.cnt.matrix.<extn>` for output file with extension `*.<extn>`. Contains read counts for the reference and alternative allele at each locus and in each cell. Rows represent loci, columns represent cell barcodes. Values are `refcnt;altcnt`. 

## SCReadCounts VAF matrix

File with extension `*.vaf.matrix.<extn>` for output file with extension `*.<extn>`. Contains variant allele frequency (VAF) for loci with at least the minimum required number of good reads (Min. Reads or `-m` option) at each locus and in each cell. Rows represent loci, columns represent cell barcodes. Loci and cells without sufficient reads are indicated with `NA`.


