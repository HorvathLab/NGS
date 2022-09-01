
# varLoci

The SCReadCounts release contains a new tool `varLoci` for enumerating *potential* variant loci directly from the indexed BAM file. It outputs SNV loci and reference and alternate nucleotides suitable for SCReadCounts SNVs input file. 

## Usage

```
% varLoci <bam_file>.bam <region> <min_var_read_count> > snv_loci.txt
```

The ```<bam_file>.bam``` must be indexed and ```<bam_file>.bam.bai``` should be in the same directory. 

```<region``` is the samtools format region specifier. Use ```<chrom>``` or ```<chrom>:<start>-<end>``` to specify a chromosome or a chromosomal region, respectively. Use ```-``` to indicate no region constraint. 

```<min_var_read_count>``` is the minimum number of variant reads required for a locus to be considered a putative variant loci.

Note that this is an agressive enumeration of *potential* variant sites for *much* more careful analysis by SCReadCounts or other tools. It will, however, provide a much more limited set of loci than an exhaustive enumeration of loci in a genomic region. It will also provide the opportunity for *de novo* discovery of loci that have not otherwise been annotated elsewhere.
