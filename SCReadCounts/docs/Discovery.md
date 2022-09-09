
# varLoci

The SCReadCounts release contains a new tool `varLoci` for enumerating *potential* variant loci directly from the indexed BAM file. It outputs SNV loci and reference and alternate nucleotides suitable for SCReadCounts SNVs input file. 

## Usage

```
% varLoci <bam_file>.bam <min_var_read_count> [ <region> ] > snv_loci.txt
```

The ```<bam_file>.bam``` must be sorted (by alignment start position) and indexed. ```<bam_file>.bam.bai``` should be placed in the same directory as the BAM file. 

```<min_var_read_count>``` is the minimum number of reads with the variant allele required for a locus to be considered a putative variant loci.

```<region>``` (optional) is the samtools format region specifier. Use ```<chrom>``` or ```<chrom>:<start>-<end>``` to specify a chromosome or a chromosomal region, respectively. Default: no region constraint. 

Notes:
* **varLoci** is an agressive enumeration of *potential* variant sites for *much* more careful analysis by SCReadCounts or other tools. It will, however, provide a much more limited set of loci than an exhaustive enumeration of loci in a genomic region. It will also provide the opportunity for *de novo* discovery of loci that have not otherwise been annotated elsewhere. Finally, it will avoid the unnecessary examination of annotated loci that are not supported by the data available in the BAM file alignments. 
* **varLoci** requires the BAM file have the ```MD``` tag for all alignments in order to determine the reference nucleotide at each locus. Some aligners do not output this tag by default, so look for options that ensure the ```MD``` (and ```NM```) tags are output in the BAM file. Alternatively, the following ```samtools``` command will add the ```MD``` (and ```NM```) tags:
    ```
    % samtools calmd -b aln.bam ref.fasta > aln_md.bam
    % samtools index aln_md.bam
    ```
