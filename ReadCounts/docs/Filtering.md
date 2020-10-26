# ReadCounts Read Filtering

The aligned reads covering a specific locus may, for a variety of
reasons, not be a reliable observation of the reference or variant
alleles at that site. The readCounts suite provides a flexible filtering 
system that can be tailored as needed to ensure reliable allele counts.

The available read-filtering strategies are defined in the file
"filter.ini" in the readCounts distribution. New or modified
read-filtering strategies can be created, using the same format in
a filter.ini in the current working directory. Named filtering
strategies in the current working directory override those with the
same name in the readCounts distribution.

## Read-Filtering Strategies

### Basic
> Filter out alignments flagged as SECONDARY, DUPLICATE, UNMAPPED, or QCFAIL, and those with a gap or indel at SNV locus. This is the minimal recommended filtering strategy.

### ReadCounts:Defaults
> Original ReadCounts alignment filtering strategy with default parameters.

### ReadCounts:Parameterized
> ReadCounts alignment filtering strategy with explicit parameter values.

### Explicit:Defaults
> Filter in which every filter rule is applied with explicit default parameters. Unmodified, this filtering strategy is equivalent to the Basic filter.

### MPileup
> Filter reads in a manner similar to samtools/vcftools/bcftools mpileup command-line tool. Filters implement the --ff, -Q, -A, and -x mpileup options.

## Read Filters

### BasicFilter
> Parameters: skip_duplicate=True skip_secondary=True skip_qcfail=True skip_unmapped=True

> Optionally filter out alignments flagged as SECONDARY, DUPLICATE, UNMAPPED, or QCFAIL, and those with a gap or indel at SNV locus.  

### SNVPileupReadFilter
> Parameters: minpad=3 minsubstdist=3 maxedits=1 maxsegments=1 minlength=45 maxhits=1 mapq=4

> Filter implementing the original ReadCounts filtering strategy. Deprecated.

### BaseQualityFilter
> Parameters: min_base_quality=None

> Minimum required base quality at the specific locus. Default: no minimum base quality. 

### MappingQualityFilter
> Parmeters: min_mapping_quality=None

> Minimum required mapping quality at the specific locus. Default: no minimum mapping quality. 

### ReadLengthFilter
> Parameters: min_length=None

> Minimum required length for the read. Default: No minimum length.

### EditsFilter
> Parameters: max_edits=None

> Maximum number of edits for the alignment. Includes the substitution observed in reads for the variant allele. References the NM tag. Default: No maximum number of edits. 

### HitsFilter
> Parameters: max_hits=None

> Maximum number of alignments placements for the read. References the NH tag. Default: No maximum number of hits. 

### EditPositionFilter
> Parameters: min_edge_dist=None min_subst_dist=None max_other_edits=None

> Filter reads based on a) the proximity of the site of interest to the beginning or end of the read, b) the proximity of the nearest edit to the site of interest, or c) total number of edits other than at the site of interest. Default: Keep all reads. 

### SegementsFilter
> Parameters: max_segments=None

> Maximum number of segments (separated by gaps) in the read alignment. Default: No maximum number of segments. 

### OverlapFilter
> Pameters: remove=False

> When a site of interest is covered by both reads of a read-pair, just one is retained. Default: Retain all reads. 

### OrphanFilter
> Parameters: remove=False

> When a paired-end read's mate is unmapped, do not count the read. Default: Ratain all reads. 

### UniqueReads
> Parameters: remove_dups=False

> Ensure exactly duplicate reads are counted exactly once. Default: Do not remove duplicate reads. 

## Examples

### Basic

```
[Basic]
Description:    (Optionally) Filter out alignments flagged as SECONDARY, DUPLICATE,
		UNMAPPED, or QCFAIL, and those with a gap or indel at
		SNV locus. This is the minimal recommended filtering strategy.
BasicFilter: skip_duplicate=True skip_secondary=True skip_qcfail=True skip_unmapped=True
```

### MPileup

```

[MPileup]
Description: Filter reads in a manner similar to samtools/vcftools/bcftools
             mpileup command-line tool. Filters implement the --ff, -Q, -A,
	     and -x mpileup options.
BasicFilter:          skip_duplicate=True
                      skip_secondary=True
                      skip_qcfail=True
                      skip_unmapped=True
BaseQualityFilter:    min_base_quality=13
OverlapFilter:        remove=True
OrphanFilter:         remove=True
```
