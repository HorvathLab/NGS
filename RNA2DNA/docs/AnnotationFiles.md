# RNA2DNA Annotation Files

**All annotation files, SNV loci, and read alignments must indicate
genomic position with respect to the same specific release of a common
reference genome.**

## RefSeq Exon Coordinates

Use of exon coordinates to filter the SNVs is strongly recommended.

1. RefSeq exon coordinates are downloaded from the UCSC genome browser
and provided in the RNA2DNAlign/data directory in the files:
`UCSC_Human_hg19_RefSeq_CDS_exon_coordinates.txt` . These annotation
files can be used as provided.

2. RefSeq exon coordinates can be recreated for another organism or assembly as follows:

```
    cd data
    ./dlexons.sh hg19 > UCSC_Human_hg19_RefSeq_CDS_exon_coordinates.txt
```

   Exon coordinates should be tab-separated and sorted by chromosome number (1,2,3,...,X,Y), start position, end position, in increasing order. 

## COSMIC Mutants

1. Register with COSMIC here:

```
   https://cancer.sanger.ac.uk/cosmic/register
```

2. Download the COSMIC mutants:

```
   sftp://sftp-cancer.sanger.ac.uk//files/grch38/cosmic/v74/CosmicMutantExport.tsv.gz
```

3. This file is used in its published format.

## DARNED Loci

1. Download the DARNED loci from here (NCBI37/hg19):

```
   http://darned.ucc.ie/static/downloads/hg19.txt
```

2. This file is used in its published format.

