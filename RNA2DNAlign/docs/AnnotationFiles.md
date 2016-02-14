# RNA2DNAlign Annotation Files

**All annotation files, SNV loci, and read alignments must indicate
genomic position with respect to the same specific release of a common
reference genome.**

The following files and instructions are for hg19/grch37/NCBI37.

## RefSeq Exon Coordinates from UCSC Genome Browser

Use of exon coordinates to filter the SNVs is strongly recommended.

1. RefSeq exon coordinates are downloaded from the UCSC genome browser
and provided in the RNA2DNAlign/data directory in the file:
`UCSC_Human_hg19_RefSeq_CDS_exon_coordinates.txt`. The annotation
file can be used as provided.

2. RefSeq exon coordinates can be recreated for another organism or assembly as follows:

```
    cd data
    ./dlexons.sh hg19 > UCSC_Human_hg19_RefSeq_CDS_exon_coordinates.txt
```

   Exon coordinates should be tab-separated and sorted by chromosome number (1,2,3,...,X,Y), start position, end position, in increasing order. 

## [COSMIC](http://cancer.sanger.ac.uk/cosmic) Mutants

1. COSMIC mutants are downloaded from the COSMIC website and provided in the RNA2DNAlign/data directory in the file: `CosmicMutantExport_hg19.tsv.gz`.  The annotation
file can be used as provided.

2. COSMIC mutants can be downloaded for another organism or assembly as follows: 
    * Register with COSMIC
    ```
    https://cancer.sanger.ac.uk/cosmic/register
    ```
    * Download the COSMIC mutants:
    ```
    sftp "login"@sftp-cancer.sanger.ac.uk:/cosmic/grch37/cosmic/v75/CosmicMutantExport.tsv.gz
    ```
    * COSMIC mutant annotations can be used in its downloaded format. 

## [DARNED](http://darned.ucc.ie/) Loci

1. DARNED loci are downloaded form the DARNED website and provided in
the RNA2DNAlign/data directory in the file: `DARNED_hg19.txt`.  The
annotation file can be used as provided.

2. DARNED loci can be downloaded for another organism or assembly as follows: 
    * Download the DARNED loci:
    ```
    http://darned.ucc.ie/static/downloads/hg19.txt
    ```
    * DARNED loci can be used in its downloaded format.

[RNA2DNAlign Home](..), [Usage](Usage.md), [Input Files](InputFiles.md), [Output Files](OutputFiles.md)
