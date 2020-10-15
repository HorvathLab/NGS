
ScReadCounts is a tool for cell-level estimation of reference and variant read counts from scRNA-seq data. ScReadCounts utilizes barcode information from pooled single cell alignments and, provided a list of SNV sites, estimates the variant and reference read counts, calculates variant allele fraction at user-defined treshold for minimum number of reads, and formats counts and derived values as cell-SNV matrices. The cell-SNV matrices can be used as inputs for a number of downstream analyses. 

ScReadCounts is available as a self-contained binary package for 64-bit
Linux systems and as Python source. The pysam package, plus a variety
of common third-party python packages including numpy and scipy must
be installed to use in Python source form. See the install
instructions for more details. The self-contained binary package is
appropriate for most Linux users.

The SNVS (-s) and ALIGNMENTS (-r) options are used to provide the genome aligned reads and SNV loci to analyze. OUTPUT(-o) is required for the output matrix.

* -s SNVS, --snvs=SNVS 
  *	Single-Nucleotide-Variant files. Required.

* -r ALIGNMENTS, --readalignments=ALIGNMENTS
  * Read alignment files in indexed BAM format. Required.

* -o OUTPUT, --output=OUTPUT
  * Output file. Required.

The remaining options provide detailed settings to better narrow the results and provide desired output formats. They are optional.

* -f FILTER, --alignmentfilter=FILTER
  * Alignment filtering strategy. Default: Basic.

* -m MINREADS, --minreads=MINREADS
  * Minimum number of good reads at SNV locus per alignment file. Default=3.

* -M MAXREADS, --maxreads=MAXREADS
  * Scale read counts at high-coverage loci to ensure at most this many good reads at SNV locus per alignment file. Values greater than 1 indicate absolute read counts, otherwise the value indicates the coverage distribution percentile. Default=No maximum.

* -G READGROUP, --readgroup=READGROUP
  * Additional read grouping based on read name/identifier strings or BAM-file RG. Default: UMITools cell barcodes ("UMITools").

* -t TPB, --threadsperbam=TPB
  * Worker threads per alignment file. Indicate no threading with 0. Default=0.

* -F, --full
  * Output extra diagnostic read count fields. Default=False.

* -U, --uniquereads
  * Consider only distinct reads.

* -q, --quiet
  * Quiet.

* -d, --debug
  * Debug.
  
## Download & Installation ##
#### Download directly from the git repository: ####
```
$ git clone https://github.com/HorvathLab/NGS.git
```
#### Download from Nathan Edwards Lab: ####
[Nathan Edwards Lab](http://edwardslab.bmcb.georgetown.edu/software/downloads/HorvathLab/)



## Tutorial ##
#### To get the read counts of RNA seq “singlecell_chr17.bam” regarding the SNV loci in file “singlecell_222_5_chr17.txt” in the single cell setting and save the output in file “singlecell-output.tsv”. ####

### Python ###
Go to SCReadCounts directory (../NGS/SCReadCounts):

```
$ python3 src/scReadCounts.py -s data/singlecell_222_5_chr17.txt -r data/singlecell_chr17.bam -o data/singlecell-output.tsv
```
To access the VAF matrix.
```
$ less data/singlecell-output.vaf.matrix.tsv
```

## Help ##
For full documentation see ()

#### To get help on ScReadCount ####
Go to ScReadCounts directory (../NGS/SCReadCounts):
```
$ python3 src/scReadCounts.py –help
```
