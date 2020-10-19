
ScReadCounts is a computational tool for a cell-level assessment of the read counts bearing a particular nucleotide at genomic positions of interest from single cell RNA sequencing (scRNA-seq) data. 

ScReadCounts is available as a self-contained binary package for 64-bit
Linux systems and as Python source. The pysam package, plus a variety
of common third-party python packages including numpy and scipy must
be installed to use in Python source form. See the install
instructions for more details. The self-contained binary package is
appropriate for most Linux users.

Currently, scReadCounts has two programs. The program readCounts requires two input files: a pooled single cell alignment and a list of genomic positions of interest. readCounts utilizes the barcode information from the pooled single cell alignments and outputs the variant and reference read counts (n_var and n_ref, respectively), for each barcode (cell), in a comma separated text file. This file is then used as an input for the second program - readCountsMatrix - which, upon providing an output directory, generates two outputs: (1) a cell-position matrix with n_var and n_ref estimates, and (2) a cell-position matrix with the expressed variant allele fraction (VAF_RNA = n_var / (n_var + n_ref)) estimated at a user-defined threshold of minimum required sequencing reads (minR). readCountsMatrix is very fast and can be re-run multiple times at various minR thersholds.


#### To get help on ScReadCount ####
Go to ScReadCounts directory (../NGS/SCReadCounts):
```
$ python3 src/scReadCounts.py –help
```

Setup:
* [Download](https://github.com/HorvathLab/NGS/releases/tag/SCReadCounts-1.0.0)
* [Install](docs/Installation.md)

Usage:
* [SCReadCounts](docs/Usage.md)

File Formats:
* [Input Files](docs/InputFiles.md)
* [Output Files](docs/OutputFiles.md)

[Examples](docs/Examples.md)
