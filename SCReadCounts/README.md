
SCReadCounts is a computational tool for a cell-level assessment of the read counts bearing a particular nucleotide at genomic positions of interest from single cell RNA sequencing (scRNA-seq) data. 

SCReadCounts is available as a self-contained binary package for 64-bit
Linux systems, as Python source, and MacOS (Darwin). The pythonic version requires pysam, numpy and scipy along with some other packages (See the install instructions for more details). The self-contained binary package is appropriate for most Linux users.

Currently, SCReadCounts has two programs. The program readCounts requires two input files: a pooled single cell alignment and a list of genomic positions of interest. readCounts utilizes the barcode information from the pooled single cell alignments and outputs the variant and reference read counts (n_var and n_ref, respectively), for each barcode (cell), in a tab separated text file. This file is then used as an input for the second program - readCountsMatrix - which, upon providing an output prefix, generates two outputs: (1) a cell-position matrix with n_var and n_ref estimates, and (2) a cell-position matrix with the expressed variant allele fraction (VAF_RNA = n_var / (n_var + n_ref)). VAF_RNA is estimated at a user-defined threshold of minimum required sequencing reads (minR); default minR = 5. readCountsMatrix is time-efficient and can be re-run multiple times at various minR thresholds.


Setup:
* [Download](https://github.com/HorvathLab/NGS/releases/tag/SCReadCounts-1.1.1)
* [Install](docs/Installation.md)

Usage:
* [SCReadCounts](docs/Usage.md)

File Formats:
* [Input Files](docs/InputFiles.md)
* [Output Files](docs/OutputFiles.md)

[Examples](docs/Examples.md)
