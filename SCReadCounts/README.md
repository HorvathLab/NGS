
ScReadCounts is a tool for cell-level estimation of reference and variant read counts from scRNA-seq data. ScReadCounts utilizes barcode information from pooled single cell alignments and, provided a list of SNV sites, estimates the variant and reference read counts, calculates variant allele fraction at user-defined treshold for minimum number of reads, and formats counts and derived values as cell-SNV matrices. The cell-SNV matrices can be used as inputs for a number of downstream analyses. 

ScReadCounts is available as a self-contained binary package for 64-bit
Linux systems and as Python source. The pysam package, plus a variety
of common third-party python packages including numpy and scipy must
be installed to use in Python source form. See the install
instructions for more details. The self-contained binary package is
appropriate for most Linux users.

Setup:
* [Download](https://github.com/HorvathLab/NGS/releases/)
* [Install](docs/Installation.md)

Usage:
* [RNA2DNAlign](docs/Usage.md)

File Formats:
* [Input Files](docs/InputFiles.md)
* [Output Files](docs/OutputFiles.md)
* [Annotation Files](docs/AnnotationFiles.md)

[Examples](docs/Examples.md)
