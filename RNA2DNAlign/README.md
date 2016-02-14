
RNA2DNAlign evaluates evidence for asymmetric allele distribution in
next-gen sequencing reads of DNA and RNA samples from the same
individual. RNA2DNAlign requires, as input: genome aligned reads and SNV
loci to analyze. Reads from each analysis type and sample must be
aligned to the same version of the human genome reference. SNVs may be
derived from the reads directly, using, for example, Samtools, or they
may be derived from independent sources, such as lists of known
annotated variants. Variant positions must correspond to genomic
coordinates of the reference genome used for the alignment.

RNA2DNAlign is available as a self-contained binary package for 64-bit
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
