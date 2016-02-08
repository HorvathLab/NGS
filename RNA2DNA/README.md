
RNA2DNA evaluates evidence of allelic imbalance and asymmetry in next-gen
sequencing reads of exomes and RNA from normal and tumor samples from
the same individual. RNA2DNA requires, as input: genome aligned reads
and SNV loci to analyze. Reads from each analysis type and sample must
be aligned to the same genomic coordinates. SNVs may be derived from
the reads directly, using, for example, TopHat2 and samtools, or they
may be derived from independent sources.

RNA2DNA is available as a self-contained binary package for 64-bit Linux
systems and as Python source. The
[pysam](https://github.com/pysam-developers/pysam) module (version >= 0.8.1) must be
installed to use in Python source form. We recommend the binary package
for x86_64 Linux systems.

[Download](https://github.com/HorvathLab/NGS/releases/) and [Install](docs/Installation.md).

Usage:
* [RNA2DNA](docs/Usage.md)

File Formats:
* [Input Files](docs/InputFiles.md)
* [Output Files](docs/OutputFiles.md)
* [Annotation Files](docs/AnnotationFiles.md)

Examples:
* [Minimal Example](docs/ExampleAnalysis.md)
* [TCGA Analysis](docs/TCGAExampleAnalysis.md)
