
readCounts a computational framework for assessing the read counts bearing
particular nucleotides at genomic positions of interest, following a statistical
test to recognize the allelic read-count patterns that show little deviation
from expected behavior.

readCounts requires, as input: genome aligned reads and SNV loci to analyze. Reads from each analysis type and sample must be aligned to the same version of the human genome reference. SNVs may be derived from the reads directly, using, for example, Samtools, or they may be derived from independent sources, such as lists of known annotated variants. Variant positions must correspond to genomic coordinates of the reference genome used for the alignment.

Setup:
* [Download](https://github.com/HorvathLab/NGS/releases/)
* [Install](docs/Installation.md)

Usage:
* [ReadsCount](docs/Usage.md)

File Formats:
* [Input Files](docs/InputFiles.md)
* [Output Files](docs/OutputFiles.md)

[Examples](docs/Examples.md)
