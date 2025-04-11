SCExecute generates cell-barcode specific BAM files from aligned,
aggregate single-cell sequencing data, executing a user-provided command
on each barcode-stratified BAM file. Unlike other tools, SCExecute
generates cell-barcode specific BAM files in batches to avoid file-system
and memory constraints, and manages the execution of the user-provided
commands on multiple processes to improve throughput. Cell-barcodes can
be [extracted](docs/Barcodes.md) from read-names or BAM-file tags populated by a variety
of tools, included STARsolo and UMI-tools, and can be restricted to
barcodes of interest. Implemented in Python3 using the PySAM package and
[distributed][Current version] for Linux and Python environments, SCExecute builds
on other NGS tools from the Horvath lab, including [SCReadCounts](../SCReadCounts#readme).

SCExecute is [available][Current version] as a self-contained binary package for 64-bit
Linux systems and as Python source. The self-contained
binary package is appropriate for most Linux users. The pythonic
version requires pysam, numpy and scipy along with other packages (See
the [install](docs/Installation.md) instructions for more details).

SCExecute provides explicit configuration for alignments barcoded
through STARsolo and UMItools. Additional cellular barcode
extraction logic can be configured software, based
on BAM file tags or RNA sequence name and delimited tokens or regular
expressions (see [Cell Barcode](docs/Barcode.md) documentation). 

**Setup:**
* [Download][Current version]
* [Install](docs/Installation.md)

**Usage:**
* [SCExecute](docs/Usage.md)

**Other:**
* [Input Files](docs/InputFiles.md)
* [Cell Barcodes](docs/Barcodes.md)
* [Command/Template Substitution](docs/CommandSubst.md)

**[Examples](docs/Examples.md)**

**Please cite:**
* Edwards, N., Dillard, C., Prashant, N.M., Liu, H., Yang, M., Ulianova, E., and Horvath, A. [SCExecute: custom cell barcode-stratified analyses of scRNA-seq data](https://doi.org/10.1093/bioinformatics/btac768). Bioinformatics (2022).

[Current version]: https://github.com/HorvathLab/NGS/releases/tag/SCExecute-1.3.3
