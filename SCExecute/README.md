
SCExecute generates cell-barcode specific BAM files from aligned,
aggregate single-cell sequencing data, executing a user-provided command
on each barcode-stratified BAM file. Unlike other tools, SCExecute
generates cell-barcode specific BAM files in batches to avoid file-system
and memory constraints, and manages the execution of the user-provided
commands on multiple processes to improve throughput. Cell-barcodes can
be [extracted](Barcodes.md) from read-names or BAM-file tags populated by a variety
of tools, included STARsolo and UMI-tools, and can be restricted to
barcodes of interest. Implemented in Python3 using the PySAM package and
[distributed](https://github.com/HorvathLab/NGS/releases) for Linux, MacOS, and Python environments, SCExecute builds
on other NGS tools from the Horvath lab, including [SCReadCounts](../SCReadCounts#readme).

SCExecute is [available](https://github.com/HorvathLab/NGS/releases) as a self-contained binary package for 64-bit
Linux systems, as Python source, and MacOS (Darwin). The self-contained
binary package is appropriate for most Linux and MacOS users. The pythonic
version requires pysam, numpy and scipy along with other packages (See
the [install](docs/Installation.md) instructions for more details).

SCExecute provides explicit configuration for alignments barcoded
through STARsolo and UMItools. Additional cellular barcode
extraction logic can be configured software, based
on BAM file tags or RNA sequence name and delimited tokens or regular
expressions (see [Cell Barcode](Barcode.md) documentation). 

**Setup:**
* [Download](https://github.com/HorvathLab/NGS/releases/tag/SCExecute-1.3.1)
* [Install](docs/Installation.md)

**Usage:**
* [SCExecute](docs/Usage.md)

**Other:**
* [Input Files](docs/InputFiles.md)
* [Cell Barcodes](docs/Barcodes.md)
* [Command/Template Substitution](docs/CommandSubst.md)

**[Examples](docs/Examples.md)**

