
SCReadCounts is a computational tool for a cell-level assessment of the read counts bearing a particular nucleotide at genomic positions of interest from single cell RNA sequencing (scRNA-seq) data. 

SCReadCounts is available as a self-contained binary package for 64-bit Linux systems, as Python source, and MacOS (Darwin). The self-contained binary package is appropriate for most Linux and MacOS users. The pythonic version requires pysam, numpy and scipy along with other packages (See the install instructions for more details). 

Currently, SCReadCounts has three programs. The program scReadCounts manages the sequential execution of programs readCounts and readCountsMatrix, collects the necessary arguments for successful execution, and avoids unnecessary execution of the expensive readCounts tool if possible. readCounts requires three input files: a pooled single cell alignment, a list of genomic positions of interest, and the barcodes file (barcodes.tsv) . Optionally, readCounts can be user-configured for read filtering and cell-barcode handling, including restriction to barcodes of interest (achieved through specifying barcodes file - barcodes.tsv). readCounts utilizes the barcode information from the pooled single cell alignments and outputs the variant and reference read counts (n<sub>var</sub> and n<sub>ref</sub>, respectively), for each barcode (cell), restricted to those present in the barcodes file, in a tab separated text file. This file is then used as an input for the second program - readCountsMatrix - which, upon providing an output prefix, generates two outputs: (1) a cell-position matrix with absolute n<sub>var</sub> and n<sub>ref</sub> counts, and (2) a cell-position matrix with the expressed VAF<sub>RNA</sub>. VAF<sub>RNA</sub> is estimated at a user-defined threshold of minimum required sequencing reads (minR); default minR = 5. readCountsMatrix is time-efficient and can be re-run multiple times at various minR thresholds. Together, these tools facilitate single-cell level assessment of read counts.

SCReadcounts provides explicit configuration for alignments barcoded through STARsolo, CellRanger and UMItools. Additional cellular barcode extraction logic can be configured for SCReadCounts software, based on BAM file tags or RNA sequence name and delimited tokens or regular expressions (see the “SCReadCounts Read Grouping” documentation).

SCReadCounts is a wrapper around readCounts and readCountsMatrix function to facilitate single-cell level assessment of read counts.

**Setup:**
* [Download](https://github.com/HorvathLab/NGS/releases/tag/SCReadCounts-1.1.10)
* [Install](docs/Installation.md)

**Usage:**
* [SCReadCounts](docs/Usage.md)
* [Discovery Mode](docs/Discovery.md)

**File Formats:**
* [Input Files](docs/InputFiles.md)
* [Output Files](docs/OutputFiles.md)

**[Examples](docs/Examples.md)**

NOTE: For any issues you are facing kindly provide the command line ouput that the software produces, and if possible, a small reproducible set of your inputs. It will help us to identify and resolve issues efficiently.
