# readCounts Usage

## Synopsis

### Graphical User Interface:

    readCounts.py

### Command-line:

    readCounts.py [options]

## Description

readCounts a computational framework for assessing the read counts bearing particular nucleotides at genomic positions of interest, following a statistical test to recognize the allelic read-count patterns that show little deviation from expected behavior.

## Graphical User Interface

![Options](readCounts2.jpg)

Click the help icon (question mark) at the top right of the GUI and
then an input field for help. Multiple files can be selected in the
file-chooser using Ctrl-Click or Shift-Click. Fields can be reset to
their default values using the Reset button. Click OK to execute
readCounts.

Additional GUI option tabs are documented below.

## Options

SNVs, -s SNVS, --snvs=SNVS

> Single-nucleotide-polymophisms (SNVs). Tabular and VCF format SNVs
> are supported. Multiple files are specified inside quotes, separated
> by spaces, and by using file globbing. See [Input
> Files](InputFiles.md) for more information. Required.

Read Alignment Files, -r ALIGNMENTS, --readalignments=ALIGNMENTS

> Read alignments files in indexed BAM format, with extension
> `.bam`. BAM index with extension `.bam.bai` must be located in the
> same directory. Multiple files are specified inside quotes,
> separated by spaces, and by using file globbing. See [Input
> Files](InputFiles.md) for more information. Required.

Output Folder, -o OUTPUT, --output=OUTPUT

> Output file. Will be created if necessary. See [Output Files](OutputFiles.md) for more information on output files. Required. 

--version

>Show program's version number and exit. 

-h, --help

>Show program help and exit.

### Advanced

![Advanced](readCounts3.jpg)

Min. Reads, -m MINREADS, --minreads=MINREADS

> Minimum number of good reads at each SNV locus per alignment file. Default=10.   

Max. Reads, -m MAXREADS, --maxreads=MAXREADS

> Scale read counts at high-coverage loci to ensure at
                        most this many good reads at SNV locus per alignment
                        file. Values greater than 1 indicate absolute read
                        counts, otherwise the value indicates the coverage
                        distribution percentile. Default=No maximum.


All Fields, -F, --full

> Output extra diagnostic read count fields.Default=False.


Filter Alignments, -f, --alignmentfilter

> (Turn off) alignment filtering by length, edits, etc.

Unique Reads, -U, --uniquereads   

> Consider only distinct reads.  

Threads/BAM, -t TPB, --threadsperbam=TPB                   

> Worker threads per alignment file. Indicate no threading with 0. Default=1.

Quiet, -q, --quiet

> Do not show readCounts progress.

