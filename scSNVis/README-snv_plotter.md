# snv\_plotter.r (Updated: 11/28/2022)
## Introduction
snv\_plotter.r plots information about SNVs provided by the user onto a
dimensionally reduced representation of the data. Currently, the script
generates a set of plots that represents basic statistics and properties of the
SNV's identified in the dataset, which includes:
- Histogram of different SNVs for each cell
  - This figure shows how many different SNVs were identified in each cell.
    This information is represented in the file names and in the script as 'n'.
- Histogram of total reads of SNV for each cell
  - This figure shows how many total reads of SNVs were identified in each
    cell. The number will be inherently larger than the total different SNVs
    for each cell. The number of reads is represented as 'reads' in the file
    names and in the script. 
- UMAP representation of different SNVs for each cell.
  - This figure shows the number of different SNV's found in each cell, where
    the cells with the highest count is labeled with dark red while the cells
    with the lowest count is labeled with yellow.
- UMAP representation of total reads of SNV for each cell.
  - This figure is similar to the previous figure except the figure shows the
    number of reads for SNVs found in each cell.
- UMAP representation with variable allelic frequency (VAF)
  - For each SNV, this figure shows the VAF for each SNV. VAF is the percentage
    of observed SNV reads from the total number of reads.

## Input
The script accepts three inputs:
- .RDS file (-r)  - This is a Seurat object that was written into an RDS
  format. The Seurat object can either be an integrated dataset from multiple
  runs or a single run. Seurat's DefaultAssay function should be used before it
  is loaded into the script to designate which data should be used for
  processing.
- .txt file containing paths to SCReadCounts files (-t) - This is a .txt file
  that contains the path to all of the SCReadCounts output files. This should
  be in the same order in which the datasets were loaded into Seurat for any
  integration, as this order is assumed to be the order which the data appears
  in the Seurat object. A sample name should also be supplied in the second
  column.
- .txt file containing list of SNVs (chromosome, position) (-l, optional) -
  This is a .txt file that contains a list of individual SNV's for which
  individual plots will be generated. These individual plots are UMAP
  representations with data points showing which cells contain the SNV's.

The files passed to the script are all tab-separated files except the .RDS
file. The text file containing a list of SCReadCounts should be constructed
such that the first column contains the path to the SCReadCounts output file
and the second column contains a name for the dataset.

The text file containg the SNV positions for which individual plots should be
generated should also be a tab-separated file. The first column should contain
the chromosome, and the second column should contain the position.

See sample files (snv.txt and short\_list.txt) for reference.

## Example
`RScript snv_plotter.r -r comb.rds -t snv.txt -l short_list.txt`
