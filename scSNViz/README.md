## Introduction
The scSNViz scripts plot information about SNVs provided by the user onto a
dimensionally reduced representation of the data (either tsne, pca, or umap,
depending on the dimensionality reduction technique selected selected). 
Currently, the script generates a set of plots that represent basic statistics 
and properties of the SNVs identified in the dataset and utilize existing packages,
such as Seurat, Slingshot, and scType.
 
## Input
The script accepts three inputs:
- STAR solo output directory (-m) that contains features.tsv.gz, barcodes.tsv.gz, and matrix.mtx.gz OR a Seurat object (-r) that was saved as a .RDS
- SCReadCounts output file (-t) that is tab-delimited and either a .tsv or .txt

See sample files for reference.

## Default Output
The produced figures include:
- Histogram of mean VAF per SNV per cell
- Histogram of the number of SNVs per cell
- Histogram of the number of Variant Reads per cell
- Histogram of the Total VAF per cell (VARreads/(VARreads+REFreads) per cell)
  
- 3D UMAP representation of mean VAF for each cell
- 3D UMAP representation of median VAF for each cell
- 3D UMAP representation of number of number of SNVs for each cell
- 3D UMAP representation of number of Variant Reads for each cell
- 3D UMAP representation of number of Reference Reads for each cell
- 3D UMAP representation of Total VAF per cell (VARreads/(VARreads+REFreads) per cell)

- Individual SNV plots
  - VAF per cell
  - Number of Reference Reads per cell
  - Number of Variant Reads per cell

- a text file of summary statistics per cell

## Installation

Download the file according to the desired dimensionality reduction technique: 
- [scSNViz.r](https://raw.githubusercontent.com/HorvathLab/NGS/master/scSNViz/scSNViz.r) (Updated: 05/28/2024)

The following CRAN R packages are required:
- optparse, stringr, openxlsx, HGNChelper, Seurat, ggplot2, dplyr, plotly, htmlwidgets.

The following Bioconductor R packages are required:
- slingshot.

Note too that the matrixStats package (a dependancy of Seurat) needs to be downgraded to version 1.1.0:
- `> remotes::install_version("matrixStats", version="1.1.0")`

## Examples
```
% Rscript scSNViz.r -t sample_SNVs.txt -m SAMNXX_wasp_Solo.out/Gene/filtered/
```
```
% Rscript scSNViz.r -t sample_SNVs.txt -m SAMNXX_wasp_Solo.out/Gene/filtered/ --dimensionality-reduction=umap \
                         --th-vars=1 --th-reads=10 \
                         -i --tissue-type=Liver -c -d -e 
```
```
% Rscript scSNViz.r -t sample_SNVs.txt -r sample_Seurat_object.rds \
                         --dimensionality-reduction=tsne \
                         --th-vars=1 --th-reads=10 \
                         -i --tissue-type=Immunesystem \
                         -c -d -e --color-scale=YlOrRd 
```





