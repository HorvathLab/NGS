## Introduction
The scSNViz script provides a comprehensive visualization of single cell-specific expressed SNVs (sceSNVs) by mapping them onto a dimensionally reduced representation of the data (e.g., UMAP, t-SNE, or PCA, depending on the selected technique). It generates a series of plots that represent the basic statistics and properties of the SNVs in the dataset. scSNViz integrates existing packages, including Seurat, Slingshot, and scType.
 
## Input
The script requires two inputs:
- STAR solo output directory (-m) that contains features.tsv.gz, barcodes.tsv.gz, and matrix.mtx.gz OR a Seurat object (-r) that was saved as a .RDS
- SCReadCounts output file (-t) that is tab-delimited and either a .tsv or .txt
  
Additionally, lists of SNVs with cell-barcode information not processed through SCReadCounts can be submitted in similar format using the (-t) option.

See sample files for reference.

## Default Output
scSNViz generates outputs for the set of the sceSNVs and for each individual sceSNV, as follows:

For the set of the sceSNVs, the produced figures include:
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
- a text file of summary statistics per cell
  
For the individual sceSNVs in the set, scSNViz generates gradient color representation plots of the following:
  - VAF_RNA per cell
  - Number of Reference Reads per cell (N_REF)
  - Number of Variant Reads per cell (N_VAR)

## Installation

Download the R file: 
- [scSNViz.r](https://raw.githubusercontent.com/HorvathLab/NGS/master/scSNViz/scSNViz.r) (Updated: 06/14/2024)

The following CRAN packages are required:
- optparse, stringr, openxlsx, HGNChelper, Seurat, ggplot2, dplyr, plotly, htmlwidgets.

The following Bioconductor packages are required:
- glmGamPoi, slingshot.

Note too that the matrixStats package (a dependancy of Seurat) needs to be downgraded to version 1.1.0:
- `> remotes::install_version("matrixStats", version="1.1.0")`

## Examples
```
% Rscript scSNViz.r -t sample_SNVs.txt -m SAMNXX_wasp_Solo.out/Gene/filtered/
```
```
% Rscript scSNViz.r -t sample_SNVs.txt -m SAMNXX_wasp_Solo.out/Gene/filtered/
                         --dimensionality-reduction=umap \
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





