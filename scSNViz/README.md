## Introduction
The scSNViz scripts plot information about SNVs provided by the user onto a
dimensionally reduced representation of the data (either tsne, pca, or umap,
depending on the script selected). Currently, the script generates a set of 
plots that represent basic statistics and properties of the SNVs identified
in the dataset and utilize existing packages, such as Seurat, Slingshot, and scType.
 
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
- [scSNViz_tsne.r](https://raw.githubusercontent.com/HorvathLab/NGS/master/scSNViz/scSNViz_tsne.r) (Updated: 05/28/2024)
- [scSNViz_pca.r](https://raw.githubusercontent.com/HorvathLab/NGS/master/scSNViz/scSNViz_pca.r) (Updated: 05/28/2024)
- [scSNViz_umap.r](https://raw.githubusercontent.com/HorvathLab/NGS/master/scSNViz/scSNViz_umap.r) (Updated: 05/28/2024)

## Examples

Rscript scSNViz_tsne.r -t sample_SNVs.txt -m SAMNXX_wasp_Solo.out/Gene/filtered/

Rscript scSNViz_umap.r -t sample_SNVs.txt -m SAMNXX_wasp_Solo.out/Gene/filtered/ --th-vars=1 --th-reads=10 --tissue-type=Liver -c -d -e

Rscript scSNViz_umap.r -t sample_SNVs.txt -r sample_seurat.rds --th-vars=1 --th-reads=10 --tissue-type=Immunesystem -c -d -e --color-scale=YlOrRd





