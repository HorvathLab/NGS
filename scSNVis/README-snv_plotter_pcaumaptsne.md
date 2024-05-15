# snv\_plotter_tsne.r (Updated: 05/14/2024)
# snv\_plotter_pca.r (Updated: 05/14/2024)
# snv\_plotter_umap.r (Updated: 05/14/2024)
## Introduction
snv\_plotter_tsne\/pca\/umap.r scripts plot information about SNVs provided by the user onto a
dimensionally reduced representation of the data (either tsne, pca, or umap,
depending on the script selected). Currently, the script generates a set of 
plots that represent basic statistics and properties of the SNVs identified
in the dataset and utilize existing packages, such as Seurat, Slingshot, and scType.
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
 
## Input
The script accepts three inputs:
- STAR solo output directory (-m) that contains features.tsv.gz, barcodes.tsv.gz, and matrix.mtx.gz OR a Seurat object (-r) that was saved as an RDS
- SCReadCounts output file (-t) saved as a .tsv

See sample files for reference.

## Examples
ml R
Rscript 240508_multivar_plotter_tsne_blRed.r -t SNV_file.tsv -m SAMNXX_wasp_Solo.out/Gene/filtered/
Rscript 240508_multivar_plotter_tsne_blRed.r -t SNV_file.tsv -m SAMNXX_wasp_Solo.out/Gene/filtered/ --th-vars=1 --th-reads=10 --tissue-type=Liver -c -d -e

