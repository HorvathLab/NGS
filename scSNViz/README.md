## Introduction
The scSNViz script provides a comprehensive visualization of single cell-specific expressed SNVs (sceSNVs) by mapping them onto a dimensionally reduced representation of the data (e.g., UMAP, t-SNE, or PCA, depending on the selected technique). It generates a series of plots that represent the basic statistics and properties of the SNVs in the dataset. scSNViz integrates existing packages, including Seurat, Slingshot, and scType.
 
## Input
The script requires two inputs:
- STAR solo output directory (-m) that contains features.tsv.gz, barcodes.tsv.gz, and matrix.mtx.gz OR a Seurat object (-r) that was saved as a .RDS
- SCReadCounts output file (-t) that is tab-delimited and either a .tsv or .txt
  
Additionally, lists of SNVs with cell-barcode information not processed through SCReadCounts can be submitted in similar format using the (-t) option.

See sample files for reference.

## Options
	-h, --help
		Show this help message and exit

	-r RDS-FILE, --rds-file=RDS-FILE
		RDS file containing Seurat object.

	-m COUNTSMATRIX-FILE, --countsmatrix-file=COUNTSMATRIX-FILE
		folder containing STARsolo output folder name that contains
                     the following files:
                     barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

	-t SNV-FILE, --snv-file=SNV-FILE
		scReadCounts file

	-w DIMENSIONALITY-REDUCTION, --dimensionality-reduction=DIMENSIONALITY-REDUCTION
		options include tSNE, PCA, UMAP. Default=UMAP.

	-x TH-VARS, --th-vars=TH-VARS
		Threshold for number of sceSNVs. Default=0 (display cells with N_SNVs > 0).

	-y TH-READS, --th-reads=TH-READS
		Threshold for number of variant reads (N_VAR). Default=0 (consider as sceSNV positions covered with N_VAR > 0).

	-c, --enable-title
		Enable title. Default=T

	-d, --disable-ind-plots
		Disable individual SNV plots. Default=F.

	-e, --disable-3d-axis
		Disable axes in 3D plots. Default=F.

	-g, --disable-slingshot
		Disable slingshot curves in 3D plots. Default=F.

	-i, --enable-sctype
		Enable sctype to run. Default=F.

	-j TISSUE-TYPE, --tissue-type=TISSUE-TYPE
		tissue type for scType; options include:
                     Immunesystem, Pancreas, Liver, Eye, Kidney, Brain,
                     Lung, Adrenal, Heart, Intestine, Muscle, Placenta,
                     Spleen, Stomach, Thymus

	-k COLOR-SCALE, --color-scale=COLOR-SCALE
		if you would like to change the default color settings with
                     these options, you may use Blues, Reds, YlOrRd, YlGnBu, plasma, RdBu

	-p, --enable-cell-border
		Enable cell border. Default=F

	-q, --enable-dynamic-cell-size
		Enable cell size to depend on number of reads. Default=F

## Output
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
  - Expressed Variant Allele Fraction per cell (VAF_RNA)
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
                         -c -d -e --color-scale=YlGnBu 
```





