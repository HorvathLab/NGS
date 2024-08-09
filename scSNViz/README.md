## Introduction
The scSNViz script provides a comprehensive visualization of single cell-specific expressed SNVs (sceSNVs) by mapping them onto a dimensionally reduced representation of the data (e.g., UMAP, t-SNE, or PCA, depending on the selected technique). It generates a series of plots that represent the basic statistics and properties of the SNVs in the dataset. scSNViz integrates existing packages, including Seurat, Slingshot, scType, and CopyKat.
 
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

	-c, --disable-title
		Disable title for individual SNV plots. Default=F

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

	-b, --enable-cell-border
		Enable cell border. Default=F

	-q, --enable-dynamic-cell-size
		Enable cell size to depend on number of reads. Default=F

  	-u, --enable-copykat
   		Enable CopyKat for displaying CNVs. Default=F

     	-s, --save-each-plot
      		Save plots as separate HTML files. Default=F

## Output

### Directory structure

scSNViz generates outputs for the set of the sceSNVs and for each individual sceSNV, as follows:

```
sample_SNVs_dimensionality_reduction_xry/
│
├── Exploratory_Combined_Plots.html
├── sample_SNVs-summary.txt
└── Figures_Individual_Plots_HTML/   (optional)
    ├── Cell_types_scType.html       (optional)
    ├── CNVs_CopyKat.html            (optional)
    ├── Median_VAF_RNA.html
    ├── Mean_VAF_RNA.html
    ├── N_REFreads.html
    ├── Total_VAF_RNA.html
    ├── N_sceSNVs.html
    ├── N_VARreads.html
    ├── Histogram_N_SNV.png
    ├── Histogram_N_VARreadsCounts.png
    ├── Histogram_MeanSNVsVAF.png
    ├── Histogram_TotalVAF.png
    └── Individual_sceSNVs/
        ├── VARreads/
        │   └── 3D VARreads plot HTML files for each sceSNV
        ├── REFreads/
        │   └── 3D REFreads plot HTML files for each sceSNV
        └── VAF/
            └── 3D VAF plot HTML files for each sceSNV
```

**Note:** The `Figures_Individual_Plots_HTML` directory and its contents are optional and will only be generated if the user requests it using the **-s** option. <br>
In `sample_SNVs_dimensionality_reduction_xry`, *dimensionality_reduction* indicates the choice of dimensionality reduction selected. *x* and *y* indicate the threshold for sceSNVs and the variant reads, respectively.

## Description

For the set of the sceSNVs, the separately produced figures show the following:

&nbsp;&nbsp;&nbsp;&nbsp;**MeanSNVsVAF**: Histogram of mean VAF per SNV per cell<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**N_SNV**: Histogram of the number of SNVs per cell<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**N_VARreadsCounts**: Histogram of the number of Variant Reads per cell<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**TotalVAF**: Histogram of the Total VAF per cell (VARreads/(VARreads + REFreads) per cell)<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**Cell_types_scType**: 3D UMAP/t-SNE/PCA representation of the cell types identified by scType<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**CNVs_CopyKat.html**: 3D UMAP/t-SNE/PCA representation of the copy number variations identified by CopyKat<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**Mean_VAF_RNA**: 3D UMAP/t-SNE/PCA representation of mean VAF for each cell<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**Median_VAF_RNA**: 3D UMAP/t-SNE/PCA representation of median VAF for each cell<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**N_sceSNVs**: 3D UMAP/t-SNE/PCA representation of number of number of SNVs for each cell<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**N_VARreads**: 3D UMAP/t-SNE/PCA representation of number of Variant Reads for each cell<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**N_REFreads**: 3D UMAP/t-SNE/PCA representation of number of Reference Reads for each cell<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**Total_VAF_RNA**: 3D UMAP/t-SNE/PCA representation of Total VAF per cell<br><br>

**sample_SNVs-summary**: a text file of the summary statistics per cell

**Individual_sceSNVs**: contains 3D dmensionality reduction plots for _individual sceSNV_ of the following:
  - Expressed Variant Allele Fraction per cell (VAF_RNA)
  - Number of Reference Reads per cell (REFreads)
  - Number of Variant Reads per cell (VARreads)

**Exploratory_Combined_Plots**: displays all the separately generated plots above into one single html for modularity

<br>

## Installation

Download the R file: 
- [scSNViz.r](https://raw.githubusercontent.com/HorvathLab/NGS/master/scSNViz/scSNViz.r) (Updated: 07/27/2024)

The following CRAN packages are required:
- optparse, stringr, openxlsx, HGNChelper, Seurat, ggplot2, dplyr, plotly, htmlwidgets, htmltools, jsonlite.

The following Bioconductor packages are required:
- glmGamPoi, slingshot.

Other packages:
- copykat 

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
                         -i --tissue-type=Liver -c -d -u  
```
```
% Rscript scSNViz.r -t sample_SNVs.txt -r sample_Seurat_object.rds \
                         --dimensionality-reduction=tsne \
                         --th-vars=1 --th-reads=10 \
                         -i --tissue-type=Immunesystem \
                         -c -d -e --color-scale=YlGnBu 
```
