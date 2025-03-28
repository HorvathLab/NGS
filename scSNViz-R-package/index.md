# scSNViz: visualization and analysis of Cell-Specific Expressed SNVs
Understanding expressed genetic variation at the single-cell level is crucial for insights into
transcriptional heterogeneity and gene expression regulation, but there is a scarcity of tools for visualizing
and analyzing cell-level genetic variants. scSNViz is a dedicated tool for the visualization and analysis of cell-specific expressed
single nucleotide variants (sceSNVs) in cell-barcoded single-cell RNA sequencing (scRNA-seq) data. It
integrates established frameworks such as Seurat, Slingshot, scType, and CopyKat into a specialized
workflow for sceSNV expression analysis. scSNViz enables quantitative 3D visualization of variant allele
frequencies (VAF_RNA), reference (N_REF) and variant (N_VAR) read counts, while also supporting SNV
clustering based on shared transcriptional activity. Beyond visualization, scSNViz provides estimation,
summarization, and graphical representation of SNV expression metrics, facilitating the study of allelic
dynamics, somatic, germline, and RNA-originating SNVs, transcriptional bursting and lineage-specific
expression, and. Additionally, scSNViz supports both individual and set-based sceSNV analyses, as well as
comparative assessments across multiple samples, making it a powerful tool for understanding SNV-driven
regulatory mechanisms and cellular heterogeneity.

#<img src='https://github.com/HorvathLab/NGS/blob/scSNViz_R_v1.0.0/scSNViz/docs/scSNViz_PanelA.png' width=50% height=50%>
![scSNViz_PanelA](scSNViz-R-package/assets/scSNViz_PanelA.png)

## Installation

#### Install scSNViz from GitHub

```
library(devtools)
install_github("HorvathLab/NGS", ref = "scSNViz_R_v1.0.0", subdir = "scSNViz")
```
If the above fails due to rate limits, try generating a GitHub Personal Access Token (PAT), add it into your environment and then run again. 

Another way to do this is to configure R to use the Windows Internet API for download: 

```
options(download.file.method = "wininet")   # can try other methods such as 'libcurl', 'wget', etc.
install_github("HorvathLab/NGS", ref = "scSNViz_R_v1.0.0", subdir = "scSNViz")
```

## Quickstart for Beginners

#### Load libraries, define paths to input files, and define the output directory.
The input files are located in the input folder on github. The snv file is an output from SCReadCounts. The user may provide a .tsv file that is not from SCReadCounts as long as it is also a .tsv and contains the following columns: CHROM, POS, REF, ALT, ReadGroup, SNVCount, RefCount.
```
load.lib<-c("scSNViz","SingleCellExperiment", "stringr", "HGNChelper", "Matrix", "umap", "Rtsne", "Seurat", "sctransform", "ggplot2", "readr",
            "dplyr", "plotly", "htmlwidgets", "htmltools", "jsonlite", "glmGamPoi", "slingshot", "listviewer","openxlsx","randomcoloR") # the installation of ("glmGamPoi") is highly recommended

install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

#CopyKat is an optional tool in analysis and must be installed separately
library(devtools)
install_github("navinlabcode/copykat")

snv_file <- 'input/sample1_SNVs.tsv'
srt_obj_file <- 'input/sample1_Seurat_object.rds'
output_dir = "output"    # or output directory of your choice
```

#### Read in the counts matrix (from either an .RDS file of an existing Seurat object or a counts matrix)
```
#gene.matrix <- Read10X(data.dir = countsmatrix_file) # for reading in a countsmatrix, the data.dir may also be the directory for that contains barcodes.tsv, genes.tsv and matrix.mtx, such as: /user/filtered_gene_bc_matrices/hg19/
#sample1 <- CreateSeuratObject(counts = gene.matrix, min.cells = 3, min.features = 200, project = 'Sample1')
sample1 <- readRDS(srt_obj_file)
sample1@project.name = 'Sample1' # set the project name
sample1$orig.ident = 'Sample1' # set the project name
```

#### Quality Control: Filter data and perform scaling and normalization
```
# define the percentage of counts per cell that originate from mitochondrial genes 
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'

# plot the number of features (or genes) per cell, the number of counts per cell and the percent.mt per cell
VlnPlot(sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # As described in Seurat introductory Vignettes

# filter the Seurat object based on the violin plot
sample1 <- subset(sample1, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA <50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot
```

The below unfiltered and filtered violin plots show the quality control process, filtering cells out based on percent mitochondria, number of features and number of counts.

<img src='https://github.com/HorvathLab/NGS/blob/scSNViz_R_v1.0.0/scSNViz/docs/prefilt_filt_vln.png' width=50% height=50%>


#### Scale and normalize the data. Then, run a PCA.
```
sample1 <- SCTransform(object = sample1, vst.flavor = "v2", method = "glmGamPoi",
           vars.to.regress = "percent.mt", verbose = F)
sample1 <- RunPCA(sample1)
sample1 <- FindNeighbors(sample1, dims = 1:10)
sample1 <- FindClusters(sample1, resolution = 0.5)

```

#### Preprocess the SNV data and incorporate the Seurat object into the workflow
```
# preprocessing the SNV data
processed_data <- preprocess_snv_data(rds_obj = sample1,
                                      snv_file = snv_file,
                                      dimensionality_reduction = "UMAP",
                                      th_vars = 0,
                                      th_reads = 0,
                                      enable_sctype = TRUE, #to classify cell types using sctype
                                      tissue_type = "Immunesystem", # other tissue options include: Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus
                                      generate_statistics = TRUE,
                                      output_dir = output_dir)
```

#### Generate 3d dimensionality reduction plots
```
plots <- plot_snv_data(seurat_object = processed_data$SeuratObject,
                       processed_data$ProcessedSNV,
                       processed_data$AggregatedSNV,
                       processed_data$PlotData,
                       output_dir = output_dir,
                       include_histograms = TRUE,  
                       dimensionality_reduction = "umap",
                       include_cell_types = TRUE,
                       include_copykat = FALSE, # CNV metrics produced by copykat; this may significantly increase processing time depending on the size of gene counts matrix provided
                       include_snv_dim_red = FALSE, # IF, set to TRUE, this function transposes the SNVxBarcode matrix and generates a dimensionality reduction plot to view similarity between SNVs.
                       slingshot = TRUE,
                       color_scale = "YlOrRd",
                       cell_border = 0,
                       save_each_plot = TRUE)
```

<img src='https://github.com/HorvathLab/NGS/blob/scSNViz_R_v1.0.0/scSNViz/docs/sample_outputs.png'>

#### Generate individual SNV plots
```
#Individual SNV's plottable capped at 50 unique.
ind_snv_plots <- individual_snv_plots(seurat_object = processed_data$SeuratObject,
                                      processed_snv = processed_data$ProcessedSNV,
                                      output_dir = output_dir,
                                      slingshot = TRUE,
                                      save_each_plot = TRUE,
                                      dimensionality_reduction = "UMAP",
                                      dynamic_cell_size = FALSE)
```

<img src='https://github.com/HorvathLab/NGS/blob/scSNViz_R_v1.0.0/scSNViz/docs/individual_snv_plots.png'>

#### Plot single SNV
```
one_snv_plot <- single_snv_plot(
       seurat_object = processed_data$SeuratObject,
       processed_snv = processed_data$ProcessedSNV,
       snv_of_choice = "1:155169447:C:T",
       output_dir = "output/1_155169447_C_T",
       slingshot = TRUE,
       dimensionality_reduction = "UMAP",
       dynamic_cell_size = FALSE,
       save_each_plot = TRUE
     )

```

#### Generate exploratory combined plots report
```
generate_report(plot_object = plots,
                ind_snv_object = ind_snv_plots,
                hide_ind_plots = TRUE, # Set this to FALSE in order to see plots for each individual SNV.
                output_dir = output_dir)
```


<img src='https://github.com/HorvathLab/NGS/blob/scSNViz_R_v1.0.0/scSNViz/docs/Exploratory_combined_plots.png'>


#### Generate exploratory combined plot for single SNV of interest
```
generate_report(plot_object = plots,
                ind_snv_object = one_snv_plot,
                hide_ind_plots = FALSE,
                output_dir = output_dir)
```





## Workflow for Advanced Users with Multiple Samples to Integrate

####
```
load.lib<-c("scSNViz","SingleCellExperiment", "stringr", "HGNChelper", "Matrix", "umap", "Rtsne", "Seurat", "sctransform", "ggplot2", "readr",
            "dplyr", "plotly", "htmlwidgets", "htmltools", "jsonlite", "glmGamPoi", "slingshot", "copykat", "listviewer","openxlsx","randomcoloR") # the installation of ("glmGamPoi") is highly recommended

install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

output_dir = "output_advanced_wkflw"    # or output directory of your choice
```


#### Prepare integrated data.

```
sample1 <- readRDS('input/sample1_Seurat_object.rds')
sample2 <- readRDS('input/sample2_Seurat_object.rds')
sample1$orig.ident = 'sample1'
sample2$orig.ident = 'sample2'

## QC ##
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'
VlnPlot(sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # As described in Seurat introductory Vignettes
sample1 <- subset(sample1, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA <50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot


sample2[["percent.mt"]] <- PercentageFeatureSet(sample2, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'
VlnPlot(sample2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sample2 <- subset(sample2, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA <50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot
########

srt_merged <- merge(x=sample1, y=c(sample2), add.cell.ids=c('sample1','sample2'))
srt_merged <- SCTransform(srt_merged, vars.to.regress = "percent.mt", verbose = F)
srt_merged <- RunPCA(srt_merged)
srt_integrated <- IntegrateLayers(object = srt_merged, method = CCAIntegration, normalization.method = "SCT", new.reduction='integrated', verbose = F) #any of the suggested integration methods in Seurat may be applied here for the methods parameter
srt_integrated <- FindNeighbors(srt_integrated, reduction = "integrated", dims = 1:30)
srt_integrated <- FindClusters(srt_integrated, resolution = 0.6)
```

#### Prepare SNV data
```
snv_sample1 <- read.table("input/sample1_SNVs.tsv", sep = "\t", header = T)
snv_sample2 <- read.table("input/sample2_SNVs.tsv", sep = "\t", header = T)

snv_sample1$ReadGroup = paste0('sample1_',snv_sample1$ReadGroup) # these IDs must match the added cell IDs from above
snv_sample2$ReadGroup = paste0('sample2_',snv_sample2$ReadGroup) # these IDs must match the added cell IDs from above

snv_file <- rbind(snv_sample1,snv_sample2)
write_tsv(snv_file,'snv_file_integrated.tsv')
```


#### Preprocess the SNV data and incorporate the integrated Seurat object into the workflow. Generate plots.
```
processed_data <- preprocess_snv_data(rds_obj = srt_integrated,
                                      snv_file = "snv_file_integrated.tsv",
                                      dimensionality_reduction = "UMAP", #you may only generate UMAP plots with integrated samples
                                      th_vars = 0,
                                      th_reads = 0,
                                      enable_integrated = TRUE,
                                      integrated_reduction_name = 'integrated',
                                      enable_sctype = TRUE, # to classify cell types using scType
                                      tissue_type = "Immunesystem", #other tissue options include: Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus
                                      generate_statistics = TRUE,
                                      output_dir = output_dir)
```

#### Generate 3d dimensionality reduction plots
```
plots <- plot_snv_data(seurat_object = processed_data$SeuratObject,
                       processed_data$ProcessedSNV,
                       processed_data$AggregatedSNV,
                       processed_data$PlotData,
                       output_dir = output_dir,
                       include_histograms = TRUE,  
                       dimensionality_reduction = "umap",
                       include_cell_types = TRUE,
                       include_copykat = FALSE, # this is not recommended for integrated objects
                       slingshot = FALSE, # this is default, no slingshot option for integrated data
                       color_scale = "YlOrRd",
                       cell_border = 0,
                       enable_integrated = TRUE,
                       save_each_plot = TRUE)
```

<img src='https://github.com/HorvathLab/NGS/blob/scSNViz_R_v1.0.0/scSNViz/docs/integrated_plot.png'>

#### Generate individual SNV plots
```
ind_snv_plots <- individual_snv_plots(seurat_object = processed_data$SeuratObject,
                                      processed_snv = processed_data$ProcessedSNV,
                                      output_dir = output_dir,
                                      slingshot = FALSE,
                                      save_each_plot = TRUE,
                                      dimensionality_reduction = "UMAP",
                                      dynamic_cell_size = FALSE,
                                      enable_integrated = TRUE)
```

#### Generate exploratory combined plots report
```
generate_report(plot_object = plots,
                ind_snv_object = ind_snv_plots,
                hide_ind_plots = TRUE, # individual plots for each SNV are hidden
                output_dir = output_dir)
```

<img src='https://github.com/HorvathLab/NGS/blob/scSNViz_R_v1.0.0/scSNViz/docs/integrated_output_example.png'>


## Output

### Directory structure

scSNViz generates outputs for the set of the scSNVs and for each individual scSNV, as follows:

```
output_dir/
│
├── Exploratory_Combined_Plots.html
├── significant_SNVs.txt             (optional)
├── SNV_Statistics.txt               (optional)
└── SNV_data_plots/   
    ├── Cell_types_scType.html       (optional)
    ├── CNVs_CopyKat.html            (optional)
    ├── Median_VAF_RNA.html
    ├── Mean_VAF_RNA.html
    ├── N_REFreads.html
    ├── Total_VAF_RNA.html
    ├── N_sceSNVs.html
    ├── N_VARreads.html
    ├── SAMPLE_ID_plot.html          (optional)
    ├── Transposed_SNV_Matrix.html
    ├── Histogram_N_SNV.png
    ├── Histogram_N_VARreadsCounts.png
    ├── Histogram_MeanSNVsVAF.png
    ├── Histogram_TotalVAF.png
└── Individual_sceSNVs/
    ├── VARreads/
    │   └── 3D N_VAR plot HTML files for each scSNV
    ├── REFreads/
    │   └── 3D N_REF plot HTML files for each scSNV
    └── VAF/
        └── 3D VAF plot HTML files for each scSNV
```

## Description

For the set of the scSNVs, the separately produced figures show the following:

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
&nbsp;&nbsp;&nbsp;&nbsp;**SAMPLE_ID_plot**: 3D UMAP/t-SNE/PCA representation of the sample IDs<br><br>
&nbsp;&nbsp;&nbsp;&nbsp;**Transposed_SNV_Matrix**: 3D UMAP/t-SNE/PCA representation of transposed SNV-Cellbarcode matrix<br><br>

**significant_SNVs.txt**: a text file of significant SNVs (p < 0.05) identified from the statistical significance test (if requested)

**SNV_Statistics.txt**: a text file displaying the results of the statistical significance test (if requested)

**Individual_sceSNVs**: contains 3D dmensionality reduction plots for _individual SNV_ of the following:
  - Expressed Variant Allele Fraction per cell (VAF_RNA)
  - Number of Reference Reads per cell (N_REF)
  - Number of Variant Reads per cell (N_VAR)

**Exploratory_Combined_Plots**: displays all the separately generated plots above into one single html for modularity

<br>
