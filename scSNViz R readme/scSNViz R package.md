# Quickstart for Beginners

#### Sample input data is in the input directory. Load libraries, define paths to input files, and define the output directory. 
```
load.lib<-c("scSNViz","SingleCellExperiment", "stringr", "HGNChelper", "Matrix", "umap", "Rtsne", "Seurat", "sctransform", "ggplot2", "readr",
            "dplyr", "plotly", "htmlwidgets", "htmltools", "jsonlite", "glmGamPoi", "slingshot", "copykat", "listviewer","openxlsx") # the installation of ("glmGamPoi") is highly recommended

install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

snv_file <- 'input/sample1_SNVs.tsv'

output_dir = "output"    # or output directory of your choice
```

#### Read in the counts matrix (from either an .RDS file of an existing Seurat object or a counts matrix)
```
#gene.matrix <- Read10X(data.dir = countsmatrix_file) # for reading in a countsmatrix, the data.dir may also be the directory for that contains barcodes.tsv, genes.tsv and matrix.mtx, such as: /user/filtered_gene_bc_matrices/hg19/
#srt <- CreateSeuratObject(counts = gene.matrix, min.cells = 3, min.features = 200)
srt <- readRDS('input/sample1_Seurat_object.rds')
```

#### Quality Control: Filter data and perform scaling and normalization
```
# define the percentage of counts per cell that originate from mitochondrial genes 
srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'

# plot the number of features (or genes) per cell, the number of counts per cell and the percent.mt per cell
VlnPlot(srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # As described in Seurat introductory Vignettes

# filter the Seurat object based on the violin plot
srt <- subset(srt, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA <50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot
```
#### Quality Control: Examples of how filtering impacts the violin plots

##### Unfiltered
[Unfiltered Violin Plot](docs/prefilt_vln.png "UNFILTERED")

##### Filtered
[Filtered Violin PLot](docs/filt_vln.png "FILTERED")


#### Scale and normalize the data. Then, run a PCA.
```
srt <- SCTransform(object = srt, vst.flavor = "v2", method = "glmGamPoi",
           vars.to.regress = "percent.mt", verbose = F)
srt <- RunPCA(srt)
srt <- FindNeighbors(srt, dims = 1:10)
srt <- FindClusters(srt, resolution = 0.5)

```

#### Preprocess the SNV data and incorporate the Seurat object into the workflow
```
# preprocessing the SNV data
processed_data <- preprocess_snv_data(rds_obj = srt,
                                      snv_file = snv_file,
                                      dimensionality_reduction = "UMAP",
                                      th_vars = 0,
                                      th_reads = 0,
                                      enable_sctype = TRUE, #to classify cell types using sctype
                                      tissue_type = "Immunesystem", # other tissue options include: Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus
                                      generate_statistics = TRUE,
                                      output_dir = output_dir)
```

#### Generate the different 3d dimensionality reduction plots
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
                       slingshot = TRUE,
                       color_scale = "YlOrRd",
                       cell_border = 0,
                       save_each_plot = TRUE)
```

#### Generate the report
```
ind_snv_plots <- individual_snv_plots(seurat_object = processed_data$SeuratObject,
                                      processed_snv = processed_data$ProcessedSNV,
                                      output_dir = output_dir,
                                      slingshot = TRUE,
                                      save_each_plot = TRUE,
                                      dimensionality_reduction = "UMAP",
                                      dynamic_cell_size = FALSE)

generate_report(plot_object = plots,
                ind_snv_object = ind_snv_plots,
                hide_ind_plots = TRUE, # Set this to FALSE in order to see plots for each individual SNV
                output_dir = output_dir)
```







# Workflow for Advanced Users with Integrated Samples

####
```
load.lib<-c("scSNViz","SingleCellExperiment", "stringr", "HGNChelper", "Matrix", "umap", "Rtsne", "Seurat", "sctransform", "ggplot2", "readr",
            "dplyr", "plotly", "htmlwidgets", "htmltools", "jsonlite", "glmGamPoi", "slingshot", "copykat", "listviewer","openxlsx") # the installation of ("glmGamPoi") is highly recommended

install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

output_dir = "output"    # or output directory of your choice
```


#### Prepare integrated data.

```
srt1 <- readRDS('sample1_Seurat_object.rds')
srt2 <- readRDS('sample2_Seurat_object.rds')
srt1$orig.ident = 'srt1'
srt2$orig.ident = 'srt2'

## QC ##
srt1[["percent.mt"]] <- PercentageFeatureSet(srt1, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'
VlnPlot(srt1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # As described in Seurat introductory Vignettes
srt1 <- subset(srt1, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA <50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot


srt2[["percent.mt"]] <- PercentageFeatureSet(srt2, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'
VlnPlot(srt2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
srt2 <- subset(srt2, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA <50000 & percent.mt < 15) # Modify numbers appropriate to your violin plot
########

srt_merged <- merge(x=srt1, y=c(srt2), add.cell.ids=c('srt1','srt2'))
srt_merged <- SCTransform(srt_merged, vars.to.regress = "percent.mt", verbose = F)
srt_merged <- RunPCA(srt_merged)
srt_integrated <- IntegrateLayers(object = srt_merged, method = CCAIntegration, normalization.method = "SCT", new.reduction='integrated', verbose = F) #any of the suggested integration methods in Seurat may be applied here for the methods parameter
srt_integrated <- FindNeighbors(srt_integrated, reduction = "integrated", dims = 1:30)
srt_integrated <- FindClusters(srt_integrated, resolution = 0.6)
```

#### Prepare SNV data
```
snv_srt1 <- read.table("sample1_SNVs.tsv", sep = "\t", header = T)
snv_srt2 <- read.table("sample2_SNVs.tsv", sep = "\t", header = T)

snv_srt1$ReadGroup = paste0('srt1_',snv_srt1$ReadGroup) # these IDs must match the added cell IDs from above
snv_srt2$ReadGroup = paste0('srt2_',snv_srt2$ReadGroup) # these IDs must match the added cell IDs from above

snv_file <- rbind(snv_srt1,snv_srt2)
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

#### Generate the different 3d dimensionality reduction plots
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

#### Generate a report
```
ind_snv_plots <- individual_snv_plots(seurat_object = processed_data$SeuratObject,
                                      processed_snv = processed_data$ProcessedSNV,
                                      output_dir = output_dir,
                                      slingshot = TRUE,
                                      save_each_plot = TRUE,
                                      dimensionality_reduction = "UMAP",
                                      dynamic_cell_size = FALSE)

generate_report(plot_object = plots,
                ind_snv_object = ind_snv_plots,
                hide_ind_plots = TRUE, # individual plots for each SNV are hidden
                output_dir = output_dir)
```


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
&nbsp;&nbsp;&nbsp;&nbsp;**Transposed_SNV_Matrix**: 3D UMAP/t-SNE/PCA representation of transposed SNV-Cellbarcode matrix<br><br>

**significant_SNVs.txt**: a text file of significant SNVs (p < 0.05) identified from the statistical significance test (if requested)

**SNV_Statistics.txt**: a text file displaying the results of the statistical significance test (if requested)

**Individual_sceSNVs**: contains 3D dmensionality reduction plots for _individual SNV_ of the following:
  - Expressed Variant Allele Fraction per cell (VAF_RNA)
  - Number of Reference Reads per cell (N_REF)
  - Number of Variant Reads per cell (N_VAR)

**Exploratory_Combined_Plots**: displays all the separately generated plots above into one single html for modularity

<br>
