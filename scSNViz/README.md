# Quickstart for Beginners

#### Sample input data is in the input directory. Load libraries, define paths to input files, and define the output directory. 
```
load.lib<-c("scSNViz","SingleCellExperiment", "stringr", "HGNChelper", "Matrix", "umap", "Rtsne", "Seurat", "ggplot2","readr",
            "dplyr", "plotly", "htmlwidgets", "htmltools", "jsonlite", "glmGamPoi", "slingshot", "copykat", "listviewer") # the installation of ("glmGamPoi") is highly recommended

install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

snv_file <- "sample1_SNVs.tsv"

output_dir = "output"    # or output directory of your choice
```

#### Read in the counts matrix (from either an .RDS file of an existing Seurat object or a counts matrix)
```
#gene.matrix <- Read10X(data.dir = countsmatrix_file) # for reading in a countsmatrix, the data.dir may also be the directory for that contains barcodes.tsv, genes.tsv and matrix.mtx, such as: /user/filtered_gene_bc_matrices/hg19/
#srt <- CreateSeuratObject(counts = gene.matrix, min.cells = 3, min.features = 200)
srt <- readRDS("input/sample1_Seurat_object.rds")
```

#### Quality Control: Filter data and perform scaling and normalization
```
# define the percentage of counts per cell that originate from mitochondrial genes 
srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'

# plot the number of features (or genes) per cell, the number of counts per cell and the percent.mt per cell
VlnPlot(srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # As described in Seurat introductory Vignettes

# filter the Seurat object based on the violin plot
srt <- subset(srt, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) # As described in Seurat introductory Vignettes
```

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

#### Generate report
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
load.lib<-c("scSNViz","SingleCellExperiment", "stringr", "HGNChelper", "Matrix", "umap", "Rtsne", "Seurat", "ggplot2",
            "dplyr", "plotly", "htmlwidgets", "htmltools", "jsonlite", "glmGamPoi", "slingshot", "copykat", "listviewer") # the installation of ("glmGamPoi") is highly recommended

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
srt1 <- subset(srt1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) # As described in Seurat introductory Vignettes

srt2[["percent.mt"]] <- PercentageFeatureSet(srt2, pattern = "^MT-") # this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'
VlnPlot(srt2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
srt2 <- subset(srt2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
########

srt_merged <- merge(x=srt1, y=c(srt2), add.cell.ids=c('srt1','srt2'))
srt_merged[['RNA']]<-split(srt_merged[['RNA']],f=srt_merged$orig.ident)
srt_merged <- SCTransform(srt_merged, vars.to.regress = "percent.mt", verbose = F)
srt_merged <- RunPCA(srt_merged)
srt_integrated <- IntegrateLayers(object = srt_merged, method = CCAIntegration, normalization.method = "SCT", new.reduction='integrated', verbose = F) #any of the suggested integration methods in Seurat may be applied here for the methods parameter
srt_integrated <- FindNeighbors(srt_integrated, reduction = "integrated", dims = 1:30)
srt_integrated <- FindClusters(srt_integrated, resolution = 0.6)
```

#### Prepare SNV data
```
snv_srt1 <- read.table("snv_file_srt1.tsv", sep = "\t", header = T)
snv_srt2 <- read.table("snv_file_srt2.tsv", sep = "\t", header = T)

snv_srt1$ReadGroup = paste0('srt1_',snv_srt1$ReadGroup) # these IDs must match the added cell IDs from above
snv_srt2$ReadGroup = paste0('srt2_',snv_srt2$ReadGroup) # these IDs must match the added cell IDs from above

snv_file <- rbind(snv_srt1,snv_srt2)
write_tsv(snv_file,paste0(output_dir,'/snv_file.tsv'))
```


#### Preprocess the SNV data and incorporate the integrated Seurat object into the workflow. Generate plots.
```
processed_data <- preprocess_snv_data(rds_obj = srt_integrated,
                                      snv_file = paste0(output_dir,'/snv_file.tsv'),
                                      dimensionality_reduction = "UMAP", #you may only generate UMAP plots with integrated samples
                                      th_vars = 0,
                                      th_reads = 0,
                                      enable_integrated = TRUE,
                                      integrated_reduction_name = 'integrated',
                                      enable_sctype = TRUE, # to classify cell types using scType
                                      tissue_type = "Immunesystem", #other tissue options include: Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, Thymus
                                      generate_statistics = TRUE,
                                      output_dir = output_dir)

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
