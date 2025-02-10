## Quickstart for Beginners

#### Sample input data is in the input directory. Load libraries, define paths to input files, and define the output directory. 
```
load.lib<-c("SingleCellExperiment", "stringr", "HGNChelper", "Matrix", "umap", "Rtsne", "Seurat", "ggplot2",
            "dplyr", "plotly", "htmlwidgets", "htmltools", "jsonlite", "glmGamPoi", "slingshot", "copykat", "listviewer") # the installation of ("glmGamPoi") is highly recommended

install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

library(scSNViz)

rds_file <- "path/to/file/.RDS"
snv_file <- "path/to/file/snv.txt"

output_dir = "output"    # or output directory of your choice
```

#### Read in the counts matrix (from either an .RDS file of an existing Seurat object or or a counts matrix)
```
srt <- readRDS(rds_file)

#gene.matrix <- Read10X(data.dir = countsmatrix_file) # for reading in a countsmatrix, the data.dir may also be the directory for that contains barcodes.tsv, genes.tsv and matrix.mtx, such as: /user/filtered_gene_bc_matrices/hg19/
#srt <- CreateSeuratObject(counts = gene.matrix, project = "Sample", # for reading in 
                  min.cells = 3, min.features = 200)
```

#### Quality Control: Filter data and perform scaling and normalization
```
# define the percentage of counts per cell that originate from mitochondrial genes 
srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-") #this is for Homo sapiens. If the organism is Mus musculus, then: pattern = '^mt-'

# plot the number of features (or genes) per cell, the number of counts per cell and the percent.mt per cell
VlnPlot(srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #As described in Seurat introductory Vignettes

# filter the Seurat object based on the violin plot
srt <- subset(srt, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #As described in Seurat introductory Vignettes
```

#### Scale and normalize the data
```
srt <- SCTransform(object = srt, vst.flavor = "v2", method = "glmGamPoi",
           vars.to.regress = "percent.mt", verbose = F)

```

#### Preprocess the SNV data and incorporate the Seurat object into the workflow
```
# preprocessing the SNV data
processed_data <- preprocess_snv_data(rds_file = srt,
                                      snv_file = snv_file,
                                      dimensionality_reduction = "UMAP",
                                      th_vars = 0,
                                      th_reads = 0,
                                      enable_sctype = TRUE, #to classify cell types using sctype
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
                       include_copykat = FALSE, #CNV metrics produced by copykat; this may significantly increase processing time depending on the size of gene counts matrix provided
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

## Workflow for Advanced Users with Integrated Samples

#### Use integrated data to generate plots; the Default Assay for this example is 'integrated'. Suggested tutorial: ()

```
srt1 <- CreateSeuratObject(counts = gene.matrix1, min.cells = 3, min.features = 200)
srt2 <- CreateSeuratObject(counts = gene.matrix2, min.cells = 3, min.features = 200)
srt3 <- CreateSeuratObject(counts = gene.matrix3, min.cells = 3, min.features = 200)

srt_merged <- merge(x=srt1, y=c(srt2,srt3), add.cell.ids=c('srt1','srt2','srt3'))
srt_integrated <- IntegrateLayers(object = srt_merged, method = CCAIntegration, normalization.method = "SCT", new.reduction='integrated', verbose = F)
srt_integrated <- FindNeighbors(srt_integrated, reduction = "integrated", dims = 1:30)
srt_integrated <- FindClusters(srt_integrated, resolution = 0.6)
```

#### Preprocess the SNV data and incorporate the integrated Seurat object into the workflow. Generate plots.
```
processed_data <- preprocess_snv_data(rds_file = integrated,
                                      snv_file = snv_file,
                                      dimensionality_reduction = "UMAP",
                                      th_vars = 0,
                                      th_reads = 0,
                                      enable_integrated = TRUE,
                                      enable_sctype = TRUE, #to classify cell types using scType
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
                       slingshot = TRUE,
                       color_scale = "YlOrRd",
                       cell_border = 0,
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
