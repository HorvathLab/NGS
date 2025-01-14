## Sample testing snippet for the modules - preprocess_snv_data, plot_snv_data, individual_snv_plots, and generate_report. 
#### The sample input data is in the input directory. Provide an optional output directory (argument output_dir) for the package output. 

```
load.lib<-c("SingleCellExperiment", "stringr", "HGNChelper", "Matrix", "umap", "Rtsne", "Seurat", "ggplot2",
            "dplyr", "plotly", "htmlwidgets", "htmltools", "jsonlite", "glmGamPoi", "slingshot", "copykat", "listviewer")

install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

library(sceSNViz)

rds_file_path <- system.file("input", "sample_Seurat_object.rds", package = "sceSNViz")
snv_file_path <- system.file("input", "sample_SNVs.txt", package = "sceSNViz")

#if above does not work:
rds_file_path <- "path/to/file/.RDS"
snv_file_path <- "path/to/file/snv.txt"

output_dir = "output"    # or output directory of your choice

# preprocessing the snv data
processed_data <- preprocess_snv_data(rds_file = rds_file_path,
                                      snv_file = snv_file_path,
                                      dimensionality_reduction = "UMAP",
                                      th_vars = 0,
                                      th_reads = 0,
                                      enable_sctype = TRUE,
                                      tissue_type = "Immunesystem",
                                      generate_statistics = TRUE,
                                      output_dir = output_dir)

# generate the different 3d dimensionality reduction plots
plots <- plot_snv_data(seurat_object = processed_data$SeuratObject,
                       plot_data = processed_data$PlotData,
                       processed_snv = processed_data$ProcessedSNV,
                       output_dir = output_dir,
                       include_histograms = TRUE,  
                       dimensionality_reduction = "umap",
                       include_cell_types = TRUE,
                       include_copykat = FALSE,
                       slingshot = TRUE,
                       color_scale = "YlOrRd",
                       save_each_plot = TRUE)

# generate individual snv plots
ind_snv_plots <- individual_snv_plots(seurat_object = processed_data$SeuratObject,
                                      processed_snv = processed_data$ProcessedSNV,
                                      output_dir = output_dir,
                                      slingshot = TRUE,
                                      save_each_plot = TRUE,
                                      dimensionality_reduction = "UMAP",
                                      dynamic_cell_size = FALSE)

#generate a report without individual snv plots
generate_report(plot_object = plots)

#generate a report with individual snv plots
generate_report(plot_object = plots,
                snv_plot_object = ind_snv_plots,
                hide_ind_plots = FALSE,
                output_dir = output_dir)
```
