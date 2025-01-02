## Sample testing snippet for the modules - preprocess_snv_data, plot_snv_data, individual_snv_plots. 
#### The sample input data is in the input directory. 

load.lib<-c("optparse", "stringr", "openxlsx", "HGNChelper", "Seurat", "ggplot2", "dplyr","plotly", "htmlwidgets", "htmltools", "jsonlite", "glmGamPoi", "slingshot")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

```
library(sceSNViz)

rds_file_path <- system.file("input", "sample_Seurat_object.rds", package = "sceSNViz")
snv_file_path <- system.file("input", "sample_SNVs.txt", package = "sceSNViz")

cat("\nRunning preprocessing...\n")

processed_data <- preprocess_snv_data(rds_file = rds_file_path,
                                      snv_file = snv_file_path,
                                      dimensionality_reduction = "UMAP",
                                      th_vars = 0,
                                      th_reads = 0)


cat("\nGenerating plots...\n")

plots <- plot_snv_data(seurat_object = processed_data$SeuratObject,
                       df_3dplot = processed_data$PlotData,
                       save_each_plot = TRUE,
                       output_dir = "output/plots",
                       include_histograms = TRUE,  
                       dimensionality_reduction = "umap",
                       include_cell_types = TRUE,
                       slingshot = TRUE,
                       color_scale = "YlOrRd",
                       cell_size = 0)


cat("\nGenerating individual SNV plots...\n")

snv_plots <- individual_snv_plots(seurat_object = processed_data$SeuratObject,
                                  processed_snv = processed_data$ProcessedSNV,
                                  output_dir = "output/individual_plots",
                                  slingshot = TRUE,
                                  save_each_plot = TRUE,
                                  dimensionality_reduction = "UMAP",
                                  dynamic_cell_size = FALSE)

#generate a report without individual snv plots
generate_report(plot_object = plots)

#generate a report with individual snv plots
generate_report(plot_object = plots,
                snv_plot_object = snv_plots,
                hide_ind_plots = FALSE)
```
