#' Plot SNV Data
#'
#' This function generates interactive 3D dimensionality reduction plots, histograms, Cell Types Plot, and CopyKat Plot, with the option to save plots individually.
#'
#' @importFrom Seurat Embeddings as.SingleCellExperiment
#' @importFrom ggplot2 ggplot aes geom_histogram xlab ylab theme element_text element_blank
#' @importFrom dplyr %>% filter group_by
#' @importFrom plotly add_markers plot_ly add_trace layout as_widget
#' @importFrom htmlwidgets saveWidget
#'
#' @param seurat_object Processed Seurat object with SNV metadata.
#' @param df_3dplot Data frame for 3D plotting.
#' @param output_dir Directory to save the plots (required if save_each_plot is TRUE).
#' @param include_histograms Logical; whether to generate histograms. Default: TRUE.
#' @param include_cell_types Logical; whether to generate cell types plot. Default: FALSE.
#' @param include_copykat Logical; whether to generate CopyKat plot. Default: FALSE.
#' @param dimensionality_reduction Dimensionality reduction method used ('UMAP', 'PCA', or 'tSNE').
#' @param slingshot Logical; whether to include slingshot trajectories in 3D plots. Default: TRUE.
#' @param color_scale Color scale for plots. Default: 'YlOrRd'.
#' @param cell_size Numeric; size of points in the plot. Default: 3.
#' @param disable_3d_axis Logical; whether to disable 3D axis labels. Default: FALSE.
#' @param save_each_plot Logical; whether to save each plot individually. Default: FALSE.
#' @return A list of generated plots.
#' @details
#' The function generates various visualizations for SNV data including:
#' - **3D Plots**: Visualizes metrics such as SNV.N, SNVCount, and TotalVAF.
#' - **Histograms**: Distribution of metrics including TotalVAF, MeanSNVsVAF, and N_VARreadCounts.
#' - **Cell Types Plot**: Shows cell types (e.g., custom classifications) with optional slingshot trajectories.
#' - **CopyKat Plot**: Depicts Copy Number Variations (CNVs) using CopyKat analysis.
#'
#' The output_dir parameter specifies where the plots will be saved if save_each_plot is enabled.
#'
#' @examples
#' # Example usage:
#' plots <- plot_snv_data(
#'   seurat_object = processed_data$SeuratObject,
#'   df_3dplot = processed_data$PlotData,
#'   output_dir = "output/plots",
#'   include_histograms = TRUE,
#'   include_cell_types = TRUE,
#'   include_copykat = FALSE,
#'   dimensionality_reduction = "UMAP",
#'   slingshot = TRUE,
#'   color_scale = "YlOrRd",
#'   cell_size = 3,
#'   disable_3d_axis = FALSE,
#'   save_each_plot = TRUE
#' )
#'
#' # Access individual plots:
#' plot_vaf <- plots$VAF
#' plot_snvcount <- plots$SNVCount
#'
#' @export
#'
plot_snv_data <- function(seurat_object, df_3dplot, output_dir = NULL, include_histograms = TRUE,
                          include_cell_types = FALSE, include_copykat = FALSE, dimensionality_reduction = "UMAP",
                          slingshot = TRUE, color_scale = "YlOrRd", cell_size = 3,
                          disable_3d_axis = FALSE, save_each_plot = FALSE) {

  valid_reductions <- c("umap", "pca", "tsne")
  dimensionality_reduction <- tolower(dimensionality_reduction)
  if (!dimensionality_reduction %in% valid_reductions) {
    stop("Invalid dimensionality_reduction method. Please use one of: 'umap', 'pca', 'tsne'.")
  }
  scale_list <- c("Blues", "Reds", "YlOrRd", "YlGnBu", "plasma",
                  "RdBu")
  reversescale_options <- c(TRUE, FALSE, TRUE, TRUE, FALSE,
                            FALSE)
  if (!color_scale %in% scale_list) {
    stop("Invalid color_scale. Please choose from: ", paste(scale_list,
                                                            collapse = ", "))
  }
  reversescale_option <- reversescale_options[which(scale_list ==
                                                      color_scale)]
  histogram_scale1 <- "tomato"
  histogram_scale2 <- "dodgerblue2"
  color_undetected <- "lightgrey"
  plots <- list()
  if (save_each_plot && !is.null(output_dir)) {
    dir.create(file.path(output_dir, "Figures_Individual_Plots_HTML"),
               showWarnings = FALSE, recursive = TRUE)
  }
  curves <- NULL
  if (slingshot) {
    sce <- as.SingleCellExperiment(seurat_object)
    sce <- slingshot(sce, clusterLabels = "seurat_clusters",
                     reducedDim = toupper(dimensionality_reduction))
    curves <- slingCurves(sce, as.df = TRUE)
    colnames(curves)[1:3] <- c(paste0(toupper(dimensionality_reduction),
                                      "_1"), paste0(toupper(dimensionality_reduction),
                                                    "_2"), paste0(toupper(dimensionality_reduction),
                                                                  "_3"))
  }
  histograms <- list()
  if (include_histograms) {
    hist_list <- list(list(aes = aes(x = TotalVAF), xlab = "TotalVAF",
                           file_suffix = "Histogram_TotalVAF"), list(aes = aes(x = MeanVAF),
                                                                     xlab = "MeanSNVsVAF", file_suffix = "Histogram_MeanSNVsVAF"),
                      list(aes = aes(x = SNVCount), xlab = "N_VARreadCounts",
                           file_suffix = "Histogram_N_VARreadCounts"), list(aes = aes(x = SNV.N),
                                                                            xlab = "N_SNV", file_suffix = "Histogram_N_SNV"))
    for (hist_info in hist_list) {
      p <- ggplot(seurat_object@meta.data, hist_info$aes) +
        geom_histogram(fill = histogram_scale2, boundary = 0,
                       alpha = 0.6) + xlab(hist_info$xlab) + ylab("Cells") +
        theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
              panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
              axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
      histograms[[hist_info$xlab]] <- ggplotly(p)
      if (save_each_plot && !is.null(output_dir)) {
        ggsave(file = file.path(output_dir, "Figures_Individual_Plots_HTML",
                                paste0(hist_info$file_suffix, ".png")), plot = p,
               device = "png")
      }
    }
  }
  plot_descriptions <- list(SNV.N = "Number of Detected SNVs per Cell",
                            SNVCount = "Total Variant Reads per Cell", 
                            RefCount = "Total Reference Reads per Cell",
                            TotalVAF = "Total VAF (Variant Allele Fraction) per Cell",
                            MeanVAF = "Mean VAF across All Loci per Cell", 
                            MedianVAF = "Median VAF across All Loci per Cell")
  dimensionality_reduction <- tolower(dimensionality_reduction)
  dim_plotting <- toupper(dimensionality_reduction)


generate_plot <- function(metric, df_3dplot, dim_plotting, color_scale, reversescale_option, 
                          cell_size, color_undetected, curves, disable_3d_axis) {
  plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
    add_markers(
      data = df_3dplot[df_3dplot[["Undetected"]] == 0, ],
      x = ~get(paste0(dim_plotting, "_1")),
      y = ~get(paste0(dim_plotting, "_2")),
      z = ~get(paste0(dim_plotting, "_3")),
      marker = list(
        color = ~get(metric),
        colorscale = color_scale,
        reversescale = reversescale_option,
        size = cell_size,
        showscale = TRUE,
        colorbar = list(x = 0.85, y = 0.5, thickness = 20, len = 0.5),
        line = list(color = ~get(metric), width = cell_size)
      ),
      opacity = 0.8,
      name = "expressed sceSNV loci"
    ) %>%
    add_markers(
      data = df_3dplot[df_3dplot[["Undetected"]] == 1, ],
      x = ~get(paste0(dim_plotting, "_1")),
      y = ~get(paste0(dim_plotting, "_2")),
      z = ~get(paste0(dim_plotting, "_3")),
      marker = list(color = color_undetected, size = cell_size, opacity = 0.5),
      name = "undetected"
    )
  
  if (!is.null(curves)) {
    plot <- plot %>%
      add_trace(
        data = curves,
        x = ~get(paste0(dim_plotting, "_1")),
        y = ~get(paste0(dim_plotting, "_2")),
        z = ~get(paste0(dim_plotting, "_3")),
        split = ~Lineage,
        mode = "lines",
        line = list(width = 2)
      )
  }
  
  plot <- plot %>%
    layout(
      scene = list(
        xaxis = list(title = paste0(dim_plotting, "_1")),
        yaxis = list(title = paste0(dim_plotting, "_2")),
        zaxis = list(title = paste0(dim_plotting, "_3"))
      )
    )
  
  if (disable_3d_axis) {
    plot <- plot %>%
      layout(
        scene = list(
          xaxis = list(title = NULL, showticklabels = FALSE),
          yaxis = list(title = NULL, showticklabels = FALSE),
          zaxis = list(title = NULL, showticklabels = FALSE)
        )
      )
  }
  
  return(plot)
}


plots <- list() 

for (metric in names(plot_descriptions)) {
  if (!metric %in% colnames(df_3dplot))
    next
  
  plot <- generate_plot(
    metric = metric, 
    df_3dplot = df_3dplot, 
    dim_plotting = dim_plotting, 
    color_scale = color_scale, 
    reversescale_option = reversescale_option, 
    cell_size = cell_size, 
    color_undetected = color_undetected, 
    curves = curves, 
    disable_3d_axis = disable_3d_axis
  )
  
  plots[[metric]] <- plot
  
  if (save_each_plot && !is.null(output_dir)) {
    saveWidget(
      as_widget(plot),
      file = file.path(output_dir, "Figures_Individual_Plots_HTML", paste0(metric, ".html")),
      selfcontained = FALSE,
      libdir = "lib"
    )
  }
}


  if (include_cell_types && "customclassif" %in% colnames(df_3dplot)) {
    cell_type_plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
      add_trace(data = df_3dplot, x = ~get(paste0(toupper(dimensionality_reduction),
                                                  "_1")), y = ~get(paste0(toupper(dimensionality_reduction),
                                                                          "_2")), z = ~get(paste0(toupper(dimensionality_reduction),
                                                                                                  "_3")), size = 0.05, opacity = 0.5, mode = "markers",
                color = ~customclassif)
    if (slingshot && !is.null(curves)) {
      cell_type_plot <- cell_type_plot %>% add_trace(data = curves,
                                                     x = ~get(paste0(toupper(dimensionality_reduction),
                                                                     "_1")), y = ~get(paste0(toupper(dimensionality_reduction),
                                                                                             "_2")), z = ~get(paste0(toupper(dimensionality_reduction),
                                                                                                                     "_3")), mode = "lines", color = ~factor(Lineage),
                                                     line = list(width = 4))
    }
    if (disable_3d_axis) {
      cell_type_plot <- cell_type_plot %>% layout(scene = list(xaxis = list(title = NULL,
                                                                            showticklabels = FALSE, zeroline = FALSE, showline = FALSE,
                                                                            showgrid = FALSE), yaxis = list(title = NULL,
                                                                                                            showticklabels = FALSE, zeroline = FALSE, showline = FALSE,
                                                                                                            showgrid = FALSE), zaxis = list(title = NULL,
                                                                                                                                            showticklabels = FALSE, zeroline = FALSE, showline = FALSE,
                                                                                                                                            showgrid = FALSE)))
    }
    else {
      cell_type_plot <- cell_type_plot %>% layout(scene = list(title = "",
                                                               xaxis = list(title = paste0(toupper(dimensionality_reduction),
                                                                                           "_1")), yaxis = list(title = paste0(toupper(dimensionality_reduction),
                                                                                                                               "_2")), zaxis = list(title = paste0(toupper(dimensionality_reduction),
                                                                                                                                                                   "_3"))))
    }
    plots[["Cell Types Plot"]] <- cell_type_plot
    if (save_each_plot && !is.null(output_dir)) {
      cell_type_plot <- cell_type_plot %>% layout(title = "Cell types (scType)",
                                                  title = list(font = "black"), margin = list(t = 50))
      saveWidget(as_widget(cell_type_plot), file = file.path(output_dir,
                                                             "Figures_Individual_Plots_HTML", "Cell_types_scType.html"),
                 selfcontained = F, libdir = "lib")
    }
    if (slingshot) {
      s <- slingshot(as.SingleCellExperiment(seurat_object),
                     clusterLabels = "seurat_clusters", reducedDim = toupper(dimensionality_reduction))
      lineage_counter <- 1
      for (i in grep("slingPseudotime", colnames(s@colData))) {
        lineage_id <- as.integer(sub("slingPseudotime_",
                                     "", colnames(s@colData)[i]))
        curve <- curves[curves$Lineage == lineage_id,
        ]
        lineage_plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
          add_trace(data = df_3dplot, x = ~get(paste0(toupper(dimensionality_reduction),
                                                      "_1")), y = ~get(paste0(toupper(dimensionality_reduction),
                                                                              "_2")), z = ~get(paste0(toupper(dimensionality_reduction),
                                                                                                      "_3")), marker = list(color = ~s@colData[,
                                                                                                                                               i], colorscale = "YlOrRd", reversescale = reversescale_option,
                                                                                                                            colorbar = list(x = 0), line = list(color = ~SNV.N,
                                                                                                                                                                width = cell_size)))
        if (!disable_3d_axis) {
          lineage_plot <- lineage_plot %>% layout(scene = list(xaxis = list(title = paste0(toupper(dimensionality_reduction),
                                                                                           "_1")), yaxis = list(title = paste0(toupper(dimensionality_reduction),
                                                                                                                               "_2")), zaxis = list(title = paste0(toupper(dimensionality_reduction),
                                                                                                                                                                   "_3")))) %>% add_trace(data = curve, x = ~get(paste0(toupper(dimensionality_reduction),
                                                                                                                                                                                                                        "_1")), y = ~get(paste0(toupper(dimensionality_reduction),
                                                                                                                                                                                                                                                "_2")), z = ~get(paste0(toupper(dimensionality_reduction),
                                                                                                                                                                                                                                                                        "_3")), mode = "lines", showlegend = TRUE,
                                                                                                                                                                                          line = list(width = 4, color = "black"))
        }
        else {
          lineage_plot <- lineage_plot %>% layout(scene = list(xaxis = list(title = NULL,
                                                                            showticklabels = FALSE, zeroline = FALSE,
                                                                            showline = FALSE, showgrid = FALSE), yaxis = list(title = NULL,
                                                                                                                              showticklabels = FALSE, zeroline = FALSE,
                                                                                                                              showline = FALSE, showgrid = FALSE), zaxis = list(title = NULL,
                                                                                                                                                                                showticklabels = FALSE, zeroline = FALSE,
                                                                                                                                                                                showline = FALSE, showgrid = FALSE)))
        }
      }
      lineage_counter <- lineage_counter + 1
    }
  }
  if (include_copykat && "karyotype" %in% colnames(seurat_object@meta.data)) {
    cat("Running CopyKat. This may take a while...\n")
    cts <- GetAssayData(seurat_object, layer = "counts",
                        assay = "SCT")
    ckt <- copykat(cts, sam.name = "sample_name")
    seurat_object[["karyotype"]] <- "Unknown"
    seurat_object[["karyotype"]][rownames(ckt$pred), ] <- ckt$pred[,
                                                                   2]
    df.srt <- as.data.frame(seurat_object@meta.data[, c("seurat_clusters",
                                                        "karyotype")])
    df.embed <- as.data.frame(Embeddings(seurat_object, reduction = dimensionality_reduction))
    df.3dplot <- merge(df.srt, df.embed, by = 0)
    colnames(df.3dplot) <- c("filler", "seurat_clusters",
                             "karyotype", paste0(toupper(dimensionality_reduction),
                                                 "_1"), paste0(toupper(dimensionality_reduction),
                                                               "_2"), paste0(toupper(dimensionality_reduction),
                                                                             "_3"))
    rownames(df.3dplot) <- df.3dplot$Row.names
    df.3dplot <- df.3dplot[, 2:ncol(df.3dplot)]
    copykat_plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
      add_markers(data = df.3dplot, x = ~get(paste0(toupper(dimensionality_reduction),
                                                    "_1")), y = ~get(paste0(toupper(dimensionality_reduction),
                                                                            "_2")), z = ~get(paste0(toupper(dimensionality_reduction),
                                                                                                    "_3")), size = 0.05, opacity = 1, color = ~karyotype)
    if (slingshot && !is.null(curves)) {
      copykat_plot <- copykat_plot %>% add_trace(data = curves,
                                                 x = ~get(paste0(toupper(dimensionality_reduction),
                                                                 "_1")), y = ~get(paste0(toupper(dimensionality_reduction),
                                                                                         "_2")), z = ~get(paste0(toupper(dimensionality_reduction),
                                                                                                                 "_3")), split = ~Lineage, mode = "lines")
    }
    if (disable_3d_axis) {
      copykat_plot <- copykat_plot %>% layout(scene = list(xaxis = list(title = NULL,
                                                                        showticklabels = FALSE, zeroline = FALSE, showline = FALSE,
                                                                        showgrid = FALSE), yaxis = list(title = NULL,
                                                                                                        showticklabels = FALSE, zeroline = FALSE, showline = FALSE,
                                                                                                        showgrid = FALSE), zaxis = list(title = NULL,
                                                                                                                                        showticklabels = FALSE, zeroline = FALSE, showline = FALSE,
                                                                                                                                        showgrid = FALSE)))
    }
    else {
      copykat_plot <- copykat_plot %>% layout(scene = list(xaxis = list(title = paste0(toupper(dimensionality_reduction),
                                                                                       "_1")), yaxis = list(title = paste0(toupper(dimensionality_reduction),
                                                                                                                           "_2")), zaxis = list(title = paste0(toupper(dimensionality_reduction),
                                                                                                                                                               "_3"))))
    }
    plots[["CopyKat Plot"]] <- copykat_plot
    if (save_each_plot && !is.null(output_dir)) {
      copykat_plot <- copykat_plot %>% layout(title = "CopyKat (CNVs)",
                                              title = list(font = "black"), margin = list(t = 50))
      saveWidget(as_widget(copykat_plot), file = file.path(output_dir,
                                                           "Figures_Individual_Plots_HTML", "CopyKat_CNVs.html"),
                 selfcontained = FALSE, libdir = "lib")
    }
  }
  plots[["histograms"]] <- histograms
  return(plots)
}
