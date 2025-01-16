#' Plot SNV Data
#'
#' This function generates interactive 3D dimensionality reduction plots, histograms, Cell Types Plot, and CopyKat Plot, with the option to save plots individually.
#'
#' @importFrom Seurat Embeddings as.SingleCellExperiment GetAssayData
#' @importFrom copykat copykat CNA.MCMC annotateGenes.hg20 annotateGenes.mm10 baseline.GMM baseline.norm.cl convert.all.bins.hg20 heatmap.3
#' @importFrom ggplot2 ggplot ggsave aes geom_histogram xlab ylab theme element_text element_blank
#' @importFrom dplyr %>% filter group_by
#' @importFrom plotly ggplotly add_markers plot_ly add_trace layout as_widget
#' @importFrom htmlwidgets saveWidget
#' @importFrom Matrix sparseMatrix
#' @importFrom umap umap
#' @importFrom Rtsne Rtsne
#'
#' @param seurat_object Processed Seurat object with SNV metadata.
#' @param processed_snv Processed SNV data.
#' @param aggregated_snv Processed SNV data with additional computed quantities.
#' @param plot_data Data frame for 3D plotting.
#' @param output_dir Directory to save the plots (required if save_each_plot is TRUE).
#' @param include_histograms Logical; whether to generate histograms. Default: TRUE.
#' @param include_cell_types Logical; whether to generate cell types plot. Default: FALSE.
#' @param include_snv_dim_red Logical; whether to generate the SNV dimensionality reduction plot. Default: TRUE.
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
#' - **3D Plots**: Visualizes metrics such as SNV.N, SNVCount, and TotalVAF, etc.
#' - **Histograms**: Distribution of metrics including N_SNV, TotalVAF, MeanSNVsVAF, and N_VARreadCounts.
#' - **Cell Types Plot**: Shows cell types (e.g., custom classifications) with optional slingshot trajectories.
#' - **Transposed SNV Matrix Plot**: Makes a dimensionality reduction plot for SNVs instead of cells using transposed SNV matrix.
#' - **CopyKat Plot**: Depicts Copy Number Variations (CNVs) using CopyKat analysis.
#'
#' The output_dir parameter specifies where the plots will be saved if save_each_plot is enabled.
#'
#' @examples
#' # Example usage:
#' plots <- plot_snv_data(
#'   seurat_object = processed_data$SeuratObject,
#'   processed_snv = processed_data$ProcessedSNV,
#'   aggergated_snv = processed_data$AggregatedSNV
#'   plot_data = processed_data$PlotData,
#'   output_dir = "output/plots",
#'   include_histograms = TRUE,
#'   include_cell_types = TRUE,
#'   include_snv_dim_red = TRUE,
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
plot_snv_data <- function(seurat_object, processed_snv, aggregated_snv, plot_data, output_dir = NULL,
                          include_histograms = T, include_cell_types = F, include_snv_dim_red = T,
                          include_copykat = F, dimensionality_reduction = "UMAP", slingshot = T,
                          color_scale = "YlOrRd", cell_size = 3, disable_3d_axis = F, save_each_plot = F) {

  cat("\nGenerating SNV data plots...\n")

  valid_reductions <- c("umap", "pca", "tsne")
  dimensionality_reduction <- tolower(dimensionality_reduction)
  if (!dimensionality_reduction %in% valid_reductions) {
    stop("Invalid dimensionality_reduction method. Please use one of: 'umap', 'pca', 'tsne'.")
  }

  scale_list <- c("Blues", "Reds", "YlOrRd", "YlGnBu", "plasma", "RdBu")
  reversescale_options <- c(T, F, T, T, F, F)
  if (!color_scale %in% scale_list) {
    stop("Invalid color_scale. Please choose from: ", paste(scale_list, collapse = ", "))
  }

  reversescale_option <- reversescale_options[which(scale_list == color_scale)]
  histogram_scale1 <- "tomato"
  histogram_scale2 <- "dodgerblue2"
  color_undetected <- "lightgrey"

  # validate inputs i.e. processed SNV data
  if (is.null(processed_snv)) {
    stop("processed_snv must be provided for the SNV dimensionality reduction plot.")
  }

  required_cols <- c("CHROM", "POS", "REF", "ALT", "ReadGroup", "SNVCount", "RefCount")
  if (!all(required_cols %in% colnames(processed_snv))) {
    stop("processed_snv must contain the following columns: ", paste(required_cols, collapse = ", "))
  }


  plots <- list()  # list for appending all plots
  if (save_each_plot && !is.null(output_dir)) {
    dir.create(file.path(output_dir, "SNV_data_plots"),
               showWarnings = F, recursive = T)
  }

  curves <- NULL
  if (slingshot) {
    sce <- as.SingleCellExperiment(seurat_object)
    sce <- slingshot(sce, clusterLabels = "seurat_clusters",
                     reducedDim = toupper(dimensionality_reduction))
    curves <- slingCurves(sce, as.df = T)
    colnames(curves)[1:3] <- c(paste0(toupper(dimensionality_reduction), "_1"),
                               paste0(toupper(dimensionality_reduction), "_2"),
                               paste0(toupper(dimensionality_reduction), "_3"))
  }


  histograms <- list()
  if (include_histograms) {
    hist_list <- list(list(aes = aes(x = TotalVAF), xlab = "TotalVAF",
                           file_suffix = "Histogram_TotalVAF"),
                      list(aes = aes(x = MeanVAF), xlab = "MeanSNVsVAF",
                           file_suffix = "Histogram_MeanSNVsVAF"),
                      list(aes = aes(x = SNVCount), xlab = "N_VARreadCounts",
                           file_suffix = "Histogram_N_VARreadCounts"),
                      list(aes = aes(x = SNV.N), xlab = "N_SNV",
                           file_suffix = "Histogram_N_SNV"))

    for (hist_info in hist_list) {
      p <- ggplot(aggregated_snv, hist_info$aes) +
        geom_histogram(fill = histogram_scale2, boundary = 0, alpha = 0.6) +
                       xlab(hist_info$xlab) + ylab("Cells") +
        theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
              panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
              axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

      histograms[[hist_info$xlab]] <- ggplotly(p)
      if (save_each_plot && !is.null(output_dir)) {
        ggsave(file = file.path(output_dir, "SNV_data_plots",
                                paste0(hist_info$file_suffix, ".png")),
               plot = p, device = "png")
      }
    }
  }
  plots[["histograms"]] <- histograms

  plot_titles <- list(SNV.N = "N_sceSNVs",
                      SNVCount = "N_VARreads",
                      RefCount = "N_REFreads",
                      TotalVAF = "Total_VAF_RNA",
                      MeanVAF = "Mean_VAF_RNA",
                      MedianVAF = "Median_VAF_RNA")

  dimensionality_reduction <- tolower(dimensionality_reduction)
  dim_plotting <- toupper(dimensionality_reduction)


  # Function for generating all SNV related plots excluding individual SNV plots

  generate_plot <- function(metric, plot_data, dim_plotting, color_scale, reversescale_option,
                            cell_size, color_undetected, curves, disable_3d_axis, title) {
    plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
      add_markers(
        data = plot_data[plot_data[["Undetected"]] == 0, ],
        x = ~get(paste0(dim_plotting, "_1")),
        y = ~get(paste0(dim_plotting, "_2")),
        z = ~get(paste0(dim_plotting, "_3")),
        marker = list(
          color = ~get(metric), colorscale = color_scale,
          reversescale = reversescale_option, size = cell_size,
          showscale = T, colorbar = list(x = 0.85, y = 0.5, thickness = 20, len = 0.5),
          line = list(color = ~get(metric), width = cell_size)),
                      opacity = 0.8, name = "expressed sceSNV loci"
      ) %>%
      add_markers(
        data = plot_data[plot_data[["Undetected"]] == 1, ],
        x = ~get(paste0(dim_plotting, "_1")),
        y = ~get(paste0(dim_plotting, "_2")),
        z = ~get(paste0(dim_plotting, "_3")),
        marker = list(color = color_undetected, size = cell_size,
                      opacity = 0.5), name = "undetected"
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
        title = title, title = list(font = "black"), margin = list(t = 50),
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
            xaxis = list(title = NULL, showticklabels = F),
            yaxis = list(title = NULL, showticklabels = F),
            zaxis = list(title = NULL, showticklabels = F)
          )
        )
    }

    return(plot)
  }


  #plots <- list()

  for (metric in names(plot_titles)) {
    if (!metric %in% colnames(plot_data))
      next

    plot <- generate_plot(
      metric = metric,
      plot_data = plot_data,
      dim_plotting = dim_plotting,
      color_scale = color_scale,
      reversescale_option = reversescale_option,
      cell_size = cell_size,
      color_undetected = color_undetected,
      curves = curves,
      disable_3d_axis = disable_3d_axis,
      title = plot_titles[[metric]]
    )

    plots[[metric]] <- plot

    file_title <- gsub("[^a-zA-Z0-9]", "_", plot_titles[[metric]])

    if (save_each_plot && !is.null(output_dir)) {
      saveWidget(as_widget(plot), file = file.path(
        output_dir, "SNV_data_plots", paste0(file_title, ".html")),
        selfcontained = F, libdir = "lib")
    }
  }


  # CELL TYPES PLOT

  if (include_cell_types && "customclassif" %in% colnames(plot_data)) {
    cell_type_plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
      add_trace(data = plot_data, x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
                y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
                z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
                size = 0.05, opacity = 0.5, mode = "markers", color = ~customclassif)
    if (slingshot && !is.null(curves)) {
      cell_type_plot <- cell_type_plot %>%
        add_trace(
               data = curves,
               x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
               y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
               z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
               mode = "lines", color = ~factor(Lineage), line = list(width = 4)
            )
    }
    if (disable_3d_axis) {
      cell_type_plot <- cell_type_plot %>%
        layout(scene = list(xaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                         showline = F, showgrid = F),
                            yaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                         showline = F, showgrid = F),
                            zaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                         showline = F, showgrid = F)))
    }
    else {
      cell_type_plot <- cell_type_plot %>% layout(
        title = '',
        scene = list(xaxis = list(title = paste0(toupper(dimensionality_reduction), "_1")),
                     yaxis = list(title = paste0(toupper(dimensionality_reduction), "_2")),
                     zaxis = list(title = paste0(toupper(dimensionality_reduction), "_3"))))
    }

    plots[["Cell Types Plot"]] <- cell_type_plot

    if (save_each_plot) {
      cell_type_plot <- cell_type_plot %>% layout(
        title = "Cell types (scType)", title = list(font = "black"), margin = list(t = 50))
      saveWidget(as_widget(cell_type_plot), file = file.path(
        output_dir, "SNV_data_plots", "Cell_types_scType.html"),
        selfcontained = F, libdir = "lib")
    }

    if (slingshot) {
      s <- slingshot(as.SingleCellExperiment(seurat_object),
                     clusterLabels = "seurat_clusters",
                     reducedDim = toupper(dimensionality_reduction))
      lineage_counter <- 1
      for (i in grep("slingPseudotime", colnames(s@colData))) {
        lineage_id <- as.integer(sub("slingPseudotime_", "", colnames(s@colData)[i]))
        curve <- curves[curves$Lineage == lineage_id, ]
        lineage_plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
          add_trace(data = plot_data,
                    x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
                    y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
                    z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
                    marker = list(color = ~s@colData[, i], colorscale = "YlOrRd",
                                  reversescale = reversescale_option, colorbar = list(x = 0),
                                  line = list(color = ~SNV.N, width = cell_size)))
        if (!disable_3d_axis) {
          lineage_plot <- lineage_plot %>% layout(
            scene = list(xaxis = list(title = paste0(toupper(dimensionality_reduction), "_1")),
                         yaxis = list(title = paste0(toupper(dimensionality_reduction), "_2")),
                         zaxis = list(title = paste0(toupper(dimensionality_reduction), "_3")))) %>%
            add_trace(data = curve,
                      x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
                      y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
                      z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
                      mode = "lines", showlegend = T, line = list(width = 4, color = "black"))
        }
        else {
          lineage_plot <- lineage_plot %>% layout(
            scene = list(xaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                      showline = F, showgrid = F),
                         yaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                      showline = F, showgrid = F),
                         zaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                      showline = F, showgrid = F)))
        }
      }
      lineage_counter <- lineage_counter + 1
    }
  }

  # SNV dimensionality reduction plot integration
  # Plot SNVs - This is a transposed SNV-Cellbarcode matrix

  if (include_snv_dim_red) {
    processed_snv$SNV <- paste0(processed_snv$CHROM, "_", processed_snv$POS, "_",
                                processed_snv$REF, ">", processed_snv$ALT)
    unique_snvs <- unique(processed_snv$SNV)
    unique_readgroups <- unique(processed_snv$ReadGroup)

    if (length(unique_snvs) == 0 || length(unique_readgroups) == 0) {
      stop("processed_snv contains no valid SNVs or ReadGroups.")
    }

    processed_snv$snv_idx <- match(processed_snv$SNV, unique_snvs)
    processed_snv$readgroups_idx <- match(processed_snv$ReadGroup, unique_readgroups)

    snv_mat_VAR <- sparseMatrix(
      i = processed_snv$snv_idx,
      j = processed_snv$readgroups_idx,
      x = processed_snv$SNVCount
    )
    rownames(snv_mat_VAR) <- unique_snvs
    colnames(snv_mat_VAR) <- unique_readgroups

    snv_mat_REF <- sparseMatrix(
      i = processed_snv$snv_idx,
      j = processed_snv$readgroups_idx,
      x = processed_snv$RefCount
    )
    rownames(snv_mat_REF) <- unique_snvs
    colnames(snv_mat_REF) <- unique_readgroups

    snv_mat <- cbind(as.matrix(snv_mat_VAR), as.matrix(snv_mat_REF))
    snv_mat <- scale(snv_mat[, colSums(snv_mat) > 0])  # Filter and scale

    if (dimensionality_reduction == "tsne") {
      snv_mat_reduced <- Rtsne(snv_mat, dims = 3, perplexity = max(5, (nrow(snv_mat) - 2) / 3))$Y
    } else if (dimensionality_reduction == "pca") {
      snv_mat_reduced <- prcomp(snv_mat, rank. = 3)$x
    } else {
      n_neighbors <- ifelse(nrow(snv_mat) > 30, 5, 3)
      snv_mat_reduced <- umap(snv_mat, n_neighbors = n_neighbors, n_components = 3)$layout
    }

    df_3dplot_snv <- as.data.frame(snv_mat_reduced)
    colnames(df_3dplot_snv) <- c(
      paste0(toupper(dimensionality_reduction), "_1"),
      paste0(toupper(dimensionality_reduction), "_2"),
      paste0(toupper(dimensionality_reduction), "_3")
    )

    trans_snv_plot <- plot_ly(
      type = "scatter3d", mode = "lines+markers"
    ) %>% add_markers(
            data = df_3dplot_snv,
            x = ~get(colnames(df_3dplot_snv)[1]),
            y = ~get(colnames(df_3dplot_snv)[2]),
            z = ~get(colnames(df_3dplot_snv)[3]),
            marker = list(size = cell_size, color = color_scale, opacity = 0.5),
            text = rownames(df_3dplot_snv), hoverinfo = "text"
      )

    if (disable_3d_axis) {
      trans_snv_plot <- trans_snv_plot %>% layout(
        scene = list(
          xaxis = list(title = NULL, showticklabels = F),
          yaxis = list(title = NULL, showticklabels = F),
          zaxis = list(title = NULL, showticklabels = F)
        )
      )
    } else {
      trans_snv_plot <- trans_snv_plot %>% layout(
        title = '',
        scene = list(
          xaxis = list(title = colnames(df_3dplot_snv)[1]),
          yaxis = list(title = colnames(df_3dplot_snv)[2]),
          zaxis = list(title = colnames(df_3dplot_snv)[3])
        )
      )
    }

    if (save_each_plot && !is.null(output_dir)) {
      trans_snv_plot <- trans_snv_plot %>%
        layout(title = 'SNVs', title = list(font = 'black'), margin = list(t = 50))
      saveWidget(as_widget(trans_snv_plot), file = file.path(
        output_dir, "SNV_data_plots", "Transposed_SNV_Matrix.html"),
        selfcontained = F, libdir = "lib"
      )
    }

    plots[["Transposed SNV Matrix Plot"]] <- trans_snv_plot
  }


  # COPYKAT PLOT

  if (include_copykat) {
    cat("Running CopyKat. This may take a while...\n")
    cts <- GetAssayData(seurat_object, layer = "counts", assay = "SCT")
    ckt <- copykat(cts, sam.name = "sample_name")
    seurat_object[["karyotype"]] <- "Unknown"
    seurat_object[["karyotype"]][rownames(ckt$pred), ] <- ckt$pred[, 2]
    df.srt <- as.data.frame(seurat_object@meta.data[, c("seurat_clusters", "karyotype")])
    df.embed <- as.data.frame(Embeddings(seurat_object, reduction = dimensionality_reduction))
    plot_data <- merge(df.srt, df.embed, by = 0)
    colnames(plot_data) <- c("filler", "seurat_clusters", "karyotype",
                             paste0(toupper(dimensionality_reduction), "_1"),
                             paste0(toupper(dimensionality_reduction), "_2"),
                             paste0(toupper(dimensionality_reduction), "_3"))
    rownames(plot_data) <- plot_data$Row.names
    plot_data <- plot_data[, 2:ncol(plot_data)]

    copykat_plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
      add_markers(data = plot_data,
                  x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
                  y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
                  z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
                  size = 0.05, opacity = 1, color = ~karyotype)

    if (slingshot && !is.null(curves)) {
      copykat_plot <- copykat_plot %>%
        add_trace(data = curves,
                  x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
                  y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
                  z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
                  split = ~Lineage, mode = "lines")
    }
    if (disable_3d_axis) {
      copykat_plot <- copykat_plot %>%
        layout(scene = list(xaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                         showline = F, showgrid = F),
                            yaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                         showline = F, showgrid = F),
                            zaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                         showline = F, showgrid = F)))
    }
    else {
      copykat_plot <- copykat_plot %>% layout(
        title = '', scene = list(xaxis = list(title = paste0(toupper(dimensionality_reduction), "_1")),
                                 yaxis = list(title = paste0(toupper(dimensionality_reduction), "_2")),
                                 zaxis = list(title = paste0(toupper(dimensionality_reduction), "_3"))))
    }
    plots[["CopyKat Plot"]] <- copykat_plot

    if (save_each_plot) {
      copykat_plot <- copykat_plot %>%
        layout(title = "CopyKat (CNVs)", title = list(font = "black"), margin = list(t = 50))
      saveWidget(as_widget(copykat_plot), file = file.path(
        output_dir, "SNV_data_plots", "CopyKat_CNVs.html"),
        selfcontained = F, libdir = "lib")
    }
  }

  return(plots)
}
