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
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom slingshot slingshot slingCurves
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
#' @param cell_border Numeric; thickness of cells'/markers' border in the plot. Default: 0.
#' @param disable_3d_axis Logical; whether to disable 3D axis labels. Default: FALSE.
#' @param save_each_plot Logical; whether to save each plot individually. Default: FALSE.
#' @param enable_integrated Logical; whether to use an integrated Seurat object. Default: FALSE.
#' @return A list of generated plots.
#' @details
#' The function generates various visualizations for SNV data in the 'plots' object, including:
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
#' plots <- plot_snv_data(seurat_object = processed_data$SeuratObject,
#'                        processed_snv = processed_data$ProcessedSNV,
#'                        aggregated_snv = processed_data$AggregatedSNV,
#'                        plot_data = processed_data$PlotData,
#'                        output_dir = output_dir,
#'                        include_histograms = TRUE,  
#'                        dimensionality_reduction = "umap",
#'                        include_cell_types = TRUE,
#'                        include_copykat = FALSE, 
#'                        include_snv_dim_red = FALSE,
#'                        slingshot = TRUE,
#'                        color_scale = "YlOrRd",
#'                        cell_border = 0,
#'                        save_each_plot = TRUE)
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
                          color_scale = "YlOrRd", cell_border = 0, disable_3d_axis = F, save_each_plot = F, 
                          enable_integrated = F) {

  cat("\nGenerating SNV data plots...\n")

  if (enable_integrated){
    slingshot=FALSE
    include_copykat=FALSE
    if (slingshot==TRUE){print('Turning off the slingshot option - this is not available for integrated objects')}
    if (include_copykat==TRUE){print('Turning off the copykat option - this is not available for integrated objects')}
  }

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
    sce <- as.SingleCellExperiment(seurat_object, assay='SCT')
    sce <- slingshot(sce, clusterLabels = "seurat_clusters",
                     reducedDim = toupper(dimensionality_reduction))
    curves <- slingCurves(sce, as.df = T)
    colnames(curves)[1:3] <- c(paste0(toupper(dimensionality_reduction), "_1"),
                               paste0(toupper(dimensionality_reduction), "_2"),
                               paste0(toupper(dimensionality_reduction), "_3"))
  }


  histograms <- list()
  if (include_histograms) {
    options(bitmapType = "cairo-png")
    hist_list <- list(list(aes = aes(x = TotalVAF), xlab = "TotalVAF",
                           file_suffix = "Histogram_TotalVAF", binwidth = 0.05),
                      list(aes = aes(x = MeanVAF), xlab = "MeanSNVsVAF",
                           file_suffix = "Histogram_MeanSNVsVAF", binwidth = 0.05),
                      list(aes = aes(x = SNVCount), xlab = "N_VARreadCounts",
                           file_suffix = "Histogram_N_VARreadCounts", binwidth = 10),
                      list(aes = aes(x = SNV.N), xlab = "N_SNV",
                           file_suffix = "Histogram_N_SNV", binwidth = 1))

    for (hist_info in hist_list) {
    p <- ggplot(aggregated_snv, hist_info$aes)

    if (enable_integrated) {    
      n_colors <- length(unique(seurat_object$orig.ident))
      palette <- distinctColorPalette(n_colors)  

      for (i in seq_along(unique(seurat_object$orig.ident))) {
        p <- p + geom_histogram(
          data = aggregated_snv[aggregated_snv$orig.ident == unique(seurat_object$orig.ident)[i], ],
          aes(y = after_stat(count), fill = factor(orig.ident)),
          boundary = 0, alpha = 0.3, binwidth = hist_info$binwidth, show.legend = TRUE
        )
      }

      p <- p +
        scale_fill_manual(values = palette, name = "Sample") +  
        xlab(hist_info$xlab) + ylab("Cells") +
        theme_minimal()

      histograms[[hist_info$xlab]] <- ggplotly(p)

    } else {
      p <- p + geom_histogram(
        aes(y = after_stat(count)),
        fill = histogram_scale2, boundary = 0, alpha = 0.6, binwidth = hist_info$binwidth
      ) +
        xlab(hist_info$xlab) + ylab("Cells") +
        theme_minimal()

      histograms[[hist_info$xlab]] <- ggplotly(p)
    }

    if (save_each_plot && !is.null(output_dir)) {
      ggsave(p, file = file.path(output_dir, "SNV_data_plots",
                                paste0(hist_info$file_suffix, ".png")),
            device = "png")
    }
  }
    options(bitmapType = "C_X11")
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
                            cell_border, color_undetected, curves, disable_3d_axis, title) {
    plot <- plot_ly(type = "scatter3d", mode = "lines+markers") 
    max_metric_val = max(plot_data[plot_data[["Undetected"]] == 0,metric])
    for (i in 1:length(unique(seurat_object$orig.ident))){
      this.id = unique(seurat_object$orig.ident)[i]
      if (any(plot_data[plot_data[["Undetected"]] == 0 & plot_data[["orig.ident"]] == this.id, metric]==max_metric_val)){
      plot <- plot %>% add_markers(
        data = plot_data[plot_data[["Undetected"]] == 0 & plot_data[["orig.ident"]] == this.id, ],
        x = ~get(paste0(dim_plotting, "_1")),
        y = ~get(paste0(dim_plotting, "_2")),
        z = ~get(paste0(dim_plotting, "_3")),
        size = 0.05, opacity = 1.00,
        marker = list(
          color = ~get(metric), colorscale = color_scale,
          reversescale = reversescale_option, showscale = T,
          colorbar = list(x = 0.85, y = 0.5, thickness = 20, len = 0.5),
          line = list(color = ~get(metric), width = cell_border)),
          name = paste0(this.id," expressed sceSNV loci"), hovertext=plot_data["orig.ident"]
      )} else {
             plot <- plot %>% add_markers(
        data = plot_data[plot_data[["Undetected"]] == 0 & plot_data[["orig.ident"]] == this.id, ],
        x = ~get(paste0(dim_plotting, "_1")),
        y = ~get(paste0(dim_plotting, "_2")),
        z = ~get(paste0(dim_plotting, "_3")),
        size = 0.05, opacity = 1.00,
        marker = list(
          color = ~get(metric), colorscale = color_scale,
          reversescale = reversescale_option,
          line = list(color = ~get(metric), width = cell_border)),
          name = paste0(this.id," expressed sceSNV loci"), hovertext=plot_data["orig.ident"])
      }
      plot <- plot %>% add_markers(
        data = plot_data[plot_data[["Undetected"]] == 1 & plot_data[["orig.ident"]] == this.id, ],
        x = ~get(paste0(dim_plotting, "_1")),
        y = ~get(paste0(dim_plotting, "_2")),
        z = ~get(paste0(dim_plotting, "_3")),
        size = 0.05, opacity = 1.00,
        marker = list(color = color_undetected,
                      line = list(width = cell_border)),
        name = paste0(this.id," undetected"), hovertext=plot_data["orig.ident"]
      )
      }

    if (!is.null(curves)) {
      plot <- plot %>%
        add_trace(
          data = curves,
          x = ~get(paste0(dim_plotting, "_1")),
          y = ~get(paste0(dim_plotting, "_2")),
          z = ~get(paste0(dim_plotting, "_3")),
          split = ~Lineage, mode = "lines"
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
    plot <-plot %>% layout(legend= list(itemsizing='constant'))
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
      cell_border = cell_border,
      color_undetected = color_undetected,
      curves = curves,
      disable_3d_axis = disable_3d_axis,
      title = plot_titles[[metric]] #(duplicate title fix)
    )

    plot_name <- plot_titles[[metric]]
    plots[[plot_name]] <- plot
    #plots[[metric]] <- plot (changed naming of the plots in the output for easier matching in generate_report)

    file_title <- gsub("[^a-zA-Z0-9]", "_", plot_titles[[metric]])

    if (save_each_plot && !is.null(output_dir)) {
      saveWidget(as_widget(plot), file = file.path(
        output_dir, "SNV_data_plots", paste0(file_title, ".html")),
        selfcontained = F, libdir = "lib")
    }
  }

  # INTEGRATED ORIG IDENT PLOT

  if (enable_integrated){
    origident_plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
    add_trace(data = plot_data, x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
                y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
                z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
                size = 0.05, opacity = 0.5, mode = "markers", color = ~orig.ident)
    if (disable_3d_axis) {
      origident_plot <- origident_plot %>%
      layout(scene = list(xaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                         showline = F, showgrid = F),
                            yaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                         showline = F, showgrid = F),
                            zaxis = list(title = NULL, showticklabels = F, zeroline = F,
                                         showline = F, showgrid = F)))}
    else {
      origident_plot <- origident_plot %>% layout(title = '',
      scene = list(xaxis = list(title = paste0(toupper(dimensionality_reduction), "_1")),
      yaxis = list(title = paste0(toupper(dimensionality_reduction), "_2")),
      zaxis = list(title = paste0(toupper(dimensionality_reduction), "_3"))))}
    
    plots[["Sample ID"]] <- origident_plot
    
    if (save_each_plot && !is.null(output_dir)) {
      saveWidget(as_widget(origident_plot), file = file.path(
        output_dir, "SNV_data_plots", "SAMPLE_ID_plot.html"),
        selfcontained = F, libdir = "lib")
    }
    }

  # CELL TYPES PLOT

  if (include_cell_types && "customclassif" %in% colnames(plot_data)) {
    if (enable_integrated){
      plot_data['origident_customclassif'] = paste0(plot_data$orig.ident,' ',plot_data$customclassif)
      cell_type_plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
        add_trace(data = plot_data, x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
                  y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
                  z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
                  size = 0.05, opacity = 0.5, mode = "markers", color = ~origident_customclassif)
    } else {
      cell_type_plot <- plot_ly(type = "scatter3d", mode = "lines+markers") %>%
      add_trace(data = plot_data, x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
                y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
                z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
                size = 0.05, opacity = 0.5, mode = "markers", color = ~customclassif)
    }
    if (slingshot && !is.null(curves)) {
      cell_type_plot <- cell_type_plot %>%
        add_trace(
          data = curves,
          x = ~get(paste0(toupper(dimensionality_reduction), "_1")),
          y = ~get(paste0(toupper(dimensionality_reduction), "_2")),
          z = ~get(paste0(toupper(dimensionality_reduction), "_3")),
          mode = "lines", color = ~factor(Lineage)
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
    plots[["Cell types (scType)"]] <- cell_type_plot

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
                                  line = list(color = ~SNV.N, width = cell_border)))
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
  if (include_snv_dim_red){
    processed_snv$SNV <- paste0(processed_snv$CHROM, "_", processed_snv$POS, "_",
                                processed_snv$REF, ">", processed_snv$ALT)
    unique_snvs <- unique(processed_snv$SNV)
    if (length(unique_snvs)<100){
      cat("You need at least 100 unique SNVs to create a transposed SNV plot. Moving past that part now...\n")
      include_snv_dim_red = FALSE
    }
  }

  if (include_snv_dim_red) {
    processed_snv$SNV <- paste0(processed_snv$CHROM, "_", processed_snv$POS, "_",
                                processed_snv$REF, ">", processed_snv$ALT)
    unique_snvs <- unique(processed_snv$SNV)
    unique_readgroups <- unique(processed_snv$ReadGroup)

    if (length(unique_snvs) == 0 || length(unique_readgroups) == 0) {
      stop("processed_snv contains no valid SNVs or ReadGroups.")
    }
    
    processed_snv_flt = processed_snv[is.numeric(processed_snv$VAF)==1 & is.finite(processed_snv$VAF)==TRUE,]
    cell_ids = data.frame(ReadGroup = unique(processed_snv_flt$ReadGroup),cell_n=(1:length(unique(processed_snv_flt$ReadGroup))))
    snv_ids = data.frame(SNV = unique(processed_snv_flt$SNV),snv_val=(1:length(unique(processed_snv_flt$SNV))))
    df_transpose = merge(processed_snv_flt,cell_ids, on='ReadGroup',how='left')
    df_transpose = merge(df_transpose,snv_ids, on='ReadGroup',how='left')
    
    mtx = Matrix(0, nrow=length(unique(df_transpose$ReadGroup)), ncol=length(unique(df_transpose$SNV)))
    
    mtx[cbind(df_transpose$cell_n,df_transpose$snv_val)] <- df_transpose$VAF
    colnames(mtx) <- snv_ids$SNV
    rownames(mtx) <- cell_ids$ReadGroup
    srt_transposed <- CreateSeuratObject(mtx,min.features=1)
    srt_transposed <- NormalizeData(srt_transposed, verbose=F)
    srt_transposed <- FindVariableFeatures(srt_transposed, verbose=F)
    srt_transposed <- ScaleData(srt_transposed, verbose=F)
    
    if (dimensionality_reduction == "tsne") {
      srt_transposed = RunTSNE(srt_transposed, dim.embed = 3, verbose=F)
      snv_mat_reduced <- as.data.frame(Embeddings(srt_transposed, reduction = dimensionality_reduction))
    } else if (dimensionality_reduction == "pca") {
      srt_transposed = RunPCA(srt_transposed, dim.embed = 3, verbose=F)
      snv_mat_reduced <- as.data.frame(Embeddings(srt_transposed, reduction = dimensionality_reduction))
    } else {
      srt_transposed = RunPCA(srt_transposed, verbose=F)
      srt_transposed <- FindNeighbors(srt_transposed, verbose=F)
      srt_transposed <- FindClusters(srt_transposed, verbose=F)
      srt_transposed = RunUMAP(srt_transposed, n.components = 3, dims=1:20, verbose=F)
      snv_mat_reduced <- as.data.frame(Embeddings(srt_transposed, reduction = dimensionality_reduction))
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
      marker = list(size = cell_border, color = color_scale, opacity = 0.5),
      text = rownames(df_3dplot_snv), hoverinfo = "text"
    )


    if (disable_3d_axis) {
      trans_snv_plot <- trans_snv_plot %>% layout(
        showlegend = FALSE,
        scene = list(
          xaxis = list(title = NULL, showticklabels = F),
          yaxis = list(title = NULL, showticklabels = F),
          zaxis = list(title = NULL, showticklabels = F)
        )
      )
    } else {
      trans_snv_plot <- trans_snv_plot %>% layout(
        title = '',
        showlegend = FALSE,
        scene = list(
          xaxis = list(title = colnames(df_3dplot_snv)[1]),
          yaxis = list(title = colnames(df_3dplot_snv)[2]),
          zaxis = list(title = colnames(df_3dplot_snv)[3])
        )
      )
    }

    if (save_each_plot && !is.null(output_dir)) {
      trans_snv_plot <- trans_snv_plot %>%
        layout(title = "Transposed SNV matrix", title = list(font = 'black'), margin = list(t = 50))
      saveWidget(as_widget(trans_snv_plot), file = file.path(
        output_dir, "SNV_data_plots", "Transposed_SNV_Matrix.html"),
        selfcontained = F, libdir = "lib"
      )
    }

    plots[["Transposed SNV matrix"]] <- trans_snv_plot
  }


  # COPYKAT PLOT

  if (include_copykat) {
    options(bitmapType = "cairo")
    cat("Running CopyKat. This may take a while...\n")
    cts <- GetAssayData(seurat_object, layer = "counts", assay = "SCT")
    ckt <- copykat(cts, sam.name = "sample_name")
    options(bitmapType = "C_X11")
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
    plots[["CNVs (CopyKat)"]] <- copykat_plot

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
