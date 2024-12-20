#' Individual SNV Plots
#'
#' This function generates individual SNV plots for VAF, N_VAR, and N_REF. It uses processed SNV data
#' and a Seurat object to create plots with options for saving and including slingshot trajectories.
#'
#' @importFrom Seurat Embeddings as.SingleCellExperiment
#' @importFrom slingshot slingshot slingCurves
#' @importFrom SingleCellExperiment reducedDims reducedDims<-
#' @importFrom plotly plotly_json
#'
#' @param seurat_object Processed Seurat object.
#' @param processed_snv Data frame of processed SNV information, typically output from the preprocess_snv_data function.
#' @param output_dir Directory to save plots and HTML (if save_each_plot is TRUE).
#' @param slingshot Logical; whether to include slingshot trajectories. Default: TRUE.
#' @param dimensionality_reduction Dimensionality reduction method ('UMAP', 'PCA', 'tSNE'). Default: "UMAP".
#' @param dynamic_cell_size Logical; whether to scale cell size dynamically based on SNV and reference read counts. Default: FALSE.
#' @param save_each_plot Logical; whether to save each plot individually. Default: FALSE.
#' @return A list containing JSON content for VAF, N_VAR, and N_REF plots.
#' @details
#' This function generates individual SNV plots using processed SNV data (processed_snv) and the dimensional
#' reduction embeddings from a Seurat object. The plots visualize key metrics such as VAF (Variant Allele Fraction),
#' N_VAR (number of variant reads), and N_REF (number of reference reads).
#'
#' The plots can be saved individually in the specified output_dir if save_each_plot is set to TRUE.
#'
#' @examples
#' # Example usage:
#' snv_plots <- individual_snv_plots(
#'   seurat_object = processed_data$SeuratObject,
#'   processed_snv = processed_data$ProcessedSNV,
#'   output_dir = "output/individual_plots",
#'   slingshot = TRUE,
#'   dimensionality_reduction = "UMAP",
#'   dynamic_cell_size = FALSE,
#'   save_each_plot = TRUE
#' )
#' @export
#'
individual_snv_plots <- function(seurat_object, processed_snv, output_dir = NULL, slingshot = TRUE,
                                 dimensionality_reduction = "UMAP", dynamic_cell_size = FALSE,
                                 save_each_plot = FALSE) {


  valid_reductions <- c("umap", "pca", "tsne")
  dimensionality_reduction <- tolower(dimensionality_reduction)
  if (!dimensionality_reduction %in% valid_reductions) {
    stop("Invalid dimensionality_reduction method. Please use one of: 'umap', 'pca', 'tsne'.")
  }
  dim.title <- switch(dimensionality_reduction, umap = "UMAP",
                      pca = "PCA", tsne = "tSNE")
  pal <- c("#EBEBEB", "#85C1E9", "#E74C3C", "#B03A2E", "#641E16")
  df.dim <- as.data.frame(Embeddings(seurat_object, reduction = dimensionality_reduction))
  colnames(df.dim) <- c("x", "y", "z")
  df.snv <- processed_snv
  df.snv <- df.snv[c("CHROM", "POS", "REF", "ALT", "ReadGroup",
                     "SNVCount", "RefCount", "VAF")]
  snvs <- unique(df.snv[c("CHROM", "POS", "REF", "ALT")])
  snv_options <- paste(snvs$CHROM, snvs$POS, snvs$REF, snvs$ALT,
                       sep = ":")

  individual_SNV_html <- NULL
  curves <- NULL

  if (slingshot) {
    sce <- as.SingleCellExperiment(seurat_object)
    if (!dimensionality_reduction %in% names(reducedDims(sce))) {
      reducedDims(sce)[[dimensionality_reduction]] <- Embeddings(seurat_object,
                                                                 reduction = tolower(dimensionality_reduction))
    }
    sce <- slingshot(sce, clusterLabels = "seurat_clusters",
                     reducedDim = dimensionality_reduction)
    curves <- slingCurves(sce, as.df = TRUE)
    colnames(curves)[1:3] <- c(paste0(dim.title, "_1"), paste0(dim.title, "_2"), paste0(dim.title, "_3"))
  }


  generate_snv_plots <- function(selected_snv, title_color = "blue",
                                 dynamic_cell_size = FALSE) {
    selected_parts <- unlist(strsplit(selected_snv, ":"))
    df_subset <- df.snv[df.snv$CHROM == selected_parts[1] &
                          df.snv$POS == as.numeric(selected_parts[2]) & 
                          df.snv$REF ==selected_parts[3] & 
                          df.snv$ALT == selected_parts[4], ]
                          
    vaf <- df_subset$VAF[match(colnames(seurat_object), df_subset$ReadGroup)]
    snv_reads <- df_subset$SNVCount[match(colnames(seurat_object),
                                          df_subset$ReadGroup)]
    ref_reads <- df_subset$RefCount[match(colnames(seurat_object),
                                          df_subset$ReadGroup)]
    y <- data.frame(x = df.dim[, 1], y = df.dim[, 2], z = df.dim[, 3],
                    vaf = vaf, ref_reads = ref_reads, snv_reads = snv_reads)

    plots <- list()
    y$vaf_label <- "Undetected"
    y$vaf_label[y$vaf == 0] <- "0 VAF, N_REF Only"
    y$vaf_label[0 < y$vaf & y$vaf <= 0.25] <- "0<VAF<=0.25"
    y$vaf_label[0.25 < y$vaf & y$vaf <= 0.75] <- "0.25<VAF<=0.75"
    y$vaf_label[0.75 < y$vaf & y$vaf <= 1] <- "0.75<VAF<=1.00"
    y$vaf_label <- factor(y$vaf_label, levels = c("Undetected", "0 VAF, N_REF Only",
                                                  "0<VAF<=0.25", "0.25<VAF<=0.75", "0.75<VAF<=1.00"))

    # make VAF plots
    f_vaf <- plot_ly(type = "scatter3d", mode = "markers+lines")
    for (j in 1:5) {
      cur_label <- levels(y$vaf_label)[j]
      if (dynamic_cell_size) {
        f_vaf <- f_vaf %>% add_trace(data = subset(y, vaf_label == cur_label),
                                     x = ~x, y = ~y, z = ~z,
                                     size = ~((snv_reads + ref_reads) / max(c(snv_reads, ref_reads), na.rm = TRUE)) * 10, type = "scatter3d",
                                     mode = "markers", marker = list(color = pal[j], line = list(width = 0)), name = cur_label)
      }
      else {
        f_vaf <- f_vaf %>% add_trace(data = subset(y, vaf_label == cur_label), x = ~x, y = ~y, z = ~z,
                                     size = 0.05, type = "scatter3d", mode = "markers",
                                     marker = list(color = pal[j], line = list(width = 0)),
                                     name = cur_label)
      }
    }
    if (!is.null(curves)) {
      f_vaf <- f_vaf %>% add_trace(data = curves,
                                   x = ~get(paste0(dim.title, "_1")),
                                   y = ~get(paste0(dim.title, "_2")),
                                   z = ~get(paste0(dim.title, "_3")),
                                   split = ~Lineage, mode = "lines", line = list(width = 2))
    }
    f_vaf <- f_vaf %>% layout(title = list(text = "VAF_RNA", font = list(color = title_color)),
                              scene = list(xaxis = list(title = paste0(dim.title, "_1")),
                                           yaxis = list(title = paste0(dim.title, "_2")),
                                           zaxis = list(title = paste0(dim.title, "_3"))))

    plots[["VAF"]] <- f_vaf

    # make N_VAR plots
    f_varreads <- plot_ly(type = "scatter3d", mode = "markers+lines") %>%
      add_trace(data = subset(y, !is.na(vaf)),
                x = ~x, y = ~y, z = ~z,
                size = ifelse(dynamic_cell_size,
                              ~((snv_reads + ref_reads) / max(c(snv_reads, ref_reads), na.rm = TRUE)) * 10, 0.05), type = "scatter3d",
                mode = "markers", marker = list(reversescale = TRUE, color = ~snv_reads, colorscale = "YlOrRd",
                                                showscale = TRUE, opacity = 0.5, line = list(color = "#FEE5D9", width = 1), colorbar = list(len = 0.5, y = 0.2)),
                name = "Cells with N_VAR") %>% add_trace(data = subset(y, is.na(vaf)),
                                                         x = ~x, y = ~y, z = ~z, size = ifelse(dynamic_cell_size,
                                                                                               ~((snv_reads + ref_reads) / max(c(snv_reads, ref_reads),
                                                                                               na.rm = TRUE)) * 10, 0.05), type = "scatter3d",
                                                         mode = "markers", marker = list(color = "#EBEBEB",
                                                                                         line = list(width = 0), opacity = 0.5), name = "Cells without N_VAR")
    if (!is.null(curves)) {
      f_varreads <- f_varreads %>% add_trace(data = curves,
                                             x = ~get(paste0(dim.title, "_1")),
                                             y = ~get(paste0(dim.title, "_2")),
                                             z = ~get(paste0(dim.title, "_3")),
                                             split = ~Lineage, mode = "lines", line = list(width = 2))
    }
    f_varreads <- f_varreads %>% layout(title = list(text = "N_VAR", font = list(color = title_color)),
                                        scene = list(xaxis = list(title = paste0(dim.title, "_1")),
                                                     yaxis = list(title = paste0(dim.title, "_2")),
                                                     zaxis = list(title = paste0(dim.title, "_3"))))

    plots[["N_VAR"]] <- f_varreads

    # make N_REF plots
    f_refreads <- plot_ly(type = "scatter3d", mode = "markers+lines") %>%
      add_trace(data = subset(y, vaf == 0 & ref_reads > 0),
                x = ~x, y = ~y, z = ~z, size = ifelse(dynamic_cell_size,
                                                      ~(snv_reads + ref_reads) / max(c(snv_reads, ref_reads),
                                                        na.rm = TRUE) * 10, 0.05), type = "scatter3d",
                                                mode = "markers", marker = list(reversescale = TRUE,
                                                color = ~ref_reads, colorscale = "Blues", showscale = TRUE,
                                                opacity = 0.5, line = list(color = "#EFF3FF", width = 1),
                                                colorbar = list(len = 0.5, y = 0.2)), name = "Cells with N_REF") %>% add_trace(
                                                         data = subset(y, (vaf == 0 & ref_reads == 0) | (is.na(vaf) == 1) | (vaf > 0)),
                                                         x = ~x, y = ~y, z = ~z,
                                                         size = ifelse(dynamic_cell_size,
                                                                       ~(snv_reads + ref_reads) / max(c(snv_reads, ref_reads),
                                                                        na.rm = TRUE) * 10, 0.05), type = "scatter3d",
                                                         mode = "markers", marker = list(color = "#EBEBEB",
                                                         line = list(width = 0), opacity = 0.5), name = "Cells without N_REF")

    if (!is.null(curves)) {
      f_refreads <- f_refreads %>% add_trace(data = curves,
                                             x = ~get(paste0(dim.title, "_1")),
                                             y = ~get(paste0(dim.title, "_2")),
                                             z = ~get(paste0(dim.title, "_3")),
                                             mode = "lines", type = "scatter3d", split = ~Lineage)
    }
    f_refreads <- f_refreads %>% layout(title = list(text = "N_REF", font = list(color = title_color)),
                                        scene = list(xaxis = list(title = paste0(dim.title, "_1")),
                                                     yaxis = list(title = paste0(dim.title, "_2")),
                                                     zaxis = list(title = paste0(dim.title, "_3"))))

    plots[["N_REF"]] <- f_refreads
    return(plots)
  }


  # function to save plots' json
  plots_json <- lapply(snv_options, function(snv) {
    plots <- generate_snv_plots(snv)
    list(VAF = list(id = paste0("plot_VAF_", gsub(":", "_", snv)), json = plotly_json(plots[["VAF"]])),
         N_VAR = list(id = paste0("plot_N_VAR_", gsub(":", "_", snv)), json = plotly_json(plots[["N_VAR"]])),
         N_REF = list(id = paste0("plot_N_REF_", gsub(":", "_", snv)), json = plotly_json(plots[["N_REF"]])))
  })
  plots_json <- unlist(plots_json, recursive = FALSE)

  if (save_each_plot && !is.null(output_dir)) {
    for (snv in snv_options) {
      plots <- generate_snv_plots(snv)
      save_snv_plot <- function(plot, snv, plot_type) {
        snv_clean <- gsub(":", "_", snv)
        file_path <- file.path(output_dir, "Figures_Individual_Pots_HTML",
                               "Individual_sceSNVs", plot_type, paste0(plot_type, "_", snv_clean, ".html"))

        dir.create(dirname(file_path), showWarnings = FALSE,
                   recursive = TRUE)
        saveWidget(as_widget(plot), file = file_path,
                   selfcontained = FALSE, libdir = "lib")
      }

      save_snv_plot(plots[["VAF"]], snv, "VAF")
      save_snv_plot(plots[["N_VAR"]], snv, "N_VAR")
      save_snv_plot(plots[["N_REF"]], snv, "N_REF")

    }
  }
  snv_out <- list()
  snv_out[["plots_json"]] <- plots_json
  snv_out[["snv_options"]] <- snv_options
  return(snv_out)
}

