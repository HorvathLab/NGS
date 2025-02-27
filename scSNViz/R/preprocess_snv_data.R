#' Preprocess scSNViz Data
#'
#' This function preprocesses SNV and Seurat object data, performs dimensionality reduction,
#' clustering, computes various statistics for visualization, and optionally generates statistics.
#'
#' @importFrom Seurat RunPCA RunUMAP RunTSNE
#' @importFrom Seurat Embeddings FindClusters FindNeighbors DefaultAssay DefaultAssay<-
#' @importFrom Seurat CreateSeuratObject PercentageFeatureSet SCTransform
#' @importFrom HGNChelper checkGeneSymbols
#' @importFrom stringr str_to_title
#' @importFrom dplyr top_n filter arrange mutate group_by
#' @importFrom stats kruskal.test p.adjust
#'
#' @param rds_obj Processed Seurat object.
#' @param snv_file SNV file (required).
#' @param dimensionality_reduction The dimensionality reduction method (umap, pca, or tsne). Default: umap.
#' @param th_vars Threshold for number of sceSNVs for a cell to be displayed. Default: 0.
#' @param th_reads Threshold for number of variant reads (N_VAR) for a locus to be considered sceSNV. Default: 0.
#' @param enable_sctype Logical; enable scType analysis for cell types. Default: FALSE.
#' @param enable_integrated Logical; enable scType analysis for cell types. Default: FALSE.
#' @param tissue_type Tissue type for scType. Required if enable_sctype is TRUE.
#' @param color_scale Color scale for plots. Options include Blues, Reds, etc. Default: YlOrRd.
#' @param generate_statistics Logical; if TRUE, generate SNV significance statistics. Default: FALSE.
#' @param th_snv_cells Threshold for maximum percentage of cells that contain an SNV for bad reads. Default: 10
#' @param output_dir Directory where the statistics files will be saved if `generate_statistics` is TRUE. Default: Current working directory.
#' @return A Seurat object with processed metadata and embeddings, processed SNV data, aggregated SNV data, and a data frame for plotting.
#'
#' @examples
#' processed_data <- preprocess_snv_data(
#'   rds_obj = srt,
#'   snv_file = "path/to/snv_file.tsv",
#'   dimensionality_reduction = "umap",
#'   th_vars = 1,
#'   th_reads = 2,
#'   enable_sctype = T,
#'   tissue_type = "Immunesystem",
#'   generate_statistics = T,
#'   th_snv_cells = 10,
#'   output_dir = "path/to/output"
#' )
#'
#' @export
#'
preprocess_snv_data <- function(rds_obj = NULL, snv_file = NULL,
                                dimensionality_reduction = "UMAP", th_vars = 0, th_reads = 0,
                                enable_sctype = F, tissue_type = NULL, enable_integrated = F, integrated_reduction_name='integrated', color_scale = "YlOrRd",
                                generate_statistics = F, th_snv_cells = 10, output_dir = NULL) {

  # Validate inputs
  if (is.null(rds_obj)) {
    stop("A Seurat object is required for the rds_obj parameter.")
  }
  if (is.null(snv_file)) {
    stop("The scReadCounts file is required for the snv_file parameter.")
  }

  if (generate_statistics && is.null(output_dir)) {
    output_dir <- getwd()
  }

  valid_reductions <- c("umap", "pca", "tsne")
  dimensionality_reduction <- tolower(dimensionality_reduction)
  if (!dimensionality_reduction %in% valid_reductions) {
    stop("You have input a dimensionality reduction that is not supported by this package. Please use pca, umap, or tsne.")
  }
  else {
    if (dimensionality_reduction == 'tsne'){
      dim.plotting = 'TSNE'
      dim.title = 'TSNE'
    }
    if (dimensionality_reduction == 'pca'){
      dim.plotting = 'PCA'
      dim.title = 'PCA'
    }
    if (dimensionality_reduction == 'umap'){
      dim.plotting = 'UMAP'
      dim.title = 'UMAP'
    }
  }

  if (enable_integrated && dimensionality_reduction!='umap') {
    stop('You can only render UMAPs with integrated data, not PCA or TSNE plots. Please leave the dimensionality_reduction parameter as the default, which is UMAP.')
  }
  
  # Cell type analysis by scType
  if (enable_sctype) {
    if (is.null(tissue_type)) {
      stop("- Tissue type is required when enable_sctype is TRUE.")
    } else {
      # Normalize and validate tissue type
      tissue_type <- str_to_title(tolower(tissue_type))
      if (tissue_type %in% c("Immune", "Immunesystem")) {
        tissue_type <- "Immune system"
      }
      valid_tissue_types <- c("Immune system", "Pancreas", "Liver", "Eye", "Kidney",
                              "Brain", "Lung", "Adrenal", "Heart", "Intestine", "Muscle",
                              "Placenta", "Spleen", "Stomach", "Thymus")
      if (!(tissue_type %in% valid_tissue_types)) {
        stop(paste("- The tissue type you have provided does not seem to be one of the tissue types that scType accepts:",
                   paste(valid_tissue_types, collapse = ", ")))
      }
    }
  }

  
  # Read in the read counts data
  srt <- rds_obj
  # Perform dimensionality reduction
  if (enable_integrated){
      srt <- RunUMAP(srt, dims = 1:20, n.components = 3, reduction=integrated_reduction_name)
  } else {
    if (tolower(dimensionality_reduction) == "tsne") {
      srt <- RunTSNE(srt, dim.embed = 3)
    }
    else if (tolower(dimensionality_reduction) == "umap") {
      srt <- RunUMAP(srt, dims = 1:20, n.components = 3)}
  }
  
  # Read and process SNV data
  snv_statistics <- function(snv,th.vars=th_vars,th.reads=th_reads){
    if (nrow(snv) == 0) {
      stop("The SNV file is empty.")
    }
    snv.modify <- snv %>% filter(SNVCount > 0)
    num.snvs <- nrow(snv.modify)
    snv$VAF[snv$VAF == "-"] <- NA
    snv$VAF <- as.numeric(snv$VAF)

    snv.n <- aggregate(SNVCount ~ ReadGroup, data = snv, length)
    rownames(snv.n) <- snv.n$ReadGroup
    colnames(snv.n)[2] <- "SNV.N"
    snv <- merge(snv, snv.n, "ReadGroup", all.x = T)

    if (th.vars > 0) {
      snv <- snv %>% filter(SNV.N >= th.vars)
    }
    if (th.reads > 0) {
      snv <- snv %>% filter(SNVCount >= th.reads)
    }

    snv.reads <- aggregate(SNVCount ~ ReadGroup, data = snv, sum)
    if (max(snv.reads[, "SNVCount"]) == 0) {
      stop("Even though you have provided a tab-separated file containing SNV information, there aren't any SNV read counts in it.")
    }

    rownames(snv.reads) <- snv.reads$ReadGroup
    ref.reads <- aggregate(RefCount ~ ReadGroup, data = snv, sum)
    rownames(ref.reads) <- snv.reads$ReadGroup

    if (max(snv.reads[, "SNVCount"]) < th.reads) {
      stop("There aren't any SNVs that satisfy the read threshold.")
    }
    if (max(snv.n[, "SNV.N"]) < th.vars) {
      stop("There aren't any SNVs that satisfy the variant threshold.")
    }

    vaf.median <- aggregate(VAF ~ ReadGroup, data = snv, median)
    rownames(vaf.median) <- vaf.median$ReadGroup
    colnames(vaf.median)[2] <- "Median.VAF"
    vaf.mean <- aggregate(VAF ~ ReadGroup, data = snv, mean)
    rownames(vaf.mean) <- vaf.mean$ReadGroup
    colnames(vaf.mean)[2] <- "Mean.VAF"
    
    return(list(snv=snv,vaf.median=vaf.median,vaf.mean=vaf.mean,snv.reads=snv.reads,ref.reads=ref.reads))
  }
  snv_file_tble <- read.table(snv_file, sep = "\t", header = T)
  snv = data.frame()
  vaf.mean = data.frame()
  vaf.median = data.frame()
  snv.reads = data.frame()
  ref.reads = data.frame()
  if (length(unique(srt$orig.ident))>1){
      for (i in 1:length(unique(srt$orig.ident))){
      #data must be integrated
        snv_file_tble$sampleid = data.frame(do.call('rbind',strsplit(as.character(snv_file_tble$ReadGroup),'_',fixed=TRUE)))$X1
        snv_file_subset <- snv_file_tble[snv_file_tble$sampleid == unique(srt$orig.ident)[i],]
        snv_file_subset$sampleid  = NULL
        output_list <- snv_statistics(snv=snv_file_subset)
        # combine statistics with other values from other samples
        snv=rbind(snv,output_list$snv)
        vaf.median=rbind(vaf.median,output_list$vaf.median)
        vaf.mean=rbind(vaf.mean,output_list$vaf.mean)
        snv.reads=rbind(snv.reads,output_list$snv.reads)
        ref.reads=rbind(ref.reads,output_list$ref.reads)
      }} else { #data not integrated
      snv_file_subset <- snv_file_tble
      output_list <- snv_statistics(snv_file_subset)
      snv=output_list$snv
      vaf.median=output_list$vaf.median
      vaf.mean=output_list$vaf.mean
      snv.reads=output_list$snv.reads
      ref.reads=output_list$ref.reads
  }


  # Combine read counts and snv data
  srt[["SNVCount"]] <- 0
  srt[["SNVCount"]][rownames(snv.reads), 1] <- snv.reads[, "SNVCount"]
  srt[["RefCount"]] <- 0
  srt[["RefCount"]][rownames(ref.reads), 1] <- ref.reads[, "RefCount"]
  srt[["TotalVAF"]] <- 0
  srt[["TotalVAF"]] <- srt[["SNVCount"]] / (srt[["SNVCount"]] + srt[["RefCount"]])
  srt[["TotalVAF"]][!is.finite(srt[["TotalVAF"]][, 1]), 1] <- 0
  srt[["MedianVAF"]] <- 0
  srt[["MedianVAF"]][rownames(vaf.median), 1] <- vaf.median[, "Median.VAF"]
  srt[["MedianVAF"]][!is.finite(srt[["MedianVAF"]][, 1]), 1] <- 0
  srt[["MeanVAF"]] <- 0
  srt[["MeanVAF"]][rownames(vaf.mean), 1] <- vaf.mean[, "Mean.VAF"]
  srt[["MeanVAF"]][!is.finite(srt[["MeanVAF"]][, 1]), 1] <- 0
  snv.n <- unique(snv[, c("ReadGroup", "SNV.N")])
  rownames(snv.n) <- snv.n[, "ReadGroup"]
  srt[["SNV.N"]] <- 0
  srt[["SNV.N"]][rownames(snv.n), 1] <- snv.n[, "SNV.N"]
  srt[["SNV.N"]][!is.finite(srt[["SNV.N"]][, 1]), 1] <- 0
  srt[["Undetected"]] <- ifelse((srt[["SNVCount"]] == 0 & srt[["RefCount"]] == 0), 1, 0)
  srt[["HasSNV"]] <- sapply(srt[["SNVCount"]][, 1], function(x) if (x >= th_reads) 1 else 0)


generate_statistics_fnction <- function(snv,th.snv.cells=th_snv_cells){
  
    snv <- snv[!is.na(snv$VAF), ]
    # filter based on `X.BadRead` and `SNVCount`
    snv$temp <- snv
    snv$temp$BadReadFlag <- 0
    snv$temp$BadReadFlag[snv$temp$`X.BadRead` > 0] <- 1
    snv.read.filt <- aggregate(BadReadFlag ~ CHROM + POS + REF + ALT, data = snv$temp,
                               function(x) 100 * sum(x) / length(x))
    snv.read.filt <- snv.read.filt[snv.read.filt$BadReadFlag <= th.snv.cells, ]
    snv <- merge(snv, snv.read.filt, by = c("CHROM", "POS", "REF", "ALT"))
    snv$VAF[snv$SNVCount < th_reads] <- 0
    if (nrow(snv) == 0) {
      stop("There are no rows left after filtering. You may need to reset your threshold parameters.")
    }

    # cluster assignment
    df.cid <- as.data.frame(srt[["seurat_clusters"]])
    df.cid <- data.frame(ReadGroup = rownames(df.cid), ClusterID = df.cid[, 1], row.names = NULL)
    df.snv <- merge(snv, df.cid, by = "ReadGroup")

    snv_groups <- split(df.snv, paste(df.snv$CHROM, df.snv$POS, df.snv$REF, df.snv$ALT, sep = "_"))

    # Kruskal-Wallis test for each SNV across clusters
    kw_test_statistics <- c()
    kw_p_values <- c()

    for (snv_key in names(snv_groups)) {
      snv_data <- snv_groups[[snv_key]]
      if (length(unique(snv_data$ClusterID)) > 1) {
        kw_test_result <- kruskal.test(VAF ~ ClusterID, data = snv_data)
        kw_test_statistics <- c(kw_test_statistics, kw_test_result$statistic)
        kw_p_values <- c(kw_p_values, kw_test_result$p.value)
      } else {
        kw_test_statistics <- c(kw_test_statistics, NA)
        kw_p_values <- c(kw_p_values, NA)
      }
    }

    kw_p_adj <- p.adjust(kw_p_values, method = "bonferroni")

    df.final <- data.frame(
      SNV = names(snv_groups),
      H_statistic = kw_test_statistics,
      p_value = kw_p_values,
      p_adj = kw_p_adj
    )
    df.final <- na.omit(df.final)
    df.final <- df.final[order(df.final$p_adj), ]
    return(df.final)
}
  # generate SNV statistics if requested
  if (generate_statistics) {
    dir.create(output_dir, recursive = T, showWarnings = F)
    if (enable_integrated){
      for (i in 1:length(unique(srt$orig.ident))){
        snv$sampleid = data.frame(do.call('rbind',strsplit(as.character(snv$ReadGroup),'_',fixed=TRUE)))$X1
        snv_input = snv[snv$sampleid == unique(srt$orig.ident)[i],]
        df.final <- generate_statistics_fnction(snv=snv_input)
        write.table(df.final, file = file.path(
        output_dir, paste0("SNV_Statistics_",unique(srt$orig.ident)[i],".txt")), sep = "\t", row.names = F)
        significant_snvs <- df.final[df.final$p_adj < 0.05, ]
        write.table(significant_snvs, file = file.path(
        output_dir, paste0("Significant_SNVs_",unique(srt$orig.ident)[i],".txt")), sep = "\t", row.names = F)
      }
    } else {
      df.final <- generate_statistics_fnction(snv)
      write.table(df.final, file = file.path(
      output_dir, "SNV_Statistics.txt"), sep = "\t", row.names = F)
      significant_snvs <- df.final[df.final$p_adj < 0.05, ]
      write.table(significant_snvs, file = file.path(
      output_dir, "Significant_SNVs.txt"), sep = "\t", row.names = F)
    }
  }

  # run scType analysis if requested
  if (enable_sctype) {
    cat("Running scType...\n")
    db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    gs_list <- tryCatch({
      gene_sets_prepare(db_url, tissue_type)
    }, error = function(e) {
      stop("Error in scType database processing: ", e$message,
           "\nCheck if the database is accessible and matches the tissue type.")
    })
    if (length(gs_list$gs_positive) == 0 || length(gs_list$gs_negative) == 0) {
      stop("Gene sets for the selected tissue type are empty. Check tissue type or database.")
    }
    es.max <- sctype_score(scRNAseqData = srt[[DefaultAssay(srt)]]@scale.data,
                           gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    results <- do.call("rbind", lapply(unique(srt@meta.data$seurat_clusters),function(cl) {
                       es.max.cl <- sort(rowSums(es.max[, rownames(srt@meta.data[srt@meta.data$seurat_clusters == cl, ])]),
                                         decreasing = T)
                       head(data.frame(cluster = cl, type = names(es.max.cl),
                                       scores = es.max.cl, ncells = sum(srt@meta.data$seurat_clusters == cl)), 10)}))
    sctype_scores <- results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
    srt@meta.data$customclassif <- ""
    for (cl in unique(sctype_scores$cluster)) {
      cl_type <- sctype_scores[sctype_scores$cluster == cl, ]
      srt@meta.data$customclassif[srt@meta.data$seurat_clusters == cl] <- as.character(cl_type$type[1])
    }
  }

  if (enable_sctype) {
    df.snv <- as.data.frame(srt[[c("SNVCount", "RefCount", "TotalVAF", "MedianVAF",
                                   "MeanVAF", "SNV.N", "seurat_clusters", "customclassif",
                                   "HasSNV", "Undetected","orig.ident")]])
    df.dim <- as.data.frame(Embeddings(srt, reduction = tolower(dimensionality_reduction)))
    df.3dplot <- merge(df.snv, df.dim, by = 0)
    colnames(df.3dplot)[13:15] <- c(paste0(dim.plotting, "_1"),
                                    paste0(dim.plotting, "_2"),
                                    paste0(dim.plotting, "_3"))
  }
  else {
    df.snv <- as.data.frame(srt[[c("SNVCount", "RefCount", "TotalVAF", "MedianVAF", "MeanVAF",
                                   "SNV.N", "seurat_clusters", "HasSNV", "Undetected","orig.ident")]])
    df.dim <- as.data.frame(Embeddings(srt, reduction = tolower(dimensionality_reduction)))
    df.3dplot <- merge(df.snv, df.dim, by = 0)
    colnames(df.3dplot)[12:14] <- c(paste0(dim.plotting, "_1"),
                                    paste0(dim.plotting, "_2"),
                                    paste0(dim.plotting, "_3"))
  }

  rownames(df.3dplot) <- df.3dplot$Row.names
  df.3dplot <- df.3dplot[, 2:ncol(df.3dplot)]
  plot_data <- df.3dplot

  return(list(SeuratObject = srt, ProcessedSNV = snv,
              AggregatedSNV = df.snv, PlotData = plot_data))
}
