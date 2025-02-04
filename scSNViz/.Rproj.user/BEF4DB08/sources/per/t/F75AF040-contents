#' Preprocess scSNViz Data
#'
#' This function preprocesses SNV and Seurat object data, performs dimensionality reduction,
#' clustering, and computes various statistics for visualization.
#'
#' @importFrom Seurat RunPCA RunUMAP RunTSNE
#' @importFrom Seurat Embeddings FindClusters FindNeighbors DefaultAssay DefaultAssay<-
#' @importFrom Seurat CreateSeuratObject PercentageFeatureSet SCTransform
#'
#' @param rds_file Path to the RDS file containing a Seurat object (required if countsmatrix_file is not provided).
#' @param countsmatrix_file Path to the folder containing STARsolo output files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz).
#' @param snv_file Path to the SNV file (required).
#' @param dimensionality_reduction The dimensionality reduction method (umap, pca, or tsne). Default: umap.
#' @param th_vars Threshold for the number of SNVs per cell. Default: 0.
#' @param th_reads Threshold for the number of variant reads per locus. Default: 0.
#' @param enable_sctype Logical; enable scType analysis. Default: FALSE.
#' @param tissue_type Tissue type for scType. Required if enable_sctype is TRUE.
#' @param color_scale Color scale for plots. Options include Blues, Reds, etc. Default: YlOrRd.
#' @return A Seurat object with processed metadata and embeddings.
#' @details
#' This function reads and processes data from an RDS file or STARsolo counts matrix, along with an SNV file.
#' - **RDS File**: A Seurat object is read and processed if provided.
#' - **STARsolo Counts Matrix**: If provided, it processes raw counts into a Seurat object.
#' - **Dimensionality Reduction**: Options include UMAP (umap), PCA (pca), or t-SNE (tsne).
#' - **Thresholding**: Filters cells based on th_vars (number of SNVs) and th_reads (variant reads per locus).
#' - **scType Analysis**: If enable_sctype is TRUE, cell-type annotation is performed using the specified tissue_type.
#'
#' @examples
#' # Example usage:
#' processed_data <- preprocess_snv_data(
#'   rds_file = "path/to/seurat.rds",
#'   snv_file = "path/to/snv_file.tsv",
#'   dimensionality_reduction = "umap",
#'   th_vars = 1,
#'   th_reads = 2,
#'   enable_sctype = TRUE,
#'   tissue_type = "Immune system"
#' )
#'
#' # Accessing the results:
#' seurat_object <- processed_data$SeuratObject
#' plot_data <- processed_data$PlotData
#' processed_snv <- processed_data$ProcessedSNV
#'
#' @export
#'
preprocess_snv_data <- function(rds_file = NULL, countsmatrix_file = NULL, snv_file,
                                dimensionality_reduction = "UMAP", th_vars = 0, th_reads = 0,
                                enable_sctype = FALSE, tissue_type = NULL, color_scale = "YlOrRd") {

  # validate inputs
  if (is.null(rds_file) & is.null(countsmatrix_file)) {
    stop("The rds_file and countsmatrix_file must be provided.")
  }
  if (is.null(snv_file)) {
    stop("The SNV file is required.")
  }

  valid_reductions <- c("umap", "pca", "tsne")
  dimensionality_reduction <- tolower(dimensionality_reduction)
  if (!dimensionality_reduction %in% valid_reductions) {
    stop("Invalid dimensionality_reduction method. Please use one of: 'umap', 'pca', 'tsne'.")
  }

  dim.plotting <- switch(dimensionality_reduction, umap = "UMAP",
                         pca = "PCA", tsne = "tSNE")

  if (enable_sctype && is.null(tissue_type)) {
    stop("Tissue type is required when enable_sctype is TRUE.")
  }
  if (!is.null(tissue_type)) {
    tissue_type <- str_to_title(tolower(tissue_type))
    if (tissue_type %in% c("Immune", "Immunesystem")) {
      tissue_type <- "Immune system"
    }
  }

  # validate tissue types
  valid_tissue_types <- c("Immune system", "Pancreas", "Liver",
                          "Eye", "Kidney", "Brain", "Lung", "Adrenal", "Heart",
                          "Intestine", "Muscle", "Placenta", "Spleen", "Stomach",
                          "Thymus")

  if (enable_sctype && !(tissue_type %in% valid_tissue_types)) {
    stop("Invalid tissue type. Available tissues are: ",
         paste(valid_tissue_types, collapse = ", "))
  }
  if (!is.null(rds_file)) {
    srt <- readRDS(rds_file)
  }
  else {
    gene.matrix <- Read10X(data.dir = countsmatrix_file)
    srt <- CreateSeuratObject(counts = gene.matrix, project = "Sample",
                              min.cells = 3, min.features = 200)
    srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-")
    srt <- SCTransform(object = srt, vst.flavor = "v2", method = "glmGamPoi",
                       vars.to.regress = "percent.mt", verbose = FALSE)
    DefaultAssay(srt) <- "SCT"
  }


  # read the SNV file and perform filtering
  snv <- read.table(snv_file, sep = "\t", header = TRUE)

  if (nrow(snv) == 0) {
    stop("The SNV file is empty.")
  }

  snv <- snv %>% filter(SNVCount > 0)
  snv$VAF[snv$VAF == "-"] <- NA
  snv$VAF <- as.numeric(snv$VAF)
  snv.n <- aggregate(SNVCount ~ ReadGroup, data = snv, length)
  rownames(snv.n) <- snv.n$ReadGroup
  colnames(snv.n)[2] <- "SNV.N"
  snv <- merge(snv, snv.n, "ReadGroup", all.x = TRUE)

  if (th_vars > 0) {
    snv <- snv %>% filter(SNV.N >= th_vars)
  }
  if (th_reads > 0) {
    snv <- snv %>% filter(SNVCount >= th_reads)
  }


  snv.reads <- aggregate(SNVCount ~ ReadGroup, data = snv, sum)
  if (max(snv.reads[, "SNVCount"]) == 0) {
    stop("Even though you have provided a tab-separated file containing SNV information, there aren't any SNV read counts in it.")
  }

  rownames(snv.reads) <- snv.reads$ReadGroup
  ref.reads <- aggregate(RefCount ~ ReadGroup, data = snv, sum)
  rownames(ref.reads) <- snv.reads$ReadGroup

  if (max(snv.reads[, "SNVCount"]) < th_reads) {
    stop("There aren't any SNVs that satisfy the read threshold.")
  }
  if (max(snv.n[, "SNV.N"]) < th_vars) {
    stop("There aren't any SNVs that satisfy the variant threshold.")
  }

  vaf.median <- aggregate(VAF ~ ReadGroup, data = snv, median)
  rownames(vaf.median) <- vaf.median$ReadGroup
  colnames(vaf.median)[2] <- "Median.VAF"
  vaf.mean <- aggregate(VAF ~ ReadGroup, data = snv, mean)
  rownames(vaf.mean) <- vaf.mean$ReadGroup
  colnames(vaf.mean)[2] <- "Mean.VAF"
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


  if (tolower(dimensionality_reduction) == "tsne") {
    srt <- RunPCA(srt)
    srt <- RunTSNE(srt, dim.embed = 3)
    dim.plotting <- "tSNE"
  }
  else if (tolower(dimensionality_reduction) == "pca") {
    srt <- RunPCA(srt)
    srt <- RunUMAP(srt, dims = 1:20, n.components = 3)
    dim.plotting <- "PCA"
  }
  else {
    srt <- RunPCA(srt)
    srt <- RunUMAP(srt, dims = 1:20, n.components = 3)
    dim.plotting <- "UMAP"
  }


  # clustering
  srt <- FindNeighbors(srt, dims = 1:10)
  srt <- FindClusters(srt, resolution = 0.5)


  # run scType
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
                       es.max.cl <- sort(rowSums(es.max[, rownames(srt@meta.data[srt@meta.data$seurat_clusters == cl, ])]), decreasing = TRUE)
                       head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(srt@meta.data$seurat_clusters == cl)), 10)}))
    sctype_scores <- results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
    srt@meta.data$customclassif <- ""
    for (cl in unique(sctype_scores$cluster)) {
      cl_type <- sctype_scores[sctype_scores$cluster == cl, ]
      srt@meta.data$customclassif[srt@meta.data$seurat_clusters == cl] <- as.character(cl_type$type[1])
    }
  }
  srt[["HasSNV"]] <- sapply(srt[["SNVCount"]][, 1], function(x) if (x >= th_reads) 1 else 0)


  if (enable_sctype) {
    df.snv <- as.data.frame(srt[[c("SNVCount", "RefCount", "TotalVAF", "MedianVAF", "MeanVAF",
                                   "SNV.N", "seurat_clusters", "customclassif", "HasSNV", "Undetected")]])
    df.dim <- as.data.frame(Embeddings(srt, reduction = tolower(dimensionality_reduction)))
    df.3dplot <- merge(df.snv, df.dim, by = 0)
    colnames(df.3dplot)[12:14] <- c(paste0(dim.plotting, "_1"), paste0(dim.plotting, "_2"), paste0(dim.plotting, "_3"))
  }
  else {
    df.snv <- as.data.frame(srt[[c("SNVCount", "RefCount", "TotalVAF", "MedianVAF", "MeanVAF",
                                   "SNV.N", "seurat_clusters", "HasSNV", "Undetected")]])
    df.dim <- as.data.frame(Embeddings(srt, reduction = tolower(dimensionality_reduction)))
    df.3dplot <- merge(df.snv, df.dim, by = 0)
    colnames(df.3dplot)[11:13] <- c(paste0(dim.plotting, "_1"), paste0(dim.plotting, "_2"), paste0(dim.plotting, "_3"))
  }

  rownames(df.3dplot) <- df.3dplot$Row.names
  df.3dplot <- df.3dplot[, 2:ncol(df.3dplot)]


  return(list(SeuratObject = srt, PlotData = df.3dplot, ProcessedSNV = snv))
}
