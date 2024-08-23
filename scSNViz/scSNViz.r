# Aug 20, 2024
suppressPackageStartupMessages({
  library('optparse')
  library('stringr')
})

script.desc <-
'This script calculates and plots basic statistics, 2-dimensional and
3-dimensional tSNEs of SNV reads, reference reads, and VAF for each cell. The
user must supply the RDS of a Seurat object (-r) which can contain multiple
integrated datasets. The user must also supply an SNV file (-t) that contains
at least one SNV. The RDS file needs to have been dimensionally reduced and
clustered. Dark gray values in plots represent NAs.'

parser <- OptionParser(description=script.desc)

parser <- add_option(parser, c('-r', '--rds-file'),
                     type='character',
                     help='RDS file containing Seurat object.')
parser <- add_option(parser, c('-m', '--countsmatrix-file'),
                     type='character',
                     help='folder containing STARsolo output folder name that contains
                     the following files:
                     barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz')
parser <- add_option(parser, c('-t', '--snv-file'),
                     type='character',
                     help='scReadCounts file')
parser <- add_option(parser, c('-w', '--dimensionality-reduction'),
                     type='character', default='umap',
                     help='options include tSNE, PCA, UMAP')
parser <- add_option(parser, c('-x', '--th-vars'),
                     type='integer', default=0,
                     help='Threshold for number of sceSNVs for a cell to be displayed. Default=0 (display cells with N_SNVs > 0).')
parser <- add_option(parser, c('-y', '--th-reads'),
                     type='integer', default=0,
                     help='Threshold for number of variant reads (N_VAR) for a locus to be considered sceSNV. Default=0 (consider as sceSNV positions covered with N_VAR > 0).')
parser <- add_option(parser, c('-c', '--disable-title'),
                     type='logical', default=F, action='store_true',
                     help='Disable title for individual SNV plots. Default=F')
parser <- add_option(parser, c('-d', '--disable-ind-plots'),
                     type='logical', default=F, action='store_true',
                     help='Disable individual SNV plots. Default=F.')
parser <- add_option(parser, c('-e', '--disable-3d-axis'),
                     type='logical', default=F, action='store_true',
                     help='Disable axes in 3D plots. Default=F.')
parser <- add_option(parser, c('-g', '--disable-slingshot'),
                     type='logical', default=F, action='store_true',
                     help='Disable slingshot curves in 3D plots. Default=F.')
parser <- add_option(parser, c('-i', '--enable-sctype'),
                     type='logical', default=F, action='store_true',
                     help='Enable scType to run. Default=F.')
parser <- add_option(parser, c('-j', '--tissue-type'),
                     type='character',
                     help='tissue type for scType; options include:
                     Immunesystem, Pancreas, Liver, Eye, Kidney, Brain,
                     Lung, Adrenal, Heart, Intestine, Muscle, Placenta,
                     Spleen, Stomach, Thymus')
parser <- add_option(parser, c('-k', '--color-scale'),
                     type='character', default='YlOrRd',
                     help='if you would like to change the default color settings with
                     these options, you may use Blues, Reds, YlOrRd, YlGnBu, plasma, RdBu')
parser <- add_option(parser, c('-b', '--enable-cell-border'),
                     type='logical', default=F, action='store_true',
                     help='Enable cell border. Default=F')
parser <- add_option(parser, c('-q', '--enable-dynamic-cell-size'),
                     type='logical', default=F, action='store_true',
                     help='Enable cell size to depend on number of reads. Default=F')
parser <- add_option(parser, c('-u', '--enable-copykat'),
                     type='logical', default=F, action='store_true',
                     help='Enable CopyKat for displaying CNVs. Default=F')
parser <- add_option(parser, c('-s', '--save-each-plot'),
                     type='logical', default=F, action='store_true',
                     help='Save plots as separate HTML files. Default=F')
args <- parse_args(parser)

error.msg <- NULL


# Check if the required argument (-r) is passed
if (is.null(args$`rds-file`) & is.null(args$`countsmatrix-file`))
  error.msg <- paste(error.msg, '- Seurat RDS object (-r) or STARsolo output directory for counts matrix is required.', sep='\n')
if (is.null(args$`snv-file`))
  error.msg <- paste(error.msg, '- scReadCounts file (-t) is required.', sep='\n')
if (!(is.null(args$`tissue-type`)) & !(args$`enable-sctype`)){
   error.msg <- paste(error.msg, 'You have provided a tissue type. However, to run sctype, you must also type -i into the command line')}

if (args$`enable-sctype`) {
  if (is.null(args$`tissue-type`)) {
    error.msg <- paste(error.msg, '- tissue type is required')
  } 
  else {
    tissue.type <- args$`tissue-type`
    tissue.type <- str_to_title(tolower(tissue.type))
    if (!(tissue.type %in% list('Immunesystem','Immune', 'Pancreas', 'Liver', 'Eye', 
          'Kidney', 'Brain', 'Lung', 'Adrenal', 'Heart', 'Intestine', 'Muscle',
          'Placenta', 'Spleen','Stomach','Thymus'))) {
      error.msg <- paste(error.msg, '- the tissue type you have provided does not seem to be one of the tissue types that scType accepts', sep='\n')
    }
    if (tissue.type=='Immune' | tissue.type=='Immunesystem'){
      tissue.type='Immune system'
    }
  }
}

if (!is.null(args$`color-scale`)){
  color.scale <- args$`color-scale`
  if (!(color.scale %in% list('Blues', 'Reds', 'YlOrRd', 'YlGnBu', 'plasma','RdBu'))) {
    error.msg <- paste(error.msg, '- the color scale you provided is either not in the acceptable form or is not one this program uses')
  }
}

cell.size = 0
if (args$`enable-cell-border`) {
  cell.size=3
}

if (is.null(args$`dimensionality-reduction`)){
  error.msg <- paste(error.msg,'You did not provide a dimensionality reduction selection. The default is UMAP. To change this in the future, look at option -w')
} else {
  dimensionality.reduction <- args$`dimensionality-reduction`
  dimensionality.reduction <- tolower(dimensionality.reduction)
  if (dimensionality.reduction!='pca' & dimensionality.reduction!='umap' & dimensionality.reduction!='tsne'){
    error.msg <- paste(error.msg, 'You inputted a dimensionality reduction that is not supported by this software. Please use pca, umap, or tsne.')
  } else {
  if (dimensionality.reduction=='tsne'){
    dim.plotting = 'tSNE'
    dim.title = 'tSNE'
  }
  if (dimensionality.reduction=='pca'){
    dim.plotting = 'PC'
    dim.title = 'PCA'
  }
  if (dimensionality.reduction=='umap'){
    dim.plotting = 'UMAP'
    dim.title = 'UMAP'
  }
  }
}

if (args$`enable-copykat`){
  cat("\nCopyKat has been enabled. Depending on the number of cells, it may take a while to finish executing...\n\n")
}

# Check if there are any errors
if (!is.null(error.msg)) {
  print_help(parser)
  stop(error.msg)}

######
# Application Logic
######

suppressPackageStartupMessages({
  library('Seurat')
  library('ggplot2')
  library('dplyr')
  library('openxlsx')
  library('HGNChelper')
})

# load files
snv.file <- args$`snv-file`
th.vars <- args$`th-vars`
th.reads <- args$`th-reads`

sample.name <- sub('^(.*/)*(.*)\\.([a-z0-9]+)$', '\\2', ignore.case=T, snv.file)
param.name <- paste0(paste0('_',dimensionality.reduction,'_'), th.vars, 'r', th.reads)

if (is.null(args$`countsmatrix-file`)){
    rds.file <- args$`rds-file`
    srt <- readRDS(rds.file)
} else {
    countsmatrix.file <- args$`countsmatrix-file`
    gene.matrix <- Read10X(data.dir = countsmatrix.file)
    srt <- CreateSeuratObject(counts = gene.matrix, project = sample.name, min.cells = 3, min.features = 200)
    srt[['percent.mt']] <- PercentageFeatureSet(srt, pattern= '^MT-')
    srt <- SCTransform(object = srt, vst.flavor = 'v2', method = "glmGamPoi", verbose = FALSE, vars.to.regress = "percent.mt")
    DefaultAssay(srt) <- 'SCT'
}

######
# Color Scales
######

if (is.null(args$`color-scale`)){
    color.scale3 <- 'YlOrRd'                   # color scales for 2D and 3D plots; simple alternatives include: Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu YlOrBr, YlOrRd
    histogram.scale1<-'tomato'                 # for histograms and combined plots
    histogram.scale2<-'dodgerblue2'            # for histograms and combined plots
    color.undetected <- 'tomato'               # color for cells where the SNV/s is/are undetected
    reversescale.option<-T
} else {
    color.scale3 <- color.scale                # color scales for 2D and 3D plots; simple alternatives include: Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu YlOrBr, YlOrRd
    color.undetected <- 'tomato'               # color for cells where the SNV/s is/are undetected
    histogram.scale1<-'tomato'                 # for histograms and combined plots
    histogram.scale2<-'dodgerblue2'            # for histograms and combined plots
    scale.list<-list('Blues','Reds','YlOrRd','YlGnBu','plasma','RdBu')
    reversescale.options <- list(T,F,T,T,F,F)
    reversescale.option<-reversescale.options[[which(scale.list==color.scale)]]
}

if (dimensionality.reduction=='tsne'){
  srt <- RunPCA(srt)
  srt <- RunTSNE(srt, dim.embed=3)
}
if (dimensionality.reduction=='pca'){
  srt <- RunPCA(srt)
  srt <- RunUMAP(srt, dims=1:20, n.components=3)
}
if (dimensionality.reduction=='umap'){
  srt <- RunPCA(srt)
  srt <- RunUMAP(srt, dims=1:20, n.components=3)
}


# clustering 
srt <- FindNeighbors(srt, dims = 1:10)
srt <- FindClusters(srt, resolution = 0.5)

# SNVs file filtering
snv <- read.table(snv.file, sep='\t', header=T)
snv.modify<- snv %>% filter(SNVCount>0)
num.snvs <- nrow(snv.modify)
# convert '-' to NA in VAF column (here VAF = ∞ or not defined as var, ref = 0)
snv$VAF[snv$VAF == '-'] <- NA
# convert VAF to numeric
snv$VAF <- as.numeric(snv$VAF)


# if scType option is selected, create cell type labels
if (args$`enable-sctype`) {
  cat("Running scType...\n")
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  tissue = tissue.type
  gs_list = gene_sets_prepare(db_, tissue)
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = srt[[DefaultAssay(srt)]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(srt@meta.data$seurat_clusters), function(cl) {
    es.max.cl = sort(rowSums(es.max[, rownames(srt@meta.data[srt@meta.data$seurat_clusters == cl, ])]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(srt@meta.data$seurat_clusters == cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] = "Unknown"
  srt@meta.data$customclassif = ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type = sctype_scores[sctype_scores$cluster == j, ] 
    srt@meta.data$customclassif[srt@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
}


# calculate number of detected snvs per cell
snv.n <- aggregate(SNVCount ~ ReadGroup, data = snv, length)
rownames(snv.n) <- snv.n$ReadGroup
colnames(snv.n)[2] <- 'SNV.N'

snv <- merge(snv, snv.n, 'ReadGroup', all.x = TRUE)  # added snv.n column 

# filters for thresholds (options -x and -y)
if (th.vars > 0) {
    snv <- snv %>% filter(SNV.N >= th.vars)
}
if (th.reads > 0) {
    snv <- snv %>% filter(SNVCount >= th.reads)
}


# calculate number of snv and ref reads per cell
snv.reads <- aggregate(SNVCount ~ ReadGroup, data = snv, sum)
if (max(snv.reads[, 'SNVCount']) == 0) {
  stop('Even though you have provided a tab-separated file containing SNV information, there aren\'t any SNV read counts in it.')
}
rownames(snv.reads) <- snv.reads$ReadGroup
ref.reads <- aggregate(RefCount ~ ReadGroup, data = snv, sum)
rownames(ref.reads) <- snv.reads$ReadGroup

if (max(snv.reads[, 'SNVCount']) < th.reads) {
  stop('There aren\'t any SNVs that satisfy the read threshold.')
}

if (max(snv.n[, 'SNV.N']) < th.vars) {
  stop('There aren\'t any SNVs that satisfy the variant threshold.')
}


suppressWarnings(dir.create(paste0(sample.name, param.name)))
dir.name <- paste0(sample.name, param.name, '/')

# calculate median vaf per cell
vaf.median <- aggregate(VAF ~ ReadGroup, data = snv, median)
rownames(vaf.median) <- vaf.median$ReadGroup
colnames(vaf.median)[2] <- 'Median.VAF'

# calculate mean vaf per cell
vaf.mean <- aggregate(VAF ~ ReadGroup, data = snv, mean)
rownames(vaf.mean) <- vaf.mean$ReadGroup
colnames(vaf.mean)[2] <- 'Mean.VAF'

# set meta-data
srt[['SNVCount']] <- 0
srt[['SNVCount']][rownames(snv.reads), 1] <- snv.reads[, 'SNVCount']

srt[['RefCount']] <- 0
srt[['RefCount']][rownames(ref.reads), 1] <- ref.reads[, 'RefCount']

srt[['TotalVAF']] <- 0
srt[['TotalVAF']] <- srt[['SNVCount']] / (srt[['SNVCount']] + srt[['RefCount']])
srt[['TotalVAF']][!is.finite(srt[['TotalVAF']][, 1]), 1] <- 0

srt[['MedianVAF']] <- 0
srt[['MedianVAF']][rownames(vaf.median), 1] <- vaf.median[, 'Median.VAF']
srt[['MedianVAF']][!is.finite(srt[['MedianVAF']][, 1]), 1] <- 0

srt[['MeanVAF']] <- 0
srt[['MeanVAF']][rownames(vaf.mean), 1] <- vaf.mean[, 'Mean.VAF']
srt[['MeanVAF']][!is.finite(srt[['MeanVAF']][, 1]), 1] <- 0

snv.n <- unique(snv[, c('ReadGroup', 'SNV.N')])
rownames(snv.n) <- snv.n[, 'ReadGroup']
srt[['SNV.N']] <- 0
srt[['SNV.N']][rownames(snv.n), 1] <- snv.n[, 'SNV.N']
srt[['SNV.N']][!is.finite(srt[['SNV.N']][, 1]), 1] <- 0

# adding undetected loci (both SNVcount/VARreads and REFreads = 0)
srt[['Undetected']] <- ifelse((srt[['SNVCount']] == 0 & srt[['RefCount']] == 0), 1, 0) 

# generate dataframes for plots
srt[['HasSNV']] <- sapply(srt[['SNVCount']][, 1], function (x) if (x >= th.reads) 1 else 0)
if (args$`enable-sctype`) {
  df.snv <- as.data.frame(srt[[c('SNVCount', 'RefCount', 'TotalVAF', 
  'MedianVAF', 'MeanVAF', 'SNV.N', 'seurat_clusters', 'customclassif', 'HasSNV', 'Undetected')]])
  df.dim <- as.data.frame(Embeddings(srt, reduction = dimensionality.reduction))
  df.3dplot <- merge(df.snv, df.dim, by = 0)
  colnames(df.3dplot)[12:14] <- c(paste0(dim.plotting, '_1'), paste0(dim.plotting, '_2'), paste0(dim.plotting, '_3'))
} else {
  df.snv <- as.data.frame(srt[[c('SNVCount', 'RefCount', 'TotalVAF',
   'MedianVAF', 'MeanVAF', 'SNV.N', 'seurat_clusters', 'HasSNV', 'Undetected')]])
  df.dim <- as.data.frame(Embeddings(srt, reduction = dimensionality.reduction))
  df.3dplot <- merge(df.snv, df.dim, by = 0)
  colnames(df.3dplot)[11:13] <- c(paste0(dim.plotting, '_1'), paste0(dim.plotting, '_2'), paste0(dim.plotting, '_3'))
}
rownames(df.3dplot) <- df.3dplot$Row.names
df.3dplot <- df.3dplot[, 2:ncol(df.3dplot)]


# CopyKat code
if (args$`enable-copykat`) {
  suppressPackageStartupMessages({
    library('copykat')
  })

    cat('Running CopyKat. This may take a while...\n')
    cts <- GetAssayData(srt, layer='counts', assay='SCT')
    ckt <- copykat(cts, sam.name=sample.name)

    srt[['karyotype']] <- 'Unknown'
    srt[['karyotype']][rownames(ckt$pred), ] <- ckt$pred[, 2]
}


# generate interactive 2D plots
suppressPackageStartupMessages({
  library('plotly')
  library('htmlwidgets')
  library('htmltools')
  library('jsonlite')
})

# generate slingshot trajectory
suppressPackageStartupMessages(library('slingshot'))
sce <- as.SingleCellExperiment(srt)
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = toupper(dimensionality.reduction))

curves <- slingCurves(sce, as.df = TRUE)
colnames(curves)[1:3] <- c(paste0(dim.plotting, '_1'), paste0(dim.plotting, '_2'), paste0(dim.plotting, '_3'))


# create directories (if they doesn't exist) for saving each plots as seperate html
if (args$`save-each-plot`) { 
  dir.create(file.path(dir.name, "Figures_Individual_Plots_HTML"), showWarnings = FALSE)
  if (!args$`disable-ind-plots`) {
    dir.create(file.path(dir.name, "Figures_Individual_Plots_HTML", "Individual_sceSNVs"), showWarnings = FALSE)
    dir.create(file.path(dir.name, "Figures_Individual_Plots_HTML", "Individual_sceSNVs", "VAF"), showWarnings = FALSE)
    dir.create(file.path(dir.name, "Figures_Individual_Plots_HTML", "Individual_sceSNVs", "VARreads"), showWarnings = FALSE)
    dir.create(file.path(dir.name, "Figures_Individual_Plots_HTML", "Individual_sceSNVs", "REFreads"), showWarnings = FALSE)
  }
}


### ==================================================== HISTOGRAMS ========================================================== ###

histogram_names <- c("TotalVAF", "MeanSNVsVAF", "N_VARreadCounts", "N_SNV")
histograms <- list()   # list to save the individual histograms

if (num.snvs > 1) {
  hist_list <- list(
    list(aes = aes(x = TotalVAF), xlab = 'TotalVAF', file_suffix = 'Histogram_TotalVAF'),
    list(aes = aes(x = MeanVAF), xlab = 'MeanSNVsVAF', file_suffix = 'Histogram_MeanSNVsVAF'),
    list(aes = aes(x = SNVCount), xlab = 'N_VARreadCounts', file_suffix = 'Histogram_N_VARreadCounts'),
    list(aes = aes(x = SNV.N), xlab = 'N_SNV', file_suffix = 'Histogram_N_SNV')
  )
  
  for (i in seq_along(hist_list)) {
    hist_info <- hist_list[[i]]
    p <- ggplot(df.snv, hist_info$aes) +
         geom_histogram(fill=histogram.scale2, boundary=0, alpha=0.6) +
         xlab(hist_info$xlab) + ylab('Cells') +
         theme(axis.text = element_text(size = 12),
               axis.title = element_text(size = 14),
               panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
               axis.ticks.x=element_blank(), axis.ticks.y=element_blank())

    histograms[[length(histograms) + 1]] <- ggplotly(p)

    if (args$`save-each-plot`) {
      ggsave(file = file.path(dir.name, "Figures_Individual_Plots_HTML", paste0(hist_info$file_suffix, ".png")), plot = p, device = "png")
    }
  }
}



### ==================================================== DIMENSIONALITY REDUCTION PLOTS ========================================================== ###

# descriptions for each plot title's pop-up 
plot_descriptions <- list(
  "N_sceSNVs" = "Absolute number of sceSNVs per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
  "N_VARreads" = "Sum of the absolute number of reads bearing the variant nucleotide across all sceSNVs in the submitted list per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
  "N_REFreads" = "Sum of the absolute number of reads bearing the reference nucleotide across all sceSNVs in the submitted list per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
  "Total_VAF_RNA" = "Total VAF_RNA is calculated by dividing the sum of N_VAR counts across the sceSNVs (loci) by the total reads covering the sceSNV loci (N_VAR + N_REF), and is intended to provide an assessment of SNV expression magnitude across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
  "Mean_VAF_RNA" = "Mean VAF_RNA is calculated as the mean of the VAF_RNA across all individual sceSNVs and is intended to provide an assessment of the variability of the SNV expression levels across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
  "Median_VAF_RNA" = "Median VAF_RNA is calculated as the median of the VAF_RNA across all individual sceSNVs and is intended to provide an assessment of the central tendency of the SNV expression levels across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
  "Cell types (scType)" = "Classification of individual cells into distinct types using the scType tool. scType employs predefined marker gene sets and a scoring algorithm to analyze single-cell RNA sequencing (scRNA-seq) data, assigning each cell to a specific type based on its gene expression profile. The expected tissue type for the analysis is expected to be submitted by the user.",
  "CNVs (CopyKat)" = "Copy number alterations (CNAs) detected using the CopyKat tool designed to infer genomic CNVs across individual cells using scRNA-seq data. By comparing gene expression profiles, CopyKat identifies regions of the genome that have been amplified or deleted, which can be indicative of genetic abnormalities such as those found in cancer cells.",
  "VAF_RNA" = "VAF_RNA: expressed Variant Allele Fraction. For each individual sceSNV VAF_RNA is calculated as the ratio of the number of variant reads (N_VAR) divided by the total number of reads (N_VAR + N_REF) covering the sceSNV locus (VAF_RNA = N_VAR / (N_VAR + N_REF).",
  "VARreads" = "N_VAR: absolute number of reads bearing the variant nucleotide at the sceSNV locus.",
  "REFreads" = "N_REF: absolute number of reads bearing the reference nucleotide at the sceSNV locus."
)


# function to generate the set of sceSNV 3D plots (with optional slingshot trajectory option)
generate_3d_plot <- function(data, color_column, color_scale, reversescale_option,
                             cell_size, dim_plotting, dim_title, curves, legend_name) {

  plot <- plot_ly(type = 'scatter3d', mode = 'lines+markers') %>%
    add_markers(data = data[data[['Undetected']] == 0, ], 
                x = ~eval(parse(text = paste0(dim_plotting, '_1'))), 
                y = ~eval(parse(text = paste0(dim_plotting, '_2'))), 
                z = ~eval(parse(text = paste0(dim_plotting, '_3'))),
                size = 0.05, opacity = 1.00, name = legend_name,
                marker = list(color = ~get(color_column), colorscale = color_scale,
                              showscale = TRUE, colorbar = list(x = 0.85, y = 0.5, thickness = 20, len = 0.5),  
                              reversescale = reversescale_option,
                              line = list(color = ~get(color_column), width = cell_size))) %>%
    add_markers(data = data[data[['Undetected']] == 1, ], 
                x = ~eval(parse(text = paste0(dim_plotting, '_1'))), 
                y = ~eval(parse(text = paste0(dim_plotting, '_2'))), 
                z = ~eval(parse(text = paste0(dim_plotting, '_3'))),
                size = 0.05, opacity = 1.00, color = 'lightgrey',
                name = 'undetected',
                marker = list(color = 'lightgrey', line = list(width = cell_size)))
  
  if (!args$`disable-slingshot`) {
    plot <- plot %>% add_trace(data = curves, 
                               x = ~eval(parse(text = paste0(dim_plotting, '_1'))), 
                               y = ~eval(parse(text = paste0(dim_plotting, '_2'))), 
                               z = ~eval(parse(text = paste0(dim_plotting, '_3'))),
                               split = ~Lineage, mode = 'lines')
  }
  
  if (args$`disable-3d-axis`) {
    plot <- plot %>% layout(scene = list(
        xaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
        yaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
        zaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)
      )
    )  
  } else {
    plot <- plot %>% layout(scene = list(
        xaxis = list(title = paste0(dim_title, '_1')),
        yaxis = list(title = paste0(dim_title, '_2')),
        zaxis = list(title = paste0(dim_title, '_3'))
      )
    ) 
  }
  
  return(plot)
}

plot1 <- generate_3d_plot(df.3dplot, 'SNV.N', color.scale3, reversescale.option, 
                          cell.size, dim.plotting, dim.title, curves, 'expressed sceSNV loci')
plot2 <- generate_3d_plot(df.3dplot, 'SNVCount', color.scale3, reversescale.option,
                          cell.size, dim.plotting, dim.title, curves, 'expressed sceSNV loci')
plot3 <- generate_3d_plot(df.3dplot, 'RefCount', color.scale3, reversescale.option,
                          cell.size, dim.plotting, dim.title, curves, 'expressed sceSNV loci')
plot4 <- generate_3d_plot(df.3dplot, 'TotalVAF', color.scale3, reversescale.option,
                          cell.size, dim.plotting, dim.title, curves, 'expressed sceSNV loci')
plot5 <- generate_3d_plot(df.3dplot, 'MeanVAF', color.scale3, reversescale.option,
                          cell.size, dim.plotting, dim.title, curves, 'expressed sceSNV loci')
plot6 <- generate_3d_plot(df.3dplot, 'MedianVAF', color.scale3, reversescale.option,
                          cell.size, dim.plotting, dim.title, curves, 'expressed sceSNV loci')


all_plots <- list()   # empty list to append 3D plots into!
all_plots <- append(all_plots, list(plot1, plot2, plot3, plot4, plot5, plot6))

if (args$`save-each-plot`) {
  # save above plots as individual HTML files with titles
  plot3d_names <- c("N_sceSNVs", "N_VARreads", "N_REFreads", "Total_VAF_RNA", "Mean_VAF_RNA", "Median_VAF_RNA")
  plots <- list(plot1, plot2, plot3, plot4, plot5, plot6)

  for (i in 1:length(plots)) {
    plot3d_with_title <- plots[[i]] %>% layout(title = plot3d_names[i], title = list(font = 'black'), margin = list(t = 50))
    saveWidget(as_widget(plot3d_with_title), file = file.path(dir.name, 'Figures_Individual_Plots_HTML', paste0(plot3d_names[i], '.html')), selfcontained = F, libdir = 'lib')
  }
}


### ==================================================== CELL TYPES PLOT ========================================================== ###

# Plot scType with slingshot trajectories
if (args$`enable-sctype`) {
  f <- plot_ly(type = 'scatter3d', mode = 'lines+markers')
  f <- f %>% add_trace(
    data = df.3dplot,
    x = ~eval(parse(text = paste0(dim.plotting, '_1'))),
    y = ~eval(parse(text = paste0(dim.plotting, '_2'))),
    z = ~eval(parse(text = paste0(dim.plotting, '_3'))),
    size = 0.05, opacity = 0.50, mode = 'markers', color = ~customclassif
  )
  
  if (!args$`disable-slingshot`) {
    f <- f %>% add_trace(
      data = curves,
      x = ~eval(parse(text = paste0(dim.plotting, '_1'))),
      y = ~eval(parse(text = paste0(dim.plotting, '_2'))),
      z = ~eval(parse(text = paste0(dim.plotting, '_3'))),
      mode = 'lines', color = ~factor(Lineage), line = list(width = 4)
    )
  }
  
  if (args$`disable-3d-axis`) {
    f <- f %>% layout(scene = list(
      xaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
      yaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
      zaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)
    )
  )
  } else {
    f <- f %>% layout(
      title = '',
      scene = list(
        xaxis = list(title = paste0(dim.title, '_1')),
        yaxis = list(title = paste0(dim.title, '_2')),
        zaxis = list(title = paste0(dim.title, '_3'))
      )
    )
  }
  
  all_plots <- append(all_plots, list(f))

  if (args$`save-each-plot`) {
    # save the Cell types plot as separate html 
    celltype_with_title <- f %>% layout(title = 'Cell types (scType)', title = list(font = 'black'), margin = list(t = 50))
    saveWidget(as_widget(celltype_with_title), file = file.path(dir.name, 'Figures_Individual_Plots_HTML', 'Cell_types_scType.html'), selfcontained = F, libdir = 'lib')
  }
  
  s <- slingshot(as.SingleCellExperiment(srt), clusterLabels = 'seurat_clusters', reducedDim = toupper(dimensionality.reduction))
  
  lineage_counter <- 1
  for (i in grep('slingPseudotime', colnames(s@colData))) {
    lineage.i <- as.integer(sub('slingPseudotime_', '', colnames(s@colData)[i]))
    curve <- curves[curves$Lineage == lineage.i, ]
    
    f <- plot_ly(type = 'scatter3d', mode = 'lines+markers')
    f <- f %>% add_trace(
      data = df.3dplot,
      x = ~eval(parse(text = paste0(dim.plotting, '_1'))),
      y = ~eval(parse(text = paste0(dim.plotting, '_2'))),
      z = ~eval(parse(text = paste0(dim.plotting, '_3'))),
      size = 0.05, opacity = 0.50, mode = 'markers',
      marker = list(
        color = ~s@colData[, i],
        colorscale = 'YlOrRd', reversescale = TRUE, colorbar = list(x = 0),
        line = list(color = ~SNV.N, width = cell.size)
      )
    )
    
    if (args$`disable-3d-axis`) {
      f <- f %>% layout(scene = list(
        xaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
        yaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
        zaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)
      )
    )
  } else {
      f <- f %>% layout(
        title = paste0('slingshot lineage ', lineage_counter),
        scene = list(
          xaxis = list(title = paste0(dim.title, '_1')),
          yaxis = list(title = paste0(dim.title, '_2')),
          zaxis = list(title = paste0(dim.title, '_3'))
      )
    )
  }
    
    if (!args$`disable-slingshot`) {
      f <- f %>% add_trace(
        data = curve,
        x = ~eval(parse(text = paste0(dim.plotting, '_1'))),
        y = ~eval(parse(text = paste0(dim.plotting, '_2'))),
        z = ~eval(parse(text = paste0(dim.plotting, '_3'))),
        mode = 'lines', showlegend = TRUE,
        line = list(width = 4, color = 'black')
      )
    }

    #all_plots <- append(all_plots, list(f))     # optional to save the individual slingshot lineage plots
    lineage_counter <- lineage_counter + 1
  }  
}



### ==================================================== COPY NUMBER VARIATIONS PLOT ========================================================== ###

if (args$`enable-copykat`) {
  # function to generate copykat CNA plot
  generate_copykat_plot <- function(srt, dim_plotting, dim_title, reduction) {
    sce <- as.SingleCellExperiment(srt)
    sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = toupper(reduction))
    curves <- slingCurves(sce, as.df = TRUE)
    colnames(curves)[1:3] <- c(paste0(dim_plotting, "_1"), paste0(dim_plotting, "_2"), paste0(dim_plotting, "_3")) 
    df.srt <- as.data.frame(srt[[c('seurat_clusters', 'karyotype')]])
    df.embed <- as.data.frame(Embeddings(srt, reduction = reduction))
    df.3dplot <- merge(df.srt, df.embed, by = 0)
    colnames(df.3dplot) <- c('filler', 'seurat_clusters', 'karyotype', 
                              paste0(dim_plotting, "_1"),
                              paste0(dim_plotting, "_2"),
                              paste0(dim_plotting, "_3"))
    rownames(df.3dplot) <- df.3dplot$Row.names
    df.3dplot <- df.3dplot[, 2:ncol(df.3dplot)]
    
    f <- plot_ly(type = 'scatter3d', mode = 'lines+markers')
    f <- f %>% add_markers(data = df.3dplot, x = ~get(paste0(dim_plotting, "_1")),
                                            y = ~get(paste0(dim_plotting, "_2")), 
                                            z = ~get(paste0(dim_plotting, "_3")),
                                            size = 0.05, opacity = 1.00, color = ~karyotype)

    f <- f %>% add_trace(data = curves, x = ~get(paste0(dim_plotting, "_1")),
                                        y = ~get(paste0(dim_plotting, "_2")), 
                                        z = ~get(paste0(dim_plotting, "_3")),
                                        split = ~Lineage, mode = 'lines')

    if (args$`disable-3d-axis`) {
      f <- f %>% layout(scene = list(
        xaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
        yaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE),
        zaxis = list(title = '', zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)
      ))  
    } else {
      f <- f %>% layout(scene = list(
        xaxis = list(title = paste0(dim_title, '_1')),
        yaxis = list(title = paste0(dim_title, '_2')),
        zaxis = list(title = paste0(dim_title, '_3'))
      )) 
    }
    
    return(f)
  }

  copykat_plot <- generate_copykat_plot(srt, dim.plotting, dim.title, dimensionality.reduction)
  all_plots <- append(all_plots, list(copykat_plot))

  if (args$`save-each-plot`) {
    # save copykat plot as separate html
    copykat_with_title <- copykat_plot %>% layout(title = 'CopyKat (CNVs)', title = list(font = 'black'), margin = list(t = 50))
    saveWidget(as_widget(copykat_with_title), file = file.path(dir.name, 'Figures_Individual_Plots_HTML', 'CNVs_CopyKat.html'), selfcontained = F, libdir = 'lib')
  }
}



### ==================================================== SNV SUMMARY TEXT FILE  ========================================================== ###

# SNV summary .txt file
if (args$`enable-sctype`){
  colnames(df.3dplot)[1:10] <- c('N_VARreadCounts', 'N_REFreadCounts', 'TotalVAF',
                                 'MedianSNVsVAF', 'MeanSNVsVAF', 'N_SNV', 'seurat_clusters',
                                 'customclassif','HasSNV','Undetected')
} else {
  colnames(df.3dplot)[1:9] <- c('N_VARreadCounts', 'N_REFreadCounts', 'TotalVAF',
                                'MedianSNVsVAF', 'MeanSNVsVAF', 'N_SNV',
                                'seurat_clusters','HasSNV','Undetected')                      
}
df.3dplot <- df.3dplot[, c(6, 1, 5, 4, 3, 7:ncol(df.3dplot))]
write.table(df.3dplot, paste0(dir.name, sample.name, '-summary.txt'), sep='\t', col.names=NA)



### ==================================================== INDIVIDUAL scSNV PLOTS ========================================================== ###

# make plot that contains VARReads, RefReads, dimensional reduction plots, Slingshot
if (!args$`disable-ind-plots`) {
  # Colors for VAF
  pal <- c('#EBEBEB', '#85C1E9', '#E74C3C', '#B03A2E', '#641E16')
  df.snv <- snv
  df.snv <- df.snv[c('CHROM', 'POS', 'REF', 'ALT', 'ReadGroup', 'SNVCount', 'RefCount', 'VAF')]
  snvs <- unique(df.snv[c('CHROM', 'POS', 'REF', 'ALT')])

  # defining the dropdown options for SNV
  snv_options <- paste(snvs$CHROM, snvs$POS, snvs$REF, snvs$ALT, sep = ":")

  # function to generate plots for a selected SNV
  generate_snv_plots <- function(selected_snv, title_color = 'blue') {
    selected_snv <- unlist(strsplit(selected_snv, ":"))
    df.subset <- df.snv[df.snv$CHROM == selected_snv[1] &
                        df.snv$POS == as.numeric(selected_snv[2]) &
                        df.snv$REF == selected_snv[3] &
                        df.snv$ALT == selected_snv[4], ]

    vaf <- df.subset$VAF[match(colnames(srt), df.subset$ReadGroup)]
    snv.reads <- df.subset$SNVCount[match(colnames(srt), df.subset$ReadGroup)]
    ref.reads <- df.subset$RefCount[match(colnames(srt), df.subset$ReadGroup)]

    y <- data.frame(x=df.dim[, 1], y=df.dim[, 2], z=df.dim[, 3],
                    vaf=vaf, ref.reads=ref.reads, snv.reads=snv.reads)

    plots <- list()

    # generate VAF plot
    y$vaf.label <- 'Undetected'
    y$vaf.label[y$vaf==0] <- '0 VAF, REFreads Only'
    y$vaf.label[0<y$vaf & y$vaf<=0.25] <- '0<VAF<=0.25'
    y$vaf.label[0.25<y$vaf & y$vaf<=0.75] <- '0.25<VAF<=0.75'
    y$vaf.label[0.75<y$vaf & y$vaf<=1.00] <- '0.75<VAF<=1.00'
    y$vaf.label <- factor(y$vaf.label, levels=c('Undetected',
      '0 VAF, REFreads Only', '0<VAF<=0.25', '0.25<VAF<=0.75', '0.75<VAF<=1.00'))

    f_vaf <- plot_ly(type='scatter3d', mode='markers+lines')
    for (j in 1:5) {
      cur.label <- levels(y$vaf.label)[j]
      if (args$`enable-dynamic-cell-size`) {
        f_vaf <- f_vaf %>% add_trace(data=subset(y, vaf.label==cur.label),
                                     x=~x, y=~y, z=~z,
                                     size=~(snv.reads + ref.reads)/max(c(snv.reads, ref.reads), na.rm=TRUE)*10,
                                     type='scatter3d', mode='markers',
                                     marker=list(color=pal[j], line=list(width=0)), name=cur.label)
      } else {
        f_vaf <- f_vaf %>% add_trace(data=subset(y, vaf.label==cur.label),
                                     x=~x, y=~y, z=~z,
                                     size=0.05, type='scatter3d', mode='markers',
                                     marker=list(color=pal[j], line=list(width=0)), name=cur.label)
      }
    }

      f_vaf <- f_vaf %>% layout(title='VAF_RNA', 
                                titlefont = list(color=title_color),
                                scene=list(xaxis=list(title=paste0(dim.title,'_1')),
                                            yaxis=list(title=paste0(dim.title,'_2')),
                                            zaxis=list(title=paste0(dim.title,'_3'))))
    
    plots[['VAF']] <- f_vaf

    # generate VARreads plot
    f_varreads <- plot_ly(type='scatter3d', mode='markers+lines') %>%
      add_trace(data=subset(y, is.na(vaf)==0),
                x=~x, y=~y, z=~z,
                size=ifelse(args$`enable-dynamic-cell-size`, 
                            ~(snv.reads + ref.reads)/max(c(snv.reads, ref.reads), na.rm=TRUE)*10, 0.05),
                type='scatter3d', mode='markers',
                marker=list(reversescale=T, color=~snv.reads,
                            colorscale='YlOrRd', showscale=T, opacity=0.50,
                            line=list(color='#FEE5D9', width=1),
                            colorbar=list(len=0.5, y=0.2)), 
                name='Cells with VARreads') %>%
      add_trace(data=subset(y, is.na(vaf)==1),
                x=~x, y=~y, z=~z,
                size=ifelse(args$`enable-dynamic-cell-size`, 
                            ~(snv.reads + ref.reads)/max(c(snv.reads, ref.reads), na.rm=TRUE)*10, 0.05),
                type='scatter3d', mode='markers',
                marker=list(color='#EBEBEB', line=list(width=0), opacity=0.50), 
                name='Cells without VARreads')

      f_varreads <- f_varreads %>% layout(title='VARreads', 
                                          titlefont = list(color=title_color),
                                          scene=list(xaxis=list(title=paste0(dim.title,'_1')),
                                                      yaxis=list(title=paste0(dim.title,'_2')),
                                                      zaxis=list(title=paste0(dim.title,'_3'))))
  
    plots[['VARreads']] <- f_varreads

    # generate REFreads plot
    f_refreads <- plot_ly(type='scatter3d', mode='markers+lines') %>%
      add_trace(data=subset(y, vaf==0 & ref.reads>0),
                x=~x, y=~y, z=~z,
                size=ifelse(args$`enable-dynamic-cell-size`, 
                            ~(snv.reads + ref.reads)/max(c(snv.reads, ref.reads), na.rm=TRUE)*10, 0.05),
                type='scatter3d', mode='markers',
                marker=list(reversescale=T, color=~ref.reads,
                            colorscale='Blues', showscale=T, opacity=0.50,
                            line=list(color='#EFF3FF', width=1),
                            colorbar=list(len=0.5, y=0.2)),
                name='Cells with REFreads') %>%
      add_trace(data=subset(y, (vaf==0 & ref.reads==0)|(is.na(vaf)==1)|(vaf>0)),
                x=~x, y=~y, z=~z,
                size=ifelse(args$`enable-dynamic-cell-size`, 
                            ~(snv.reads + ref.reads)/max(c(snv.reads, ref.reads), na.rm=TRUE)*10, 0.05),
                type='scatter3d', mode='markers',
                marker=list(color='#EBEBEB', line=list(width=0), opacity=0.50), 
                name='Cells without REFreads')

      f_refreads <- f_refreads %>% layout(title='REFreads', 
                                          titlefont = list(color=title_color),
                                          scene=list(xaxis=list(title=paste0(dim.title,'_1')),
                                                     yaxis=list(title=paste0(dim.title,'_2')),
                                                     zaxis=list(title=paste0(dim.title,'_3'))))
  
    plots[['REFreads']] <- f_refreads

    
    if (!args$`disable-slingshot`) {
      for (plot_type in names(plots)) {
        plots[[plot_type]] <- plots[[plot_type]] %>% add_trace(data=curves, 
                                                               x=~eval(parse(text = paste0(dim.plotting,'_1'))), 
                                                               y=~eval(parse(text = paste0(dim.plotting,'_2'))), 
                                                               z=~eval(parse(text = paste0(dim.plotting,'_3'))),
                                                               mode='lines', type='scatter3d', split=~Lineage)
      }
    }
    plots
  }

  # generating all plots and saving their JSON content into a list
  plots_json <- lapply(snv_options, function(snv) {
    plots <- generate_snv_plots(snv, title_color='blue')
    list(
      VAF = list(id = paste0("plot_VAF_", gsub(":", "_", snv)), json = plotly_json(plots[['VAF']])),
      VARreads = list(id = paste0("plot_VARreads_", gsub(":", "_", snv)), json = plotly_json(plots[['VARreads']])),
      REFreads = list(id = paste0("plot_REFreads_", gsub(":", "_", snv)), json = plotly_json(plots[['REFreads']]))
    )
  })

  plots_json <- unlist(plots_json, recursive = FALSE)

  if (args$`save-each-plot`) {
    # function to save individual sceSNV plots seprately
    save_snv_plot <- function(plot, snv, plot_type) {
      snv_clean <- gsub(":", "_", snv)
      file_name <- paste0(plot_type, "_", snv_clean, ".html")
      folder <- switch(plot_type, "VAF" = "VAF", "VARreads" = "VARreads", "REFreads" = "REFreads")
      file_path <- file.path(dir.name, 'Figures_Individual_Plots_HTML', 'Individual_sceSNVs', folder, file_name)

      if (!args$`disable-title`) {
        snv_parts <- unlist(strsplit(snv, ":"))
        snv_title_format <- paste0(snv_parts[1], ":", snv_parts[2], " ", snv_parts[3], ">", snv_parts[4])
        indsnv_title <- paste(plot_type, "<br>", snv_title_format)
        plot <- plot %>% layout(title = list(text = indsnv_title, font = list(color = 'black')), margin = list(t = 50))
      }
      
      #saveWidget(as_widget(indsnv_with_title), file = file_path, selfcontained = F, libdir = 'lib')
      saveWidget(as_widget(plot), file = file_path, selfcontained = F, libdir = 'lib')
      }

    for (snv in snv_options) {
      plots <- generate_snv_plots(snv, title_color='black')
      save_snv_plot(plots[['VAF']], snv, "VAF")
      save_snv_plot(plots[['VARreads']], snv, "VARreads")
      save_snv_plot(plots[['REFreads']], snv, "REFreads")
    }
  }

  # HTML code for individual sceSNV plots  
  individual_SNV_html <- paste0('
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Individual sceSNV plot selection</title>
        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
        body {text-align: center; margin-bottom: 50px;}
        .center-container {display: inline-block; text-align: left;}
        .plot-container {margin-top: 20px; display: grid; grid-template-columns: repeat(2, 1fr); gap: 10px; padding: 0 20px;}
        .plot-item {width: 800px; height: 500px; margin-bottom: 50px;}
        .plot-title {
            cursor: pointer;
            color: blue;
            font-weight: bold;
            text-decoration: none;
            display: inline-block;
        }
        </style>
    </head>
    <body>
        <div class="center-container">
        <h2>Individual sceSNV plot selection</h2>
        <label for="snv">Select SNV:</label>
        <select id="snv">
            ', paste0('<option value="', snv_options, '">', snv_options, '</option>', collapse = ''), '
        </select>
        <br><br>
        <button id="updatePlot">Update Plot</button>
        <p id="selectedSNV"></p>
        </div>
        <br><br>
        <div class="plot-container">
        <div id="plotContainer_VAF" class="plot-item">
            <div class="plot-title" data-plot-id="VAF_RNA">VAF</div>
        </div>
        <div id="plotContainer_VARreads" class="plot-item">
            <div class="plot-title" data-plot-id="VARreads">VARreads</div>
        </div>
        <div id="plotContainer_REFreads" class="plot-item">
            <div class="plot-title" data-plot-id="REFreads">REFreads</div>
        </div>
        </div>
        <script>
        const plots = {
            ', paste0(sapply(plots_json, function(p) {
            paste0('"', p$id, '": ', p$json)
            }), collapse = ", "), '
        };

        function updatePlot() {
            var snv = $("#snv").val();
            var plotIds = ["VAF", "VARreads", "REFreads"];
            plotIds.forEach(function(plotType) {
            var plotId = "plot_" + plotType + "_" + snv.replace(/:/g, "_");
            var plotData = plots[plotId];
            if (plotData) {
                Plotly.newPlot("plotContainer_" + plotType, plotData.data, plotData.layout, {
                responsive: true,
                autosize: true,
                width: document.getElementById("plotContainer_" + plotType).offsetWidth,
                height: document.getElementById("plotContainer_" + plotType).offsetHeight
                });
            }
            });
            $("#selectedSNV").text("Selected SNV: " + snv);
        }

        $(document).ready(function(){
            $("#updatePlot").click(updatePlot);
            $("#snv").val($("#snv option:first").val());
            updatePlot();
        });
        </script>
    </body>
    </html>'
    )
}


# the link HTML to show/hide the histograms
link_html <- '<br><button id="toggleButton" onclick="toggleHistograms()" style="font-size: 20px; padding: 10px 20px; margin: 40px 0;">View Histograms</button><br>'

# JavaScript function to toggle the visibility of the histograms
toggle_histogram_script <- '
    <script>
    function toggleHistograms() {
        var histograms = document.getElementById("histograms");
        var toggleButton = document.getElementById("toggleButton");
        if (histograms.style.display === "none") {
        histograms.style.display = "grid";
        toggleButton.textContent = "Hide Histograms";
        } else {
        histograms.style.display = "none";
        toggleButton.textContent = "View Histograms";
        }
    }
    </script>
'

# JavaScript code for plot title pop-ups
title_popup_script <- '
    <script>
    document.addEventListener("DOMContentLoaded", function() {
        document.querySelectorAll(".plot-title").forEach(function(element) {
        element.addEventListener("click", function() {
            var plotId = element.getAttribute("data-plot-id");
            var descriptions = {
            "N_sceSNVs": "Absolute number of sceSNVs per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
            "N_VARreads": "Sum of the absolute number of reads bearing the variant nucleotide across all sceSNVs in the submitted list per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
            "N_REFreads": "Sum of the absolute number of reads bearing the reference nucleotide across all sceSNVs in the submitted list per cell. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
            "Total_VAF_RNA": "Total VAF_RNA is calculated by dividing the sum of N_VAR counts across the sceSNVs (loci) by the total reads covering the sceSNV loci (N_VAR + N_REF), and is intended to provide an assessment of SNV expression magnitude across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
            "Mean_VAF_RNA": "Mean VAF_RNA is calculated as the mean of the VAF_RNA across all individual sceSNVs and is intended to provide an assessment of the variability of the SNV expression levels across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
            "Median_VAF_RNA": "Median VAF_RNA is calculated as the median of the VAF_RNA across all individual sceSNVs and is intended to provide an assessment of the central tendency of the SNV expression levels across cell populations. Only sceSNVs loci from the submitted list are considered, after applying the user-selected thresholds.",
            "Cell types (scType)": "Classification of individual cells into distinct types using the scType tool. scType employs predefined marker gene sets and a scoring algorithm to analyze single-cell RNA sequencing (scRNA-seq) data, assigning each cell to a specific type based on its gene expression profile. The expected tissue type for the analysis is expected to be submitted by the user.",
            "CNVs (CopyKat)": "Copy number alterations (CNAs) detected using the CopyKat tool designed to infer genomic CNVs across individual cells using scRNA-seq data. By comparing gene expression profiles, CopyKat identifies regions of the genome that have been amplified or deleted, which can be indicative of genetic abnormalities such as those found in cancer cells.",
            "VAF_RNA": "VAF_RNA: expressed Variant Allele Fraction. For each individual sceSNV VAF_RNA is calculated as the ratio of the number of variant reads (N_VAR) divided by the total number of reads (N_VAR + N_REF) covering the sceSNV locus (VAF_RNA = N_VAR / (N_VAR + N_REF).",
            "VARreads": "N_VAR: absolute number of reads bearing the variant nucleotide at the sceSNV locus.",  
            "REFreads": "N_REF: absolute number of reads bearing the reference nucleotide at the sceSNV locus."
            };

            var existingTooltip = document.querySelector(`.tooltip-box[data-plot-id="${plotId}"]`);
            if (existingTooltip) {
            document.body.removeChild(existingTooltip);
            return;
            }

            var tooltip = document.createElement("div");
            tooltip.className = "tooltip-box";
            tooltip.setAttribute("data-plot-id", plotId);
            tooltip.innerHTML = descriptions[plotId];

            document.body.appendChild(tooltip);

            var rect = element.getBoundingClientRect();
            var tooltipWidth = tooltip.offsetWidth;
            var tooltipHeight = tooltip.offsetHeight;
            var pageWidth = window.innerWidth;
            var pageHeight = window.innerHeight;

            // Adjust tooltip positioning to ensure it remains within page boundaries
            if (rect.right + tooltipWidth > pageWidth) {
            tooltip.style.left = (rect.left - tooltipWidth - 10) + "px";
            } else {
            tooltip.style.left = (rect.right + 10) + "px";
            }

            if (rect.bottom + tooltipHeight > pageHeight) {
            tooltip.style.top = (rect.top - tooltipHeight - 10) + "px";
            } else {
            tooltip.style.top = (rect.bottom + window.scrollY) + "px";
            }

            window.addEventListener("scroll", function() {
            var updatedRect = element.getBoundingClientRect();
            tooltip.style.left = (updatedRect.left + window.scrollX) + "px";
            tooltip.style.top = (updatedRect.bottom + window.scrollY) + "px";
            }, { passive: true });
        });
        });
    });
    </script>
'

# CSS for tooltip and plot titles
tooltip_css <- '
    <style>
    .tooltip-box {
        position: absolute;
        background-color: white;
        border: 1px solid black;
        padding: 10px;
        z-index: 1000;
        max-width: 300px;
        box-shadow: 0px 0px 10px rgba(0,0,0,0.5);
        font-size: 16px; /* Adjust text size */
    }
    .plot-title {
        cursor: pointer;
        color: blue;
        font-weight: bold;
        text-decoration: none;
        display: inline-block;
    }
    </style>
'


### ============================================== COMBINING ALL PLOTS ================================================= ###

#assigning classes to 3D plots
plot_divs <- lapply(seq_along(all_plots), function(i) {
  plot <- all_plots[[i]]
  plot_id <- names(plot_descriptions)[i]
  div(class = "grid-item plot3d", 
      div(class = "plot-title", `data-plot-id` = plot_id, plot_id),  # plot_id for title
      as_widget(plot)
  )
})

# CSS layout style
grid_html <- tags$html(
  tags$head(
    tags$style(HTML("
      .grid-container {
        display: grid;
        grid-template-columns: repeat(2, 1fr);  /* 1st arg. specifies no of columns */
        gap: 20px 5px;   /* adjust the gap: 1st value for row spacing, 2nd for column spacing */
        justify-content: center; 
      }
      .grid-item {
        padding: 0; 
      }
      .grid-item.plot3d {
        width: 100%;   /* width for snv plots */
        height: 100%;  /* height for snv plots */
        text-align: center; /* center align titles */
      }
      .grid-item.histogram {
        width: 100%;   /* width for histograms */
        height: 100%;  /* height for histograms */
      }
      #histograms {
        display: none;
        grid-template-columns: repeat(2, 1fr); /* 2-column grid for histograms */
        gap: 10px;
        justify-content: center;
      }
      .button-container {
        width: 100%;         
        display: flex;    /* flexbox is for centering */
        justify-content: center;  
        margin: 20px 0;     
      }
      #toggleButton {
        display: inline-block; /* inline-block for better control */
        font-size: 20px;
        padding: 10px 20px;
      }
      body {
        font-size: 20px; 
      }
    ")),
    HTML(tooltip_css)
  ),
  tags$body(
    div(class = "grid-container", plot_divs),
    div(class = "button-container",
      HTML(link_html)
    ),
    div(id = "histograms", 
      lapply(histograms, function(plot) {
        div(class = "grid-item histogram", as_widget(plot))
      })
    ),
    if (!args$`disable-ind-plots`) {
      HTML(individual_SNV_html)
    },
    HTML(toggle_histogram_script),
    HTML(title_popup_script)
  )
)

# saving the final html file with 3D plots, histograms, and individual sceSNV plots
output_file <- paste0(dir.name, 'Exploratory_Combined_Plots.html')
save_html(grid_html, file = output_file)
