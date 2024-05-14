# May 14, 2024
suppressWarnings(library('optparse'))
suppressWarnings(library('stringr'))

script.desc <-
'This script calculates and plots basic statistics, 2-dimensional and
3-dimensional PCAs of SNV reads, reference reads, and VAF for each cell. The
user must supply the RDS of a Seurat object (-r) which can contain multiple
integrated datasets. The user must also supply an SNV file (-t) that contains
at least one SNV.'

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
parser <- add_option(parser, c('-x', '--th-vars'),
                     type='integer', default=0,
                     help='Threshold for number of variants. Default=0.')
parser <- add_option(parser, c('-y', '--th-reads'),
                     type='integer', default=0,
                     help='Threshold for number of variant reads. Default=0.
                            An SNV needs more than 0 counts to be considered an SNV.')
parser <- add_option(parser, c('-c', '--enable-title'),
                     type='logical', default=T, action='store_true',
                     help='Enable title. Default=T')
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
                     help='Enable sctype to run. Default=F.')
parser <- add_option(parser, c('-j', '--tissue-type'),
                     type='character',
                     help='tissue type for scType; options include:
                     Immune system, Pancreas, Liver, Eye, Kidney, Brain,
                     Lung, Adrenal, Heart, Intestine, Muscle, Placenta,
                     Spleen, Stomach, Thymus')
parser <- add_option(parser, c('-k', '--color-scale'),
                     type='character',
                     help='if you would like to change the default color settings with
                     these options, you may use Blues, Reds, YlOrRd, YlGnBu, plasma. Default=\'Blues\'')
parser <- add_option(parser, c('-p', '--enable-cell-border'),
                     type='logical', default=F, action='store_true',
                     help='Enable cell border. Default=F')
parser <- add_option(parser, c('-q', '--enable-dynamic-cell-size'),
                     type='logical', default=F, action='store_true',
               help='Enable cell size to depend on number of reads. Default=F')
args <- parse_args(parser)

error.msg <- NULL

# Check if the required argument (-r) is passed
if (is.null(args$`rds-file`) & is.null(args$`countsmatrix-file`))
  error.msg <- paste(error.msg, '- Seurat RDS object (-r) or STARsolo output directory for counts matrix is required.',
                     sep='\n')
if (is.null(args$`snv-file`))
  error.msg <- paste(error.msg, '- scReadCounts file (-t) is required.',
                     sep='\n')


if (args$`enable-sctype`) {
  if (is.null(args$`tissue-type`)) {
    error.msg <- paste(error.msg, '- tissue type is required')
  } 
  else {
    tissue.type <- args$`tissue-type`
    tissue.type <- str_to_title(tolower(tissue.type))
    if (!(tissue.type %in% list('Immune system', 'Pancreas', 'Liver', 'Eye', 
          'Kidney', 'Brain', 'Lung', 'Adrenal', 'Heart', 'Intestine', 'Muscle',
          'Placenta', 'Spleen','Stomach','Thymus'))) {
      error.msg <- paste(error.msg, '- the tissue type you have provided does not seem to be one of the tissue types that scType accepts',
                         sep='\n')
    }
  }
}

if (!is.null(args$`color-scale`)){
  color.scale <- args$`color-scale`
  if (!(color.scale %in% list('Blues', 'Reds', 'YlOrRd', 'YlGnBu', 'plasma'))) {
    error.msg <- paste(error.msg, '- the color scale you provided is either not in the acceptable form or is not one this program uses')
  }
}

cell.size = 0
if (args$`enable-cell-border`) {
  cell.size=3
}

# Check if there are any errors
if (!is.null(error.msg)) {
  print_help(parser)
  stop(error.msg)}

######
# Application Logic
######

library('Seurat')
library('ggplot2')
library('dplyr')
library('openxlsx')
library('HGNChelper')

# load files
snv.file <- args$`snv-file`
th.vars <- args$`th-vars`
th.reads <- args$`th-reads`

sample.name <- sub('^(.*/)*(.*)\\.([a-z0-9]+)$', '\\2', ignore.case=T, snv.file)
param.name <- paste0('_pca_', th.vars, 'r', th.reads)

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
    #srt <- NormalizeData(srt)
    #srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
    #srt <- ScaleData(srt, vars.to.regress = "percent.mt")
}


######
# Color Scales
######

if (is.null(args$`color-scale`)){
color.scale1 <- c('#cfe2f3','tomato')      # for histograms and combined plots
color.scale2 <- c('#cfe2f3','dodgerblue2') # for histograms and combined plots
color.scale3 <- 'Blues'                    # color scales for 2D and 3D plots; simple alternatives include: Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu YlOrBr, YlOrRd
histogram.scale1<-'tomato'
histogram.scale2<-'dodgerblue2'
color.undetected <- 'tomato'               # color for cells where the SNV/s is/are undetected
reversescale.option<-T
} else {
color.scale1 <- color.scale                # for histograms and combined plots
color.scale2 <- color.scale                # for histograms and combined plots
color.scale3 <- color.scale                # color scales for 2D and 3D plots; simple alternatives include: Blues, BuGn, BuPu, GnBu, Greens, Greys, Oranges, OrRd, PuBu, PuBuGn, PuRd, Purples, RdPu, Reds, YlGn, YlGnBu YlOrBr, YlOrRd
color.undetected <- 'tomato'               # color for cells where the SNV/s is/are undetected
histogram.scale1<-'tomato'
histogram.scale2<-'dodgerblue2'
scale.list<-list('Blues','Reds','YlOrRd','YlGnBu','plasma')
reversescale.options <- list(T,F,T,T,F)
reversescale.option<-reversescale.options[[which(scale.list==color.scale)]]
}

srt <- RunPCA(srt)
srt <- RunUMAP(srt, dims=1:20, n.components=3)
srt <- FindNeighbors(srt, dims = 1:10)
srt <- FindClusters(srt, resolution = 0.5)
snv <- read.table(snv.file, sep='\t', header=T)
snv.modify<- snv %>% filter(SNVCount>0)
num.snvs <- nrow(snv.modify)

# if scType option is selected, create cell type labels
if (args$`enable-sctype`){
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  tissue = tissue.type
  gs_list = gene_sets_prepare(db_, tissue)
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = srt[[DefaultAssay(srt)]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(srt@meta.data$seurat_clusters), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(srt@meta.data[srt@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(srt@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  srt@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    srt@meta.data$customclassif[srt@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
}

# calculate number of snv and ref reads per cell
snv.reads <- aggregate(SNVCount~ReadGroup, data=snv, sum)
if (max(snv.reads[,'SNVCount'])==0){
  stop('Even though you have provided a tab-separated file containing SNV information, there aren\'t any SNV read counts in it.')
}
rownames(snv.reads) <- snv.reads$ReadGroup
ref.reads <- aggregate(RefCount~ReadGroup, data=snv, sum)
rownames(ref.reads) <- snv.reads$ReadGroup
if (max(snv.reads[, 'SNVCount']) < th.reads) {
  stop('There aren\'t any SNVs that satisfy the read threshold.')
}

# calculate number of detected snvs per cell
snv.n <- aggregate(SNVCount~ReadGroup, data=snv,
  subset=SNVCount>=th.reads & SNVCount>0, length)
rownames(snv.n) <- snv.n$ReadGroup
colnames(snv.n)[2] <- 'SNV.N'

if (max(snv.n[, 'SNV.N']) < th.vars) {
  stop('There aren\'t any SNVs that satisfy the variant threshold.')
}

snv <- merge(snv, snv.n, 'ReadGroup', all.x=T)

suppressWarnings(dir.create(paste0(sample.name, param.name)))
dir.name <- paste0(sample.name, param.name, '/')

# calculate median vaf per cell
vaf.median <- aggregate(VAF~ReadGroup, data=snv,
  subset=SNVCount>=th.reads & SNV.N>=th.vars, median)
rownames(vaf.median) <- vaf.median$ReadGroup
colnames(vaf.median)[2] <- 'Median.VAF'

# calculate mean vaf per cell
vaf.mean <- aggregate(VAF~ReadGroup, data=snv,
  subset=SNVCount>=th.reads & SNV.N>=th.vars, mean)
rownames(vaf.mean) <- vaf.mean$ReadGroup
colnames(vaf.mean)[2] <- 'Mean.VAF'

# set meta-data
srt[['SNVCount']] <- 0
srt[['SNVCount']][rownames(snv.reads), 1] <- snv.reads[, 'SNVCount']

srt[['RefCount']] <- 0
srt[['RefCount']][rownames(ref.reads), 1] <- ref.reads[, 'RefCount']

srt[['TotalVAF']] <- 0
srt[['TotalVAF']] <- srt[['SNVCount']]/(srt[['SNVCount']]+srt[['RefCount']])
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


# generate dataframes for plots
srt[['HasSNV']] <- sapply(srt[['SNVCount']][, 1],
                          function (x) if (x>=th.reads) 1 else 0)
srt[['Undetected']] <- sapply(rowSums(srt[[c('TotalVAF','HasSNV')]]),
                               function (x) if (x==0 | is.na(x)==1) 1 else 0)
if (args$`enable-sctype`){
  df.snv <- as.data.frame(srt[[c('SNVCount', 'RefCount', 'TotalVAF',
  'MedianVAF', 'MeanVAF', 'SNV.N','seurat_clusters','customclassif','HasSNV','Undetected')]])
  df.pca <- as.data.frame(Embeddings(srt, reduction='pca'))
  df.3dplot <- merge(df.snv, df.pca, by=0)
  colnames(df.3dplot)[12:14] <- c('PC_1', 'PC_2', 'PC_3')
} else {
  df.snv <- as.data.frame(srt[[c('SNVCount', 'RefCount', 'TotalVAF',
  'MedianVAF', 'MeanVAF', 'SNV.N','seurat_clusters','HasSNV','Undetected')]])
  df.pca <- as.data.frame(Embeddings(srt, reduction='pca'))
  df.3dplot <- merge(df.snv, df.pca, by=0)
  colnames(df.3dplot)[11:13] <- c('PC_1', 'PC_2', 'PC_3')
}
rownames(df.3dplot) <- df.3dplot$Row.names
df.3dplot <- df.3dplot[, 2:ncol(df.3dplot)]

# generate interactive 2D plots
library('plotly')
library('htmlwidgets')

# generate slingshot trajectory
library('slingshot')
sce <- as.SingleCellExperiment(srt)
sce <- slingshot(sce, clusterLabels='seurat_clusters', reducedDim='PCA')

curves <- slingCurves(sce, as.df=T)
colnames(curves)[1:3] <- c('PC_1', 'PC_2', 'PC_3')

if (num.snvs>1) {
  # generate histograms
  p <- ggplot(df.snv, aes(x=TotalVAF))+
       geom_histogram(fill=histogram.scale2,boundary=0,alpha=0.6)+
       xlab('TotalVAF')+ylab('cells')+
       theme(panel.grid.minor=element_blank(), axis.ticks.x=element_blank(),
             axis.ticks.y=element_blank())
  if (args$`enable-title`) {
    p <- p+ggtitle(sample.name)
  }
  p
  ggsave(paste0(dir.name, sample.name, '_hTotalVAF', param.name, '.png'))

  p <- ggplot(df.snv, aes(x=MeanVAF))+
       geom_histogram(fill=histogram.scale2, boundary=0, alpha=0.6)+
       xlab('MeanSNVsVAF')+ylab('cells')+
       theme(panel.grid.minor=element_blank(), axis.ticks.x=element_blank(),
             axis.ticks.y=element_blank())
  if (args$`enable-title`) {
    p <- p+ggtitle(sample.name)
  }
  p
  ggsave(paste0(dir.name, sample.name, '_hMeanSNVsVAF', param.name, '.png'))

  p <- ggplot(df.snv, aes(x=SNVCount))+
       geom_histogram(fill=histogram.scale1, binwidth=1, boundary=0, alpha=0.6)+
       xlab('N_VARreadCounts')+ylab('cells')+
       theme(panel.grid.minor=element_blank(), axis.ticks.x=element_blank(),
             axis.ticks.y=element_blank())
  if (args$`enable-title`) {
    p <- p+ggtitle(sample.name)
  }
  p
  ggsave(paste0(dir.name, sample.name, '_hN_VARreadCounts', param.name, '.png'))

  df.snv_mod <- rbind(df.snv$SNV.N)
  p <- ggplot(df.snv, aes(x=SNV.N))+
       geom_histogram(fill=histogram.scale1, binwidth=1, boundary=0, alpha=0.6)+
       xlab('N_SNV')+ylab('cells')+
       theme(panel.grid.minor=element_blank(), axis.ticks.x=element_blank(),
             axis.ticks.y=element_blank())
  if (args$`enable-title`) {
    p <- p+ggtitle(sample.name)
  }
  p
  ggsave(paste0(dir.name, sample.name, '_hN_SNVs', param.name, '.png'))
}

# generate 3D plots
f <- plot_ly(type='scatter3d', mode='lines+markers',colors=color.undetected)
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==0,], 
                                      x=~PC_1, y=~PC_2, z=~PC_3,
                                      size=0.05, opacity=1.00, name='cells',
                       marker=list(color=~TotalVAF, colorscale=color.scale3,
                                   showscale=T, colorbar=list(x=0),
                                   reversescale=reversescale.option,
                                   line=list(color=~TotalVAF, width=cell.size)))
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==1,], 
                                      x=~PC_1, y=~PC_2, z=~PC_3, 
                                      size=0.05, opacity=1.00, 
                                      color=color.undetected,
                                      name='undetected',
                                      marker=list(line=list(width=cell.size)))
if (!args$`disable-slingshot`) {
  f <- f %>% add_trace(data=curves, x=~PC_1, y=~PC_2, z=~PC_3,
                      split=~Lineage, mode='lines')
}
if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(sample.name,'\n TotalVAF'))
}
if (args$`disable-3d-axis`) {
  f <- f %>% layout(scene=list(xaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F),
                               yaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F),
                               zaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F)))
}
saveWidget(f, paste0(dir.name, sample.name, '_TotalVAF_Slingshot',
                     param.name, '.html'),
           selfcontained=F, libdir='lib')

f <- plot_ly(type='scatter3d', mode='lines+markers',colors=color.undetected)
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==0,], 
                                      x=~PC_1, y=~PC_2, z=~PC_3, 
                                      size=0.05, opacity=1.00, name='cells',
                       marker=list(color=~MedianVAF, colorscale=color.scale3,
                                   showscale=T, colorbar=list(x=0),
                                   reversescale=reversescale.option, 
                                   line=list(color=~MedianVAF,
                                             width=cell.size)))
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==1,],
                       x=~PC_1, y=~PC_2, z=~PC_3, size=0.05,
                       opacity=1.00, color=color.undetected,
                       name='undetected',
                       marker=list(line=list(width=cell.size)))
if (!args$`disable-slingshot`) {
  f <- f %>% add_trace(data=curves, x=~PC_1, y=~PC_2, z=~PC_3,
                     split=~Lineage, mode='lines')
}
if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(sample.name,'\n MedianVAF'))
}
if (args$`disable-3d-axis`) {
  f <- f %>% layout(scene=list(xaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F),
                               yaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F),
                               zaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F)))
}
saveWidget(f, paste0(dir.name, sample.name, '_MedianVAF_Slingshot',
                     param.name, '.html'),
           selfcontained=F, libdir='lib')

f <- plot_ly(type='scatter3d', mode='lines+markers',colors=color.undetected)
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==0,], 
                                      x=~PC_1, y=~PC_2, z=~PC_3,
                                      size=0.05, opacity=1.00, name='cells',
                       marker=list(color=~MeanVAF, colorscale=color.scale3,
                                   showscale=T, colorbar=list(x=0),
                                   reversescale=reversescale.option,
                                   line=list(color=~MeanVAF,
                                             width=cell.size)))
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==1,],
                       x=~PC_1, y=~PC_2, z=~PC_3, size=0.05,
                       opacity=1.00, color=color.undetected,name='undetected')
if (!args$`disable-slingshot`) {
  f <- f %>% add_trace(data=curves, x=~PC_1, y=~PC_2, z=~PC_3,
                     split=~Lineage, mode='lines')
}
if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(sample.name,'\n MeanVAF'))
}
if (args$`disable-3d-axis`) {
  f <- f %>% layout(scene=list(xaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F),
                               yaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F),
                               zaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F)))
}
saveWidget(f, paste0(dir.name, sample.name, '_MeanVAF_Slingshot',
                     param.name, '.html'),
           selfcontained=F, libdir='lib')

f <- plot_ly(type='scatter3d', mode='lines+markers',colors=color.undetected)
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==0,], 
                       x=~PC_1, y=~PC_2, z=~PC_3,
                       size=0.05, opacity=1.00, name='cells',
                       marker=list(color=~SNVCount, colorscale=color.scale3,
                                   showscale=T, colorbar=list(x=0),
                                   reversescale=reversescale.option,
                                   line=list(color=~SNVCount,
                                             width=cell.size)))
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==1,], 
                       x=~PC_1, y=~PC_2, z=~PC_3, size=0.05,
                       opacity=1.00, color=color.undetected,
                       name='undetected')
if (!args$`disable-slingshot`) {
  f <- f %>% add_trace(data=curves, x=~PC_1, y=~PC_2, z=~PC_3,
                       split=~Lineage, mode='lines')
}
if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(sample.name,'\n SNVCount'))
}
if (args$`disable-3d-axis`) {
  f <- f %>% layout(scene=list(xaxis=list(title='', zeroline=F, showline=F,
                                          showticklabels=F, showgrid=F),
                               yaxis=list(title='', zeroline=F, showline=F,
                                          showticklabels=F, showgrid=F),
                               zaxis=list(title='', zeroline=F, showline=F,
                                          showticklabels=F, showgrid=F)))
}
saveWidget(f, paste0(dir.name, sample.name, '_N_VARreadCounts_Slingshot',
                     param.name, '.html'),
           selfcontained=F, libdir='lib')

f <- plot_ly(type='scatter3d', mode='lines+markers',colors=color.undetected)
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==0,], 
                                      x=~PC_1, y=~PC_2, z=~PC_3,
                                      size=0.05, opacity=1.00, name='cells',
                        marker=list(color=~RefCount, colorscale=color.scale3,
                                    showscale=T, colorbar=list(x=0),
                                    reversescale=reversescale.option,
                                    line=list(color=~RefCount,
                                              width=cell.size)))
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==1,],
                       x=~PC_1, y=~PC_2, z=~PC_3, size=0.05,
                       opacity=1.00, color=color.undetected, name='undetected')
if (!args$`disable-slingshot`) {
  f <- f %>% add_trace(data=curves, x=~PC_1, y=~PC_2, z=~PC_3,
                       split=~Lineage, mode='lines')
}
if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(sample.name,'\n RefreadCounts'))
}
if (args$`disable-3d-axis`) {
  f <- f %>% layout(scene=list(xaxis=list(title='', zeroline=F, showline=F,
                                          showticklabels=F, showgrid=F),
                               yaxis=list(title='', zeroline=F, showline=F,
                                          showticklabels=F, showgrid=F),
                               zaxis=list(title='', zeroline=F, showline=F,
                                          showticklabels=F, showgrid=F)))
}
saveWidget(f, paste0(dir.name, sample.name, '_REFreadCounts_Slingshot',
                     param.name, '.html'),
           selfcontained=F, libdir='lib')

f <- plot_ly(type='scatter3d', mode='lines+markers', colors=color.undetected)
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==0,], 
                       x=~PC_1, y=~PC_2, z=~PC_3,
                       size=0.05, opacity=1.00, name='cells',
                       marker=list(color=~SNV.N, colorscale=color.scale3,
                                   showscale=T, colorbar=list(x=0),
                                   reversescale=reversescale.option,
                                   line=list(color=~SNV.N,
                                             width=cell.size)))
f <- f %>% add_markers(data=df.3dplot[df.3dplot[['Undetected']]==1,], 
                       x=~PC_1, y=~PC_2, z=~PC_3, size=0.05,
                       opacity=1.00, color=color.undetected, name='undetected')
if (!args$`disable-slingshot`) {
  f <- f %>% add_trace(data=curves, x=~PC_1, y=~PC_2, z=~PC_3,
                       split=~Lineage, mode='lines')
}
if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(sample.name,'\n N_SNV'))
}
if (args$`disable-3d-axis`) {
  f <- f %>% layout(scene=list(xaxis=list(title='', zeroline=F, showline=F,
                                          showticklabels=F, showgrid=F),
                               yaxis=list(title='', zeroline=F, showline=F,
                                          showticklabels=F, showgrid=F),
                               zaxis=list(title='', zeroline=F, showline=F,
                                          showticklabels=F, showgrid=F)))
}
saveWidget(f, paste0(dir.name, sample.name, '_N_SNV_Slingshot',
                     param.name, '.html'),
           selfcontained=F, libdir='lib')


if (args$`enable-sctype`){
ggplotColors <- function(n = 6, h = c(0, 360) + 15){
if ((diff(h)/360) < 1) h[2] <- h[2] - 360/n
hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)}

color.list <- ggplotColors(n=length(unique(srt@meta.data$seurat_clusters)))

f <- plot_ly(data=df.3dplot, x=~PC_1, y=~PC_2, z=~PC_3,
             color=~customclassif, colors=color.list)
f <- f %>% add_markers(size=0.05, opacity=1.00)
if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(sample.name,'\n scType'))
}
if (args$`disable-3d-axis`) {
  f <- f %>% layout(scene=list(xaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F),
                               yaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F),
                               zaxis=list(title='', zeroline=F, showline=F,
                                            showticklabels=F, showgrid=F)))
}
saveWidget(f, paste0(dir.name, sample.name, '_scType.html'),
           selfcontained=F, libdir='lib')
}

if (args$`enable-sctype`){
  colnames(df.3dplot)[1:10] <- c('N_VARreadCounts', 'N_REFreadCounts',
    'TotalVAF', 'MedianSNVsVAF', 'MeanSNVsVAF',
    'N_SNV','seurat_clusters','customclassif','HasSNV','Undetected')
} else {
  colnames(df.3dplot)[1:9] <- c('N_VARreadCounts', 'N_REFreadCounts',
    'TotalVAF', 'MedianSNVsVAF', 'MeanSNVsVAF',
    'N_SNV','seurat_clusters','HasSNV','Undetected')
}
df.3dplot <- df.3dplot[, c(6, 1, 5, 4, 3, 7:ncol(df.3dplot))]
write.table(df.3dplot, paste0(dir.name, sample.name, '-summary.txt'),
            sep='\t', col.names=NA)






#make plot that contains VARReads, RefReads, PCAs, Slingshot
if (!args$`disable-ind-plots`) {
suppressWarnings(dir.create(
                 paste0(dir.name,sample.name, '_individual_sceSNVs/')))
suppressWarnings(dir.create(
                 paste0(dir.name,sample.name, '_individual_sceSNVs/VAF/')))
suppressWarnings(dir.create(
                 paste0(dir.name,sample.name, '_individual_sceSNVs/VARreads/')))
suppressWarnings(dir.create(
                 paste0(dir.name,sample.name, '_individual_sceSNVs/REFreads/')))
# Colors for VAF
# Undetected, 0 VAF, 0<VAF<=0.25, 0.25<VAF<=0.75,  0.75<VAF
pal <- c('#EBEBEB', '#85C1E9', '#E74C3C', '#B03A2E', '#641E16')
df.snv <- snv
df.snv <- df.snv[c('CHROM', 'POS', 'REF', 'ALT', 'ReadGroup', 'SNVCount',
  'RefCount', 'VAF')]
snvs <- unique(df.snv[c('CHROM', 'POS', 'REF', 'ALT')])

suppressWarnings({
# Iterate through SNVs
for (i in 1:nrow(snvs)) {
  df.subset <- df.snv[df.snv$CHROM == snvs$CHROM[i] &
                      df.snv$POS == snvs$POS[i] &
                      df.snv$REF == snvs$REF[i] &
                      df.snv$ALT == snvs$ALT[i], ]

  # Determine VAF for each SNV
  vaf <- colnames(srt)
  vaf <- NA
  vaf <- df.subset$VAF[match(colnames(srt), df.subset$ReadGroup)]
  snv.reads <- df.subset$SNVCount[match(colnames(srt), df.subset$ReadGroup)]
  ref.reads <- df.subset$RefCount[match(colnames(srt), df.subset$ReadGroup)]

  y <- data.frame(x=df.pca[, 1], y=df.pca[, 2], z=df.pca[, 3],
                  vaf=vaf, ref.reads=ref.reads, snv.reads=snv.reads)

  # set colors
  y$vaf.label <- 'Undetected'
  y$vaf.label[y$vaf==0] <- '0 VAF, REFreads Only'
  y$vaf.label[0<y$vaf & y$vaf<=0.25] <- '0<VAF<=0.25'
  y$vaf.label[0.25<y$vaf & y$vaf<=0.75] <- '0.25<VAF<=0.75'
  y$vaf.label[0.75<y$vaf & y$vaf<=1.00] <- '0.75<VAF<=1.00'
  y$vaf.label <- factor(y$vaf.label, levels=c('Undetected',
    '0 VAF, REFreads Only', '0<VAF<=0.25', '0.25<VAF<=0.75', '0.75<VAF<=1.00'))

  f <- plot_ly(type='scatter3d', mode='markers+lines')
  for (j in 1:5) {
    cur.label <- levels(y$vaf.label)[j]
    if (args$`enable-dynamic-cell-size`) {
      f <- f %>% add_trace(data=subset(y,vaf.label==cur.label),
                           x=~x, y=~y, z=~z,
                           size=~(snv.reads+ref.reads)/max(snv.reads+ref.reads)*10,
                           type='scatter3d', mode='markers',
                           marker=list(color=pal[j],line=list(width=0)), name=cur.label)
    } else {
      f <- f %>% add_trace(data=subset(y,vaf.label==cur.label),
                           x=~x, y=~y, z=~z,
                           size=0.05, type='scatter3d', mode='markers',
                           marker=list(color=pal[j],line=list(width=0)), name=cur.label)
    }
  }
  if (!args$`disable-slingshot`) {
    f <- f %>% add_trace(data=curves, x=~PC_1, y=~PC_2, z=~PC_3,
                         mode='lines', type='scatter3d', split=~Lineage)
  }

  if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(snvs$CHROM[i], ':', snvs$POS[i], ' ',
                                 snvs$REF[i], '>', snvs$ALT[i], 
                                 ' (', sample.name, ')'),
                    scene=list(xaxis=list(title='PC_1'),
                               yaxis=list(title='PC_2'),
                               zaxis=list(title='PC_3')))
  }

  f
  file.name <- paste0(sample.name, '-vaf-', snvs$CHROM[i], '-',
                      snvs$POS[i], '-', snvs$REF[i], '-', snvs$ALT[i], '.html')
  saveWidget(f, paste0(dir.name,sample.name, '_individual_sceSNVs/VAF/', file.name),
             selfcontained=F, libdir='lib')

f <- plot_ly(type='scatter3d', mode='markers+lines')
f <- f %>% add_trace(data=subset(y,is.na(vaf)==0),
                       x=~x, y=~y, z=~z,
                       size=0.05, type='scatter3d', mode='markers',
                       marker=list(reversescale=T, color=~snv.reads,
                       colorscale='YlOrRd', showscale=T, opacity=0.50,
                       line=list(color='#FEE5D9', width=1),
                       colorbar=list(len=0.3, y=0.2, title='VARreads')), 
                       name='Cells with VARreads')
f <- f %>% add_trace(data=subset(y,is.na(vaf)==1),
                       x=~x, y=~y, z=~z,
                       size=0.05, type='scatter3d', mode='markers',
                       marker=list(color='#EBEBEB', line=list(width=0),opacity=0.50), 
                       name='Cells without VARreads')
if (!args$`disable-slingshot`) {
    f <- f %>% add_trace(data=curves, x=~PC_1, y=~PC_2, z=~PC_3,
                         mode='lines', type='scatter3d', split=~Lineage)
}
if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(snvs$CHROM[i], ':', snvs$POS[i], ' ',
                                 snvs$REF[i], '>', snvs$ALT[i], 
                                 ' (', sample.name, ')'),
                    scene=list(xaxis=list(title='PC_1'),
                               yaxis=list(title='PC_2'),
                               zaxis=list(title='PC_3')))
}
f
file.name <- paste0(sample.name, '-VARreads-', snvs$CHROM[i], '-',
                    snvs$POS[i], '-', snvs$REF[i], '-', snvs$ALT[i], '.html')
saveWidget(f, paste0(dir.name,sample.name, '_individual_sceSNVs/VARreads/', file.name),
           selfcontained=F, libdir='lib')


f <- plot_ly(type='scatter3d', mode='markers+lines')
f <- f %>% add_trace(data=subset(y, vaf==0 & ref.reads>0),
                       x=~x, y=~y, z=~z,
                       size=0.05, type='scatter3d', mode='markers',
                       marker=list(reversescale=T, color=~ref.reads,
                       colorscale='Blues', showscale=T, opacity=0.50,
                       line=list(color='#EFF3FF', width=1),
                       colorbar=list(len=0.3, y=0.4, title='RefReads')),
                       name='Cells with only REFreads')
f <- f %>% add_trace(data=subset(y, (vaf==0 & ref.reads==0)|(is.na(vaf)==1)|(vaf>0)),
                       x=~x, y=~y, z=~z,
                       size=0.05, type='scatter3d', mode='markers',
                       marker=list(color='#EBEBEB', line=list(width=0),opacity=0.50), 
                       name='All other cells')
if (!args$`disable-slingshot`) {
    f <- f %>% add_trace(data=curves, x=~PC_1, y=~PC_2, z=~PC_3,
                         mode='lines', type='scatter3d', split=~Lineage)
}
if (args$`enable-title`) {
  f <- f %>% layout(title=paste0(snvs$CHROM[i], ':', snvs$POS[i], ' ',
                                 snvs$REF[i], '>', snvs$ALT[i], 
                                 ' (', sample.name, ')'),
                    scene=list(xaxis=list(title='PC_1'),
                               yaxis=list(title='PC_2'),
                               zaxis=list(title='PC_3')))
}
f
file.name <- paste0(sample.name, '-REFreads-', snvs$CHROM[i], '-',
                    snvs$POS[i], '-', snvs$REF[i], '-', snvs$ALT[i], '.html')
saveWidget(f, paste0(dir.name,sample.name, '_individual_sceSNVs/REFreads/', file.name),
           selfcontained=F, libdir='lib')
}
})}
