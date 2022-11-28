# November 28, 2022

suppressWarnings(library('optparse'))

script.desc <-
'Generates plots containing SNV information.
User must supply the RDS of a Seurat object (-r) which can either contain a
single dataset or multiple integrated datasets. A file also must be submitted
with (-t) which is a file that contains a path to all the scReadCounts files
containing the SNV information. This file should be listed in the order that it
appears in the integrated Seurat object. Lastly, if desired, a third argument
can be passed which lists the individual SNVs (-l) for which visualizations can
be generated.'

parser <- OptionParser(description=script.desc)
parser <- add_option(parser, c('-r', '--rds-file'),
                     type='character',
                     help='RDS file containing Seurat object.')
parser <- add_option(parser, c('-t', '--snv-file'),
                     type='character',
                     help='Tab-delimited file containg SNVs.')
parser <- add_option(parser, c('-l', '--target-snvs'),
                     type='character',
                     help='List of files containing individual SNVs to plot.')
parser <- add_option(parser, c('-m', '--min-threshold'), 
                     type='character',
                     help='Minimum threshold number of reads.')
parser <- add_option(parser, c('-x', '--plot-width'),
                     type='integer',
                     help='Width of split plots')
parser <- add_option(parser, c('-y', '--plot-height'),
                     type='integer',
                     help='Height of split plots')
args <- parse_args(parser)

# Default parameters
min.threshold <- 3 # Minimum number of reads for VAF calculation

plot.split.width <- 14 # Default width for split plot
plot.split.height <- 7 # Default height for split plot
plot.split.units <- 'in'

error.msg <- NULL
# Check if only one of two arguments (-s, -r) is passed
if (is.null(args$`rds-file`))
  error.msg <- paste(error.msg, '- Seurat RDS object (-r) is required.',
                     sep='\n')

# Check if SNV argument is passed
if (is.null(args$`snv-file`))
  error.msg <- paste(error.msg, '- SNV data (-t) is required.', sep='\n')

# Check if there are any errors before proceeding
if (!is.null(error.msg)) {
  print_help(parser)
  stop(error.msg)
}

# If the minimum threshold is set, set it to user-supplied value
if (!is.null(args$`min-threshold`))
  min.threshold <- args$`min-threshold`

# If the plot width is set, set it to user-supplied value
if (!is.null(args$`plot-width`))
  plot.split.width <- args$`plot-width`
  
if (!is.null(args$`plot-height`))
  plot.split.height <- args$`plot-height`


######
# APPLICATION LOGIC
######
library('Seurat')
library('patchwork')
library('RColorBrewer')
library('ggplot2')

# Calculate number of different SNVs per cell from file
CalculateSNVCount <- function(raw.snv) {
  snv.data <- subset(raw.snv, raw.snv$SNVCount > 0)
  snv.data <- aggregate(snv.data$ReadGroup, list(snv.data$ReadGroup), length)
  colnames(snv.data) <- c('ReadGroup', 'Count')
  snv.data <- snv.data[order(-snv.data$Count), ]
  rownames(snv.data) <- 1:nrow(snv.data)
  snv.data
}

# Calculate total SNV reads for each cell from file
CalculateTotalSNVCount <- function(raw.snv) {
  snv.data <- subset(raw.snv, raw.snv$SNVCount > 0)
  snv.data <- aggregate(snv.data$SNVCount, list(snv.data$ReadGroup), sum)
  colnames(snv.data) <- c('ReadGroup', 'Count')
  snv.data <- snv.data[order(-snv.data$Count), ]
  rownames(snv.data) <- 1:nrow(snv.data)
  snv.data
}

# Process read groups
LoadAggregateSNVReadGroups <- function(snv.file) {
  snvs <- read.table(snv.file)
  i <- 1
  full.snv <- data.frame()
  for (f in snvs$V1) {
    raw.snv <- read.delim(f)
    new.index <- rep(paste0('_', i), length(raw.snv$ReadGroup))
    raw.snv$ReadGroup <- paste0(raw.snv$ReadGroup, new.index)
    full.snv <- rbind(full.snv, raw.snv)
    i <- i+1
  }
  full.snv
}

# Generate specific SNV information based on cell
IndividualSNVCells <- function(data, full.snv, snv.list) {
  snv.data <- subset(full.snv, full.snv$SNVCount > 0)
  snv.short <- merge(snv.list, snv.data, c('CHROM', 'POS'))
  gene.pos <- unique(snv.short[, c('CHROM', 'POS')])
  for (i in 1:nrow(gene.pos)) {
    sub.snv.data <- subset(snv.data, snv.data$CHROM == gene.pos$CHROM[i] &
                           snv.data$POS == gene.pos$POS[i])
    data[['indiv.snv']] <- 'without'
    data[['indiv.snv']][colnames(data) %in% sub.snv.data$ReadGroup, ] <- 'with'

    Idents(data) <- data[['indiv.snv']]

    # Make the dimensional reduction plots and save to appropriate directory
    DimPlot(data, order=c('without', 'with'), cols=c('black', 'gray95')) +
      plot_annotation(title=paste0(gene.pos$CHROM[i], ':', gene.pos$POS[i]))
    ggsave(paste0('plots/indiv/umap.indiv.whole.', gene.pos$CHROM[i], '.',
                  gene.pos$POS[i], '.png'))

    DimPlot(data, split.by='expt', order=c('without', 'with'),
            cols=c('black', 'gray95')) +
      plot_annotation(title=paste0(gene.pos$CHROM[i], ':', gene.pos$POS[i]))
    ggsave(paste0('plots/indiv/umap.indiv.split.', gene.pos$CHROM[i], '.',
                  gene.pos$POS[i], '.png'), width=plot.split.width,
                  height=plot.split.height, units=plot.split.units)
  }
}

# Generate specific VAF for each cell for a given SNV
IndividualVAFCells <- function(data, full.snv, snv.list) {
  snv.short <- merge(snv.list, full.snv, c('CHROM', 'POS'))
  gene.pos <- unique(snv.short[, c('CHROM', 'POS')])
  for (i in 1:nrow(gene.pos)) {
    sub.snv.data <- subset(snv.short, snv.short$CHROM == gene.pos$CHROM[i] &
                           snv.short$POS == gene.pos$POS[i])
    data[['indiv.vaf']] <-
      sub.snv.data$VAF[match(colnames(data),
                             sub.snv.data$ReadGroup, nomatch=NA)]
    data[['indiv.vaf.scaled']] <- ceiling(data[['indiv.vaf']]*4)

    # Apply humamn-readable labels
    data[['indiv.vaf.scaled']][data[['indiv.vaf.scaled']]==0] <- '0'
    data[['indiv.vaf.scaled']][data[['indiv.vaf.scaled']]==1] <- '(0, 25]'
    data[['indiv.vaf.scaled']][data[['indiv.vaf.scaled']]==2] <- '(25, 50]'
    data[['indiv.vaf.scaled']][data[['indiv.vaf.scaled']]==3] <- '(50, 75]'
    data[['indiv.vaf.scaled']][data[['indiv.vaf.scaled']]==4] <- '(75, 100]'
    Idents(data) <- data[['indiv.vaf.scaled']]
    id.order = c('(75, 100]', '(50, 75]', '(25, 50]', '(0, 25]', '0')

    # Make the dimensional reduction plots and save to appropriate directory
    DimPlot(data, order=id.order,
            cols=brewer.pal(n=5, name='YlOrRd'), na.value='gray95') +
      plot_annotation(title=paste0(gene.pos$CHROM[i], ':', gene.pos$POS[i]))
    ggsave(paste0('plots/indiv/umap.indiv.whole.vaf.', gene.pos$CHROM[i], '.',
                  gene.pos$POS[i], '.png'))

    DimPlot(data, split.by='expt', order=id.order,
            cols=brewer.pal(n=5, name='YlOrRd'), na.value='gray95') + 
      plot_annotation(title=paste0(gene.pos$CHROM[i], ':', gene.pos$POS[i]))
    ggsave(paste0('plots/indiv/umap.indiv.split.vaf.', gene.pos$CHROM[i], '.',
                  gene.pos$POS[i], '.png'), width=plot.split.width,
                  height=plot.split.height, units=plot.split.units)
  }
}

###
# DATA PREPARATION
# ---
# Loads data from files
# Calculates and scales SNV data for visualization
###

# Load the files into variables.
data <- readRDS(args$`rds-file`)
full.snv <- LoadAggregateSNVReadGroups(args$`snv-file`)
snvs <- read.table(args$`snv-file`)

# Apply thresholding
full.snv.thresh <- subset(full.snv, full.snv$GoodReads >= min.threshold)
write.table(full.snv, 'output.vaf.txt', sep='\t')

# Load SNV's into a data frame
data[['expt']] <- NA
data[['n.snv']] <- NA
data[['reads.snv']] <- NA

# Calculate SNV information for the dataset
full.n.snv <- CalculateSNVCount(full.snv)
full.reads.snv <- CalculateTotalSNVCount(full.snv)

# Label the integrated Seurat object with appropriate sample labels
i <- 1
for (f in snvs$V1) {
  data[['expt']][grep(paste0('_', i, '$'), colnames(data)), ] <- snvs$V2[i]
  i <- i+1
}

# Add SNV information to main Seurat object
# Calculate number of different SNVs in each cell
data[['n.snv']][, 1] <-
  full.n.snv[match(colnames(data), full.n.snv$ReadGroup), 2]
data[['n.snv']][is.na(data[['n.snv']])] <- 0
data[['n.snv.scaled']] <-
  ceiling(data[['n.snv']]/max(data[['n.snv']])*8)

# Calculate number of SNV reads in each cell
data[['reads.snv']][, 1] <-
  full.reads.snv[match(colnames(data), full.reads.snv$ReadGroup), 2]
data[['reads.snv']][is.na(data[['reads.snv']])] <- 0
data[['reads.snv.scaled']] <-
  ceiling(data[['reads.snv']]/max(data[['reads.snv']])*8)

# Create data frame for histogram
snv.df <- data.frame(expt=data[['expt']], n.snv=data[['n.snv']],
                     reads.snv=data[['reads.snv']])
snv.df.ord <- data.frame(reads.snv=snv.df$reads.snv[order(snv.df$reads.snv)],
                         n.snv=snv.df$n.snv[order(snv.df$n.snv)])
snv.df.ord$x <- 1:nrow(snv.df.ord)
print('test2')

###
# VISUALIZATION SECTION
###
# Make directory
suppressWarnings(dir.create('plots'))
suppressWarnings(dir.create('plots/indiv'))

# Theme changes
theme_set(theme_classic())
theme_update(plot.title=element_text(hjust = 0.5))

# Generate histograms for number of SNVs
message('\nGenerating histograms for number of SNVs...')
ggplot(snv.df)+geom_histogram(aes(x=n.snv), binwidth=1.0)+
  xlab('Number of SNVs in Cell')+ylab('Count')
ggsave('plots/hist.whole.n.png')

ggplot(snv.df)+geom_histogram(aes(x=n.snv), binwidth=1.0)+
  xlab('Number of SNVs in Cell')+ylab('Count')+facet_wrap(~ expt)
ggsave('plots/hist.split.n.png')

# Generate histograms for SNV reads
message('\nGenerating histograms for number of SNV reads...')
ggplot(snv.df)+geom_histogram(aes(x=reads.snv), binwidth=1.0)+
  xlab('Number of SNV Reads in Cell')+ylab('Count')
ggsave('plots/hist.whole.snv.png')

ggplot(snv.df)+geom_histogram(aes(x=reads.snv), binwidth=1.0)+
  xlab('Number of SNV Reads in Cell')+ylab('Count')+
  facet_wrap(~ expt)
ggsave('plots/hist.split.reads.png')

# Generate cumulative SNV plot
message('\nGenerating cumulative SNV plot...')
ggplot(snv.df.ord)+aes(x=x, y=reads.snv)+geom_point(color='blue', size=0.1)
ggsave('plots/cum.whole.reads.png')

ggplot(snv.df.ord)+aes(x=x, y=n.snv)+geom_point(color='blue', size=0.1)
ggsave('plots/cum.whole.n.png')

# Generate UMAP with gradient scale for number of SNVs
message('\nGenerating UMAP representation for number of SNVs...')
DimPlot(data, group.by='n.snv.scaled', order=8:0,
        cols=brewer.pal(n=9, name='YlOrRd'))
ggsave('plots/umap.whole.n.gr.png')

DimPlot(data, group.by='n.snv.scaled', split.by='expt',
        order=8:0, cols=brewer.pal(n=9, name='YlOrRd'))
ggsave('plots/umap.split.n.gr.png', width=plot.split.width,
       height=plot.split.height, units=plot.split.units)

# Generate UMAP with gradient scale for SNV reads
message('\nGenerating UMAP representation for number of SNV reads...')
DimPlot(data, group.by='reads.snv.scaled', order=8:0,
        cols=brewer.pal(n=9, name='YlOrRd'))
ggsave('plots/umap.whole.reads.gr.png')

DimPlot(data, group.by='reads.snv.scaled', split.by='expt',
        order=8:0, cols=brewer.pal(n=9, name='YlOrRd'))
ggsave('plots/umap.split.reads.gr.png', width=plot.split.width,
       height=plot.split.height, units=plot.split.units)

# Process individual SNV plots
# This function requires Seurat object to be labeled with 'expt' meta-data
if (!is.null(args$`target-snvs`)) {
message('\nGenerating individual UMAP representations for individual SNVs...')
  snv.list <- read.delim(args$`target-snvs`)
  IndividualSNVCells(data, full.snv, snv.list)
  IndividualVAFCells(data, full.snv.thresh, snv.list)
}
