# March 10, 2023

suppressWarnings(library('optparse'))

script.desc <-
'Calculates differential expression using DESeq2.
User must supply the RDS of a Seurat object (-r) which must contain cell identities annotated by SingleR.
User must also supply a scReadCounts output file (-t) which contains the SNV
position and number of reads.
User can supply cell type (-c) or a short list of SNVs (-l) to analyze.'

parser <- OptionParser(description=script.desc)
parser <- add_option(parser, c('-r', '--rds-file'),
                     type='character',
                     help='RDS file containing Seurat object.')
parser <- add_option(parser, c('-t', '--snv-file'),
                     type='character',
                     help='File containg SNVs.')
parser <- add_option(parser, c('-c', '--cell-types'),
                     type='character',
                     help='List of cell types to analyze.')
parser <- add_option(parser, c('-l', '--short-list'),
                     type='character',
                     help='File containg short list of SNVs.')
parser <- add_option(parser, c('-x', '--use-sctype'),
                     type='logical', default=F, action='store_true',
                     help='Use scType instead of SingleR.')
args <- parse_args(parser)

error.msg <- NULL

# Check if the required argument (-r) is passed
if (is.null(args$`rds-file`))
  error.msg <- paste(error.msg, '- Seurat RDS object (-r) is required.',
                     sep='\n')
if (is.null(args$`snv-file`))
  error.msg <- paste(error.msg, '- SNV data (-t) is required.',
                     sep='\n')

# Check if there are any errors
if (!is.null(error.msg)) {
  print_help(parser)
  stop(error.msg)
}


######
# Application Logic
######

library('Seurat')
library('EnhancedVolcano')
library('DESeq2')
srt <- readRDS(args$`rds-file`)
snv.data <- read.table(args$`snv-file`, sep='\t', header=T)

if (args$`use-sctype`) { # Hack to make it compatiable with scType
  srt[['SingleR.cluster.labels']] <- srt[['customclassif']]
}

calculateDEx <- function(srt, snv.sub, cell.type) {
  # Check if SNV list only has one SNV position
  stopifnot(nrow(unique(snv.sub[, 1:4]))==1)
  
  # Find cell types
  bc <- srt[['SingleR.cluster.labels']]

  # Differentiate cell types
  bc[bc$SingleR.cluster.labels!=cell.type, 'SNV'] <- 'Irrelevant'
  bc[bc$SingleR.cluster.labels==cell.type, 'SNV'] <- 'Undetected'

  # Find barcodes that contain 0 and non-0 VAF
  bc.ref <- snv.sub$ReadGroup[snv.sub$VAF==0]
  bc.ref <- bc.ref[!is.na(bc.ref)]
  bc.alt <- snv.sub$ReadGroup[snv.sub$VAF>0]
  bc.alt <- bc.alt[!is.na(bc.alt)]

  # Set reference or alternate
  bc[(rownames(bc) %in% bc.ref) & (bc$SNV=='Undetected'), 'SNV'] <- 'Ref'
  bc[(rownames(bc) %in% bc.alt) & (bc$SNV=='Undetected'), 'SNV'] <- 'Alt'

  srt.snv <- srt
  stopifnot(rownames(bc)==colnames(srt.snv))

  # Assign SNV information to Seurat meta-data
  srt.snv[['SNV']] <- bc$SNV
  
  if (sum(srt.snv[['SingleR.cluster.labels']]==cell.type) == 0)
    stop('Cell type not found in dataset.')
  if (sum(srt.snv[['SNV']]=='Alt') < 3 || sum(srt.snv[['SNV']]=='Ref') < 3)
    stop('Not enough Ref/Alt reads for analysis.')
  suppressMessages(
    res.dex <- FindMarkers(srt.snv, ident.1='Alt', ident.2='Ref',
                           group.by='SNV', test.use='DESeq2', assay='RNA')
  )
  res.dex <- res.dex[order(res.dex$p_val_adj, res.dex$p_val), ]

  # Print out status msg before returning value
  message(paste0('\nSNV: ', snv.sub[1, 1], ':', snv.sub[1, 2], ' ',
                            snv.sub[1, 3], '>', snv.sub[1, 4], '\n',
                 'Cell type: ', cell.type), '\n',
                 'With Ref: ', sum(bc$SNV=='Ref'), '; ',
                 'With SNV: ', sum(bc$SNV=='Alt'))
  res.dex
}


snv.list <- unique(snv.data[, 1:4])
if (!is.null(args$`cell-types`)) {
  cell.types <- unlist(strsplit(gsub('[\'\"]', '', args$`cell-types`), ','))
} else {
  cell.types <- unique(srt[['SingleR.cluster.labels']][, 1])
}
message(paste0('Analyzing: ', paste(unlist(cell.types), collapse=', ')))

# Reduce the list
if (!is.null(args$`short-list`)) {
  snv.short <- read.table(args$`short-list`, sep='\t', header=T,
                          colClasses=c('character', 'character', 'character',
                          'character')) 
  snv.list <- merge(snv.list, snv.short, by=c('CHROM', 'POS', 'REF', 'ALT'))
  snv.list <- unique(snv.list)
}


suppressWarnings(dir.create('output'))
for (i in 1:nrow(snv.list)) {
  for (cell.type in cell.types) {
    # Put together name
    data.name <- substr(gsub(' ', '', cell.type), 1, 5)
    data.name <- paste0(data.name, '_',
                        snv.list$CHROM[i], '_',
                        snv.list$POS[i], '_', 
                        snv.list$REF[i], '_', 
                        snv.list$ALT[i])

    snv.sub <- subset(snv.data, CHROM==snv.list$CHROM[i] &
                                POS==snv.list$POS[i] &
                                REF==snv.list$REF[i] &
                                ALT==snv.list$ALT[i])
    tryCatch({
      res <- calculateDEx(srt, snv.sub, cell.type)
      write.table(res, paste0('output/DE_', data.name, '.txt'),
                  sep='\t', col.names=NA)
      suppressMessages({
        EnhancedVolcano(res, lab=rownames(res), x='avg_log2FC', y='p_val',
                        title=paste0(snv.sub[1, 1], ':', snv.sub[1, 2], ' ',
                              snv.sub[1, 3], '>', snv.sub[1, 4]),
                        subtitle=NULL, pCutoff=0.05)
        ggsave(paste0('output/DE_', data.name, '_volc.png'))
      })
      },
      error=function (e) {
        warning(paste0('Issue with ', data.name, '.\n'), e$message)}
    )
  }
}
