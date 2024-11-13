
suppressWarnings(library('optparse'))

script.desc <-
'Generates statistics to determine if SNV distribution has role in clustering.
User must supply the RDS of a Seurat object (-r) which should contain a single
Seurat experiment. A file also must be submitted with (-t) which is a file that
contains a path to all the scReadCounts files containing the SNV information.'

parser <- OptionParser(description=script.desc)
parser <- add_option(parser, c('-r', '--rds-file'),
                     type='character',
                     help='RDS file containing Seurat object.')
parser <- add_option(parser, c('-t', '--snv-file'),
                     type='character',
                     help='scReadCounts output file containing SNVs.')
parser <- add_option(parser, c('-w', '--th-snv-cells'),
                     type='numeric', default=10,
                     help=paste0('Threshold for maximum percentage of cells ',
                     ' that contain an SNV for bad reads. ', 
                     '(Default: 10)'))
parser <- add_option(parser, c('-z', '--th-snv-reads'),
                     type='numeric', default=1,
                     help=paste0('Threshold for number of minimum reads to ',
                     'qualify as SNV. (Default: 1)'))
args <- parse_args(parser)

error.msg <- NULL
if (is.null(args$`rds-file`))
  error.msg <- paste(error.msg, '- Seurat RDS object (-r) is required.',
                     sep='\n')
if (is.null(args$`snv-file`))
  error.msg <- paste(error.msg, '- SNV data (-t) is required.', sep='\n')
if (args$`th-snv-reads`<=0)
  error.msg <- paste(error.msg, '- th-snv-reads needs to be greater than 0.',
                     sep='\n')
if (!is.null(error.msg)) {
  print_help(parser)
  stop(error.msg)
}


suppressWarnings(library('dplyr'))
library('Seurat')

# filtering based on arguments passed
srat <- readRDS(args$`rds-file`)
snv <- read.table(args$`snv-file`, sep= '\t', header=T)
sample.name <- gsub('(.*/)*([A-Za-z0-9]+)_.*.rds', '\\2', args$`rds-file`)

# handling missing VAFs
snv$VAF[snv$VAF == '-'] <- NA
snv$VAF <- suppressWarnings(as.numeric(snv$VAF))
snv <- snv[!is.na(snv$VAF), ]

# filtering based on `X.BadRead` and `SNVCount`
snv.temp <- snv
snv.temp['BadReadFlag'] = 0
snv.temp[snv.temp[['X.BadRead']]>0, 'BadReadFlag'] = 1
snv.read.filt <- aggregate(BadReadFlag~CHROM+POS+REF+ALT, data=snv.temp,
                           function (x) 100*sum(x)/length(x))
snv.read.filt <- snv.read.filt[snv.read.filt$BadReadFlag<=args$`th-snv-cells`, 
                               c('CHROM', 'POS', 'REF', 'ALT')]
snv <- merge(snv, snv.read.filt, by=c('CHROM', 'POS', 'REF', 'ALT'))

# filter on `SNVCount`
snv$VAF[snv$SNVCount<args$`th-snv-reads`] <- 0
if (nrow(snv) == 0)
    stop('There are no rows left after filtering.')

# cluster assignment (mergnig SNV data with seurat_clusters)
df.cid <- as.data.frame(srat[['seurat_clusters']])
df.cid <- data.frame(ReadGroup = rownames(df.cid),
                     ClusterID = df.cid[, 1], row.names = NULL)
n.cluster <- length(unique(df.cid$ClusterID))
df.snv <- merge(snv, df.cid, by = 'ReadGroup')

# pre-split the data by unique SNV combinations
snv_groups <- split(df.snv, paste(df.snv$CHROM, df.snv$POS, df.snv$REF, df.snv$ALT, sep = "_"))

# Kruskal-Wallis test for each SNV across clusters
kw_test_statistics <- c()  
kw_p_values <- c()        

for (snv_key in names(snv_groups)) {
    snv_data <- snv_groups[[snv_key]]

    # perform the test iff there are multiple clusters with data for the SNV
    if (length(unique(snv_data$ClusterID)) > 1) {
        kw_test_result <- kruskal.test(VAF ~ ClusterID, data = snv_data)
        kw_test_statistics <- c(kw_test_statistics, kw_test_result$statistic)  # H-statistic
        kw_p_values <- c(kw_p_values, kw_test_result$p.value)
    } else {
        kw_test_statistics <- c(kw_test_statistics, NA)  # NA if test is not applicable
        kw_p_values <- c(kw_p_values, NA)
    }
}

# adjustment for multiple comparisons using Bonferroni correction
kw_p_adj <- p.adjust(kw_p_values, method = 'bonferroni')

df.final <- data.frame(
    SNV = names(snv_groups),
    H_statistic = kw_test_statistics, 
    p_value = kw_p_values,
    p_adj = kw_p_adj
)

df.final <- na.omit(df.final)                  # filtering NA cases
df.final <- df.final[order(df.final$p_adj), ]  # sorting by p_adj

write.table(df.final, paste0(sample.name, '_snv_candidate_', 
            sprintf('%02d', args$`th-snv-cells`), 'n', sprintf('%02d', args$`th-snv-reads`), 
            'n.txt'), sep = '\t', row.names = FALSE)

significant_snvs <- df.final[df.final$p_adj < 0.1, ]
write.table(significant_snvs, paste0(sample.name, '_snv_candidate_', 
            sprintf('%02d', args$`th-snv-cells`), 'n', sprintf('%02d', args$`th-snv-reads`), 
            'n_q01.txt'), sep = '\t', row.names = FALSE)
