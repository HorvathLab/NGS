#plot snp sith gene
load_package <- function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, dep = TRUE)
        if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
    }
}

load_package("data.table");load_package("dplyr");load_package("ggplot2");load_package("ggpubr");load_package("tidyr");load_package("stringr");load_package("optparse");

# manages command line arguments
option_list <- list(
  make_option(c("-c", "--corr"), action = "store", type ="character", 
              default = NULL, help= "File containing correlations to plot"),
  make_option(c("-g", "--ge"), action="store", type = "character", default=NULL,
              help="Gene expression (Gene-cell) matrix file"),
  make_option(c("-v", "--vaf"), action = "store", type ="character", 
              default = NULL, help= "VAF matrix file"),
  make_option(c("-t", "--top"), action = "store", type ="numeric", 
              default = 50, help= "Provide the number of correlations to plot. DEFAULT: 50"),
  make_option(c("-o", "--outpref"), action = "store", type ="character", 
              default = "Top50", help= "Output prefix for the TIFF file with correlation plots. DEFAULT: Top<number of correlations to plot>.tiff")
)

opt <- parse_args(OptionParser(option_list=option_list, description = "-c, -g, -v, and -o options are necessary!!!!", usage = "usage: Rscript Plot_scReQTLS_inbulk.R -c <correlation_file.txt> -v <vaf_matrix_file.txt> -g <gene_expression_matrix_file.txt> -o <output_prefix> -t <top_#_correlations>"))

#read_in correlations to plot
in_file = fread(opt$corr)
in_file = in_file %>% filter(FDR < 0.05)
in_file = in_file %>% head(opt$top)

to_plot = in_file %>% dplyr::select(SNP, gene)
names(to_plot) = c("SNP", 'gene')


#read_in VAF/VAF_star matrix
VAF_star = fread(opt$vaf) %>% filter(SNV %in% to_plot$SNP)
VAF_star = VAF_star %>% gather(sample, vaf, 2:ncol(VAF_star)) %>% drop_na()
colnames(VAF_star)[1] <- 'SNP'
#read_in GE matrix
GE_star = fread(opt$ge)
GE_star = GE_star %>% filter(gene_id %in% to_plot$gene) %>% gather(sample, fpkm, 2:ncol(GE_star)) %>% drop_na()

#plot_ReQTL_star
df_star = left_join(to_plot, VAF_star, by = 'SNP')
df_star = left_join(df_star, GE_star, by = c('gene' = 'gene_id',"sample"))
df_star = df_star %>% drop_na() %>% unite(id, SNP, gene)
names(df_star) = c("id", "sample", "SNP", "gene")
p = ggscatter(df_star, x = "SNP", y = "gene",
              fill = "lightsteelblue", color = "lightsteelblue", shape = 21, size = 1.5, # Points color, shape and size
              add = "reg.line",  # Add regression line
              add.params = list(color = "dodgerblue4", fill = alpha("dodgerblue4"), 0.5), # Customize reg. line
              cor.coef = T, # Add correlation coefficient. see ?stat_cor
              cor.coef.size = 5,
              cor.coeff.args = list(method = "spearman", label.sep = "\n"))

tiff(filename = ifelse(opt$top != 50 & opt$outpref == "Top50", paste0(gsub("50", opt$top, opt$outpref), ".tiff"), paste0(opt$outpref, "_top_", opt$top, ".tiff")), width = 25, height = 40, units = "in", res = 120)
facet(p, facet.by = "id", ncol = 5, scales = "free_y")
dev.off()



