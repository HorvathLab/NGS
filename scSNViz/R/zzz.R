# Declaring global variables so that R CMD check doesn’t flag non-standard evaluations or missing bindings.

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "FindVariableFeatures", "Matrix", "MeanVAF", "NormalizeData",
    "SNV.N", "SNVCount", "ScaleData", "TotalVAF", "after_stat",
    "aggregate", "cluster", "count", "gene_sets_prepare", "head",
    "median", "na.omit", "orig.ident", "read.table", "sampleid",
    "scale_fill_manual", "scores", "sctype_score", "theme_minimal",
    "vaf_label", "write.table"
  ))
}