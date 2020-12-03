# Introduction
scReQTL is a strategy formulated when a method published in our lab&mdash;ReQTL\(RNA-eQTL\)&mdash;is applied on the data obtained through single-cell RNA sequencing technologies.

More info on ReQTL([paper](https://doi.org/10.1093/bioinformatics/btz750); [method](https://github.com/HorvathLab/ReQTL)).

# Workflow
![Pipeline](https://github.com/HorvathLab/NGS/blob/master/scReQTL/docs/pipeline.png?raw=true)
The scReQTL workflow includes three major components: scRNA-seq data processing, VAF<sub>RNA</sub> assessment, and SNV-GE correlation by cell type. Processing includes barcode and UMI processing, alignment, GE estimation and cell type classification, and can employ a variety of publicly available tools. 

The scReQTL-specific steps include the following steps:

1. Producing VAF<sub>RNA</sub> values through [SCReadCounts](https://github.com/HorvathLab/NGS/tree/master/SCReadCounts).
2. Gene-cell matrices pre-processing, clustering and assigning cell-types using [Seurat_and_SingleR](https://github.com/hliu5259/scReQTL).
3. scReQTL:
    1. Harmonize matrices [\(harmonize_matrices.R\)](https://github.com/HorvathLab/ReQTL/#harmonize_matricesr)
    2. Run matrix scReQTL [\(run_matrix_ReQTL.R\)](https://github.com/HorvathLab/ReQTL/#run_matrix_reqtlr)
    3. Annotate cis-trans scReQTLs [\(annotate_cis_trans.R\)](https://github.com/HorvathLab/ReQTL/#annotate_cis_transr)
    4. Plot the scReQTLs [\(Plot_scReQTLs_inbulk.R\)](README.md#plotting-screqtls)

# Plotting scReQTLs
To plot the scReQTLs, run the [Plot_scReQTLs_inbulk.R](https://github.com/HorvathLab/NGS/tree/master/scReQTL/docs/Plot_scReQTLs_inbulk.R) as `Rscript Plot_scReQTLs_inbulk.R -c <correlation_file.txt> -v <vaf_matrix_file.txt> -g <gene_expression_matrix_file.txt> -o <output_prefix> -t <top_#_correlations> -f <FDR_threshold>`

## Inputs
### Required Arguments
#### -c <correlation\_file.txt>:
The correlation file produced after running the [run_matrix_ReQTL.R](https://github.com/HorvathLab/ReQTL/#run_matrix_reqtlr) script.
#### -v <vaf\_matrix\_file.txt>:
The VAF matrix file produced by [SCReadCounts](https://github.com/HorvathLab/NGS/tree/master/SCReadCounts).
#### -g <gene\_expression\_matrix\_file.txt>:
The gene-cell matrix produced by featureCounts \(or a similar program/software\).

### Optional Arguments
#### -f <FDR\_threshold>:
The FDR threshold below which correlations are considered significant \(and filter out the rest\). By default, 0.05.
#### -t <top\_\#\_correlations>:
The number of significant correlations to retain \(once filtered using the [-f <FDR\_threshold>](README.md#-f-fdr_correlations)\)
#### -o <output\_prefix>:
The output prefix for the `.tiff` file produced by the script.

## Output
This script will produce 1 `.tiff` file.

Depending on the value supplied \(or not\) to [-o <output\_prefix>](README.md#-o-output_prefix) the name of the output file varies. By default, it will produce:
1. File called `Top_n_scReQTLs.tiff`, where `n` is the number of correlations to be plot \(see [-t <top\_\#\_correlations>](README.md#-t-top__correlations)\), when the output prefix is not supplied.
2. File called `<output_prefix>_top_n_scReQTLs.tiff`, where `n` is the number of correlations to be plot, when the output prefix is supplied.


# Authors and Acknowledgements

Hongyu Liu, Prashant N M, Liam Spurr, Pavlos Bousounis, Nawaf Alomran, Helen Ibeawuchi, Justin Sein, Dacian Reece-Stremtan, Piotr Słowiński, Krasimira Tsaneva-Atanasova, Muzi Li, Qianqian Zhang, Gabriel Asher, Keith A. Crandall, and Anelia Horvath

We would like to thank the Matrix eQTL team (Shabalin, et al. 2012) for their sample code and R package upon which run_matrix_ReQTL.R is based.
