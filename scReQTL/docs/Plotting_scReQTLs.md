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