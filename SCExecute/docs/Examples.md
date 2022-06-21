# SCExecute Examples

Split BAM file based on STARsolo barcodes into indexed barcode-stratified BAM files, export the last ten aligned reads in text format into specific filename. Analyze three different cell-barcode specific analyses at a time.
```
% cd data
% ../bin/scExecute -r singlecell2_117.bam -G STARsolo -i -C "samtools view {} | tail > test_{BAMBASE}_{BARCODE}_{CBINDEX}.sam" -t 3
```

## See Also

[SCExecute Home](..), [Usage](Usage.md), [Input Files](InputFiles.md), [Cell Barcodes](Barcodes.md), [Command/Template Substitution](CommandSubst.md), [Examples](Examples.md)

