# SCExecute Examples

Split BAM file based on STARsolo cell barcodes into indexed cell-specific BAM files, and for each, write the last ten aligned reads into a cell-specific SAM file. Run up to three commands at a time. 
```
% bin/scExecute -r data/singlecell2_117.bam -G STARsolo -b data/barcodes.tsv -t 3 -i -C "samtools view {} | tail > test_{BAMBASE}_{BARCODE}_{CBINDEX}.sam"
```

## See Also

[SCExecute Home](..), [Usage](Usage.md), [Input Files](InputFiles.md), [Cell Barcodes](Barcodes.md), [Command/Template Substitution](CommandSubst.md), [Examples](Examples.md)

