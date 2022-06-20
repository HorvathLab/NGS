# SCExecute Cell Barcodes

SCExecute splits the aligned reads by cell barcodes extracted from each alignment record. Cell barcodes can
be extracted from read headers by regular expression, splitting lines
according to some separator character, or directly from BAM alignment
headers. Aligned reads without a group identifier can be assigned a
specific identifier or skipped.

The available cell-barcode strategies are defined in the file
`group.ini` in the SCExecute distribution. New or modified
cell-barcode strategies can be created, using the same format, using 
a similar `group.ini` file in the current working directory. Named cell-barcode extraction rules in the current working directory override those with the
same name in the SCExecute distribution.

## Cell Barcode Extraction Rules

### UMI-tools
> Cell barcodes from read name added by UMI-tools.

### STARsolo
> Cell barcodes from the CB tag of aligned read - reads without a CB tag or with CB tag not in the accept list (default: file "barcodes.tsv" in the current directory) dropped

### CellRanger
> Cell barcodes from the CB tag of aligned read - reads without a CB tag or with CB tag not in the accept list (default: file "barcodes.tsv" in the current directory) dropped
                        
## Cell Barcode Extraction Operations

### ReadNameWord
> Parameters: field_index field_sep=_ missing=None

> Split the read name into words according to `field_sep` (default: "_"), and retain word with index `field_index` (required). Field index starts at 0. If the read name does not have enough words, use the read group identifier specified by `missing` (default: not specified). If `missing` is not specified, drop the read.

### ReadNameRegex
> Parameters: regex regexgrp=1 missing=None

> Apply (Python) regular expression `regex` (required) to the read name and extract matching group number `regexgrp` (default: 1). If the regular expression doesn't match, use the read group identifier specified by `missing` (default: not specified). If `missing` is not specified, drop the read.

### ReadTagValue
> Parameters: tag missing=None

> Use the value in the BAM tag `tag` (required). If `tag` is missing, use the read group specified by identifier `missing` (default: not specified). If `missing` is not specified, drop the read. 

### RGTag
> Parameters: missing=None

> Use the value in the BAM read-group tag "RG". If "RG" is missing, use the read group specified by identifier `missing` (default: not specified). If `missing` is not specified, drop the read. 

## Examples

### UMI-tools

```
[UMI-tools_CB]
Name: UMI-tools
Description: Cell barcodes from read name added by umi_tools.
Type: CellBarcode
ReadNameWord: field_index=1 field_sep=_
```

### STARsolo
```
[STARsolo_CB]
Name: STARsolo
Description: Cell barcodes from the CB tag of aligned read - reads without a CB tag or with CB tag not in the accept list (default: file "barcodes.tsv" in the current directory) dropped.                                                                
Type: CellBarcode                                                                                                            
ReadTagValue: tag='CB' acceptlist='barcodes.tsv'
```

### UMI-tools_Regex
```
[UMI-tools_Regex_CB]
Name: UMI-tools-RE
Description: Pick out cell barcodes from read name/identifier added by umi_tools using a regular expression.
Type: CellBarcode
ReadNameRegex: regex='_([ACGT]{16})_' regexgrp=1 missing="XXXXXXXXX"
```

## See Also

[SCExecute Home](..), [Usage](Usage.md), [Input Files](InputFiles.md), [Command Substitution](CommandSubst.md), [Examples](Examples.md)
