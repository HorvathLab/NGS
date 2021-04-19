# ReadCounts Read Grouping

The counts for aligned reads are tabulated by BAM file and, if desired,
by a group identifier extracted from each alignment record. Usecases
include cell-barcodes added for single-cell sequencing. Read groups can
be extracted from read headers by regular expression, splitting lines
according to some separator character, or directly from BAM alignment
headers. Aligned reads without a group identifier can be assigned a
specific identifier or omitted from the output.

The available read-grouping strategies are defined in the file
`group.ini` in the ReadCounts distribution. New or modified
read-grouping strategies can be created, using the same format, in
a `group.ini` file in the current working directory. Named grouping
strategies in the current working directory override those with the
same name in the ReadCounts distribution.

## Read-Grouping Strategies

### UMI-tools
> Pick out cell barcodes from read name/identifier added by umi_tools.

### STARsolo
> Cell barcodes added by STARsolo as CB tag in aligned read, reads without a CB tag dropped.

## Read-Grouping Operations

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
[UMI-tools]
Description: Pick out cell barcodes from read name/identifier added by umi_tools.
ReadNameWord: field_index=1 field_sep=_
```

### STARsolo
```
[STARsolo]
Description: Cell barcodes added by STARsolo in CB tag in aligned read, reads without a CB tag dropped.
ReadTagValue: tag='CB'
```

### UMI-tools_Regex
```
[UMI-tools_Regex]
Description: Pick out cell barcodes from read name/identifier added by umi_tools using a regular expression.
ReadNameRegex: regex='_([ACGT]{16})_' regexgrp=1
```

### UB-Tag
```
[UB-Tag]
Description: UB tag from aligned read, reads without a UB tag get value "XXXXXXXX"
ReadTagValue: tag='UB' missing='XXXXXXXX'
```

