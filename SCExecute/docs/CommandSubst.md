# SCExecute Command and Filename Template Substutition

SCExecute can replace a variety of placeholder terms in the specified command and/or filename templates to ensure cell-barcode specific BAM files can have a cell-barcode specific command execution, output capture, and working directory, as needed. The placeholders are specified inside curly braces `{` and `}` and are uppercase (case-sensitive), as shown below.

## Generated cell-barcode specific BAM file placeholders

### {}

> Full path to generated cell-barcode speicific BAM file.

### {CBPATH}

> Full path to generated cell-barcode specific BAM file.

### {CBFILE}

> Filename part of generated cell-barcode specific BAM file.

### {CBBASE}

> Filename part of generated cell-barcode specific BAM file, without the .bam extension.

## Original input BAM file placeholders

### {BAMPATH}

> Full path of original input BAM file. 

### {BAMFILE}

> Filename part of original input BAM file. 

### {BAMBASE}

>  Filename part of original input BAM file, without the .bam extension.

## Miscellaneous

### {BARCODE}

> The cell-barcode for a cell-barcode specific BAM file. 

### {BFINDEX}

> Index of original input BAM file.

### {CBINDEX}

> Index of cell-barcode from its original input BAM file.

### {WORKER}

> Processor index for a command execution.

## Example

### Filename template
```
{BAMBASE}.{BARCODE}.bam
```

### Logfile template
```
{BAMBASE}.{BARCODE}.log
```

### Command template
```
gatk HaplotypeCaller -R human.fasta -I {} -O out/{BFINDEX}.{CBINDEX}.vcf
```
