# RNA2DNA Input Files

## SNVs

Single-nucleotide-variants (SNVs) in tabular or VCF format. Tabular
formats and their required extensions include whitespace separated
text files (`.txt`), tab-separated values files
(`.tsv`), comma-separated values files (`.csv`), Excel (`.xlsx`),
and Excel 2003 (`.xls`).

Text files must have four white-space separated columns
representing the chromosome (CHROM), locus (POS), wild-type allele
nucleotide (REF), and SNV nucleotide (ALT). Other tabular formats must
provide CHROM, POS, REF, ALT headings. Extra values in tabular or VCF
format files are mapped to the output.

All SNV loci, read alignments, and annotation files *must* indicate
genomic position with respect to the same specific release of a common
reference genome.

## Read Alignment Files

Read alignment files in indexed BAM format. Filename extension `.bam`
expected with `.bam.bai` index files in the same folder. All read
alignemnts, SNV loci, and annotation files *must* indicate genomic
position with respect to the same specific release of a common
reference genome. RNA2DNA will execute faster if all BAM file
alignments are sorted in a consistent manner.

## See Also

[RNA2DNA Home](..), [RNA2DNA Usage](Usage.md), [Output Files](OutputFiles.md), [Annotation Files](AnnotationFiles.md)

