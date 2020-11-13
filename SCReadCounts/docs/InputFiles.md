# SCReadCounts Input Files

## SNVs

Single Nucleotide Genomic positions of interest. Tabular
formats and their required extensions include whitespace separated
text files (`.txt`), tab-separated values files
(`.tsv`), comma-separated values files (`.csv`), Excel (`.xlsx`),
and Excel 2003 (`.xls`).

Text files must have four white-space separated columns
representing the chromosome (CHROM), locus (POS), wild-type allele
nucleotide (REF), and SNV nucleotide (ALT). Other tabular formats must
provide CHROM, POS, REF, ALT headings. 

## Read Alignment Files

Read alignment files in indexed BAM format. Filename extension `.bam`
expected with `.bam.bai` index files in the same folder. readCounts will
execute fastest if all BAM files are sorted and indexed in a
consistent manner.

**All read alignemnts and SNV loci must indicate
genomic position with respect to the same specific release of a common
reference genome.**

