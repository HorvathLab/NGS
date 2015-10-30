# SNPlice Input Files

## Read Alignments

Read alignments in (indexed) BAM/SAM format. Filename extension `.bam` expected, with `.bam.bai` index files. 

## Junctions

Exon-intron-exon junctions in BED format. Filename extension `.bed` expected.

Junctions are represented as two exon "blocks" separated by an intron, as output by TopHat's junction output option.

## SNPs

Single-nucleotide-polymophisms (SNPs) in tabular and VCF format. Tabular formats/extensions include txt, tsv, csv, xlsx, xls.

Text files (txt) must have four white-space separated columns representing the chromosome (CHROM), locus (POS), wild-type allele nucleotide (REF), and SNP nucleotide (ALT). Other tabular formats must provide CHROM, POS, REF, ALT headings. Extra values in tabular or VCF format files are mapped to the output.

## See Also

[SNPlice Home](..), [SNPlice](SNPliceUsage.md), [SNPlice-Combine](SNPliceCombineUsage.md)

