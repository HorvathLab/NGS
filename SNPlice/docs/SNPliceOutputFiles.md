# SNPlice Output Files

## Output Files

Tabular formats including txt, tsv, csv, xlsx, xls.

Text files (txt) are tab-separated fields without headers. Other tabular formats provide headers. Fields are output in deterministic column order as shown on SNPlice Output Fields page. 

## Output Fields

Note that fields and their values passed through from the input SNP file (VCF or tabular formats) are not documented here, refer to the source of the input SNP file for this.

CHROM  
> Chromosome identifier of SNP

POS  
> Chromosome position of SNP

REF  
> Reference allele nucleotide

ALT  
> Variant allele nucleotide

NumofJuncs  
> Number of exon-intron junctions within the specified distance (-d/--distance/Distance)

Distance  
> Nucleotide distance from SNP locus to exon-intron junction.

Junctions
> Chromosome position of the exon-intron-exon junctions represented by each intron's boundaries.

SNPJuncIntronCount  
> Number of spliced spanning reads with ALT nucleotide.

SNPJuncNoIntronCount	Number of unspliced spanning reads with ALT nucleotide.
NoSNPJuncIntronCount	Number of spliced spanning reads with REF nucleotide.
NoSNPJuncNoIntronCount	Number of unspliced spanning reads with REF nucleotide.
Probability	Probability score of observed read counts.
LOD 	Log-odds of unspliced spanning reads vs spliced spanning reads with ALT nucleotide with respect to unspliced spanning reads vs spliced spanning reads.
P-Value	Fisher exact test of 2x2 contingency table for spliced vs unspliced spanning reads containing ALT vs REF nucleotides.
Bonferroni	Bonferroni multiple test correction of Fisher exact test p-value.
FDR	Benjaminni-Hochberg false-discovery-rate multiple test correction of Fisher exact test p-value.
SNPMateCount (M,F)	Number of spanning reads with the ALT nucleotide whose mate contains the exon-intron junction.
NoSNPMateCount (M,F)	Number of spanning reads with the REF nucleotide whose mate contains the exon-intron junction.
SNPCount (F)	Number of spanning reads with ALT nucleotide.
NoSNPCount (F)	Number of spanning reads with REF nucleotide.
MatesCount (M,F)	Number of spanning reads containing the SNP site whose mate contains the exon-intron junction.
NotMatesCount (M,F)	Number of spanning reads containing the SNP site and the exon-intron junction.
IntronCount (F)	Number of spliced spanning reads.
NoIntronCount (F)	Number of unspliced spanning reads.
SpanningReads (F)	Number of reads (or reads + mates) containing the SNP site and exon-intron junction.
RemovedDuplicateReads (F)	Number of duplicate reads removed.
SNPLociReads (F)	Number of reads containing SNP site.

Key:

    (M) Fields output if -M/--matepairs/Mates option is true. 

    (F) Fields output if -F/--full/All Fields option is true. 

## See Also

[SNPlice Home](..), [SNPlice](SNPliceUsage.md), [SNPlice-Combine](SNPliceCombineUsage.md), [Input Files](SNPliceInputFiles.md)


