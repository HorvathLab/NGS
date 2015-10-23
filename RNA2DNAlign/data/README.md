# RNA2DNAlign Data

Supporing datafiles required for filtering or annotation of SNP loci

## RefSeq (Human) exon coordinates from UCSC

0. The human RefSeq exon coordinates (hg19) from UCSC are provided in
the RNA2DNAlign/data directory and can be ussed as provided. Steps 1 and
2 can be used to recreate this file or for another organism or assembly.

1. Download BED format RefSeq human coordinates (hg19) from UCSC:

   wget -O table.txt 'http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=449541223_Pr2eLcTGHmShfVff6FgGWadUVVmS&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_refGene&hgta_ctDesc=table+browser+query+on+refGene&hgta_ctVis=pack&hgta_ctUrl=&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbQual=cds&fbDownBases=200&hgta_doGetBed=get+BED' 

2. Extract loci, sort, and remove duplicates. 

   cat table.txt | sed -e 's/^chr//' -e 's/^X/23/' -e 's/^Y/24/' | sort -k1n,1 -k2n,2 -k3n,3 | sed -e 's/^23/X/' -e 's/^24/Y/' | awk '$1 !~ /_/ {print $1"\t"$2"\t"$3}' | uniq > sorted_exon_coordintes.txt

   The format for sorted_exon_coordintes.txt is tab-separated chromosome number (1,2,3,...,X,Y), start position, end position, in increasing order.

## COSMIC Mutants

1. Register with COSMIC here:

   https://cancer.sanger.ac.uk/cosmic/register

2. Download the COSMIC mutants:

   sftp://sftp-cancer.sanger.ac.uk//files/grch38/cosmic/v74/CosmicMutantExport.tsv.gz

3. This file is used in its published format.

## DARNED loci

1. Download the DARNED loci from here (NCBI37/hg19):

   http://darned.ucc.ie/static/downloads/hg19.txt

2. This file is used in its published format.
