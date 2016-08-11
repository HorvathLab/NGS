#!/bin/sh
SRC=`dirname $0`/../src
SRC=`readlink -f $SRC`
EXONS="UCSC_Human_hg19_RefSeq_CDS_exon_coordinates.txt"
DARNED="DARNED_hg19.txt"
COSMIC="CosmicMutantExport_hg19.tsv.gz"

function testcmd() {
  NUMBER=$1
  BASE=.testing-output
  shift
  rm -rf $BASE-$NUMBER
  mkdir -p $BASE-$NUMBER
  echo -n "Test $NUMBER: "
  $SRC/RNA2DNAlign.py "$@" -o $BASE-$NUMBER > $BASE-$NUMBER/logfile.txt 2>&1
  if [ $? -eq 0 ]; then
    echo "Success"
    cksum `find $BASE-$NUMBER -type f | sort | grep -v logfile` > $BASE-$NUMBER.cksum.tmp
    if [ -f $BASE-$NUMBER.cksum ]; then
	# sdiff $BASE-$NUMBER.cksum $BASE-$NUMBER.cksum.tmp
	if diff -q $BASE-$NUMBER.cksum $BASE-$NUMBER.cksum.tmp >/dev/null 2>&1; then
            # echo "No change in output"
            true
        else
            echo "Unexpected change in:"
	    diff $BASE-$NUMBER.cksum $BASE-$NUMBER.cksum.tmp | fgrep '>' | awk '{print $4}' 
        fi
	rm -f $BASE-$NUMBER.cksum.tmp
    else
	mv -f $BASE-$NUMBER.cksum.tmp $BASE-$NUMBER.cksum
    fi
  else
    echo "Failed"
    echo RNA2DNAlign.py "$@"
    cat $BASE-$NUMBER/logfile.txt
  fi
  rm -rf $BASE-$NUMBER
}

testcmd 1  -r 'example-*.bam' -s 'example-*.vcf' -m 3 -e $EXONS -d $DARNED -c $COSMIC
testcmd 2  -r 'example-*.bam' -s 'example-*.vcf' -m 3 -e $EXONS 
testcmd 3  -r 'example-*.bam' -s 'example-*.vcf' -m 3 -d $DARNED -c $COSMIC
testcmd 4  -r 'example-*.bam' -s 'example-*.vcf' -m 3
testcmd 5  -r 'example-*.bam' -s example-SNV.tsv -m 3 -e $EXONS -d $DARNED -c $COSMIC
testcmd 6  -r 'example-*.bam' -s example-SNV.tsv -m 3 -e $EXONS
testcmd 7  -r 'example-*.bam' -s example-SNV.tsv -m 3 -d $DARNED -c $COSMIC
testcmd 8  -r 'example-*.bam' -s example-SNV.tsv -m 3
testcmd 9  -r 'example-*.bam' -s example-SNV.tsv -m 100
testcmd 10 -r 'example-*.bam' -s example-SNV.tsv
testcmd 11 -r 'example-GDNA.bam example-SDNA.bam' -s example-SNV.tsv -m 3
testcmd 12 -r 'example-GDNA.bam example-NRNA.bam' -s example-SNV.tsv -m 3
testcmd 13 -r 'example-SDNA.bam example-NRNA.bam example-TRNA.bam' -s example-SNV.tsv -m 3
testcmd 14 -r 'example-GDNA.bam example-SDNA.bam example-NRNA.bam example-TRNA.bam' -s example-SNV.tsv -m 3
testcmd 15 -r 'example-*.bam' -s 'example-*.vcf' -t 4
testcmd 16 -r 'example-*.bam' -s 'example-*.vcf' -t 0
testcmd 17 -r 'example-*.bam' -s 'example-*.vcf' -m 3 -M 0.5
testcmd 18 -r 'example-*.bam' -s 'example-*.vcf' -m 3 -M 100
