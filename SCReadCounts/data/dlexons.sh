#!/bin/sh
if [ "$1" == "-h" ]; then
  echo "Usage: $0 [ <EMBL-Release> ]" 1>&2
  exit 1
fi
if [ "$1" == "" ]; then
  REL1="current_gtf"
else
  REL1="release-$1/gtf"
fi
FN=`wget -q -O - "https://ftp.ensembl.org/pub/$REL1/homo_sapiens/" | sed -n 's/^.*href="\(.*.[0-9][0-9]*.gtf.gz\)".*$/\1/p'`
REL2=`echo "$FN" | tr '.'  ' ' | awk '{print $(NF-2)}'` 
ASS=`echo "$FN" | tr '.'  ' ' | awk '{print $(NF-3)}'` 
URL="https://ftp.ensembl.org/pub/$REL1/homo_sapiens/Homo_sapiens.$ASS.$REL2.gtf.gz"
echo "Downloading Homo_sapiens.$ASS.$REL2.gtf.gz" 1>&2
wget -q -O - "$URL" | gunzip -c | \
  grep -v '^#' | \
  awk '$3 == "exon" {print $1,$4,$5}' | \
  awk '$1 ~ /^([12]?[0-9]|X|Y|MT)$/' | \
  sed -e 's/^X /100 /' -e 's/^Y /101 /' -e 's/^MT /102 /' | \
  sort -k1n,1 -k2n,2 -k3n,3 | \
  sed -e 's/^100 /X /' -e 's/^101 /Y /' -e 's/^102 /MT /' | \
  tr ' ' '\t' | \
  uniq 
