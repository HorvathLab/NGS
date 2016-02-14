#!/bin/sh
if [ "$1" = "" ]; then
  echo "Usage: $0 <UCSC-assembly>"
  exit 1
fi
# Why is there no web-api for the hgTables stuff? 
URL="http://hgdownload.cse.ucsc.edu/goldenPath/$1/database/refGene.txt.gz"
wget -q -O - "$URL" | gunzip -c | \
  fgrep -w cmpl | awk '$3 ~ /^chr..?$/ {print $3,$9,$10,$11,$7,$8}' | fgrep -v chrUn | fgrep -v '_' | \
  sed -e 's/^chr//' -e 's/^X/100/' -e 's/^Y/101/' | \
  awk '{split($3,s,","); split($4,e,","); for (i=1;i<=$2;i++) {if (e[i] >= $5 && s[i] <= $6) { si = s[i]; if ($5 > si) { si = $5 }; ei = e[i]; if ($6 < ei) { ei = $6 }; print $1"\t"si"\t"ei}}}' | \
  sort -k1n,1 -k2n,2 -k3n,3 | \
  sed -e 's/^100/X/' -e 's/^101/Y/' | \
  uniq 
