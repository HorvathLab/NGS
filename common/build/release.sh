#!/bin/sh
if [ "$1" = "" ]; then
   echo "Usage: release.sh <release-tag>"
   exit 1
fi
echo md5sum $1.*tgz > $1.md5
echo gh release create $1 -F $1.txt $1.*tgz $1.md5
