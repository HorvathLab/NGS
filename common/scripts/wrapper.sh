#!/bin/sh
DIR=`dirname $0`
DIR=`readlink -f "$DIR"`
PROG=`basename $0`
LD_LIBRARY_PATH="${DIR}/_bin:${LD_LIBRARY_PATH}:"
export LD_LIBRARY_PATH
exec "$DIR/_bin/$PROG" "$@"

