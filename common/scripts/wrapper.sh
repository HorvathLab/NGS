#!/bin/sh
DIR=`dirname $0`
DIR=`readlink -f "$DIR/../lib"`
PROG=`basename $0`
LD_LIBRARY_PATH=":${LD_LIBRARY_PATH}:$DIR"
export LD_LIBRARY_PATH
exec "$DIR/$PROG" "$@"

