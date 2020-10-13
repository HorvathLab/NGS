#!/bin/sh
DIR=`dirname $0`
# DIR=`readlink -f "$DIR"`
PROG=`basename $0`
if [ ! -d ${DIR}/../src ]; then
  LD_LIBRARY_PATH="${DIR}/_bin:${LD_LIBRARY_PATH}:"
  export LD_LIBRARY_PATH
  if [ -d $DIR/_bin/Contents ]; then
    exec "$DIR/_bin/Contents/MacOS/$PROG" "$@"
  else
    exec "$DIR/_bin/$PROG" "$@"
  fi
else
  if ${PYTHON3:-python3} </dev/null >/dev/null 2>&1; then
    exec ${PYTHON3:-python3} "$DIR/../src/${PROG}.py" "$@"
  else
    echo "Please add python3 to your path or set environment variable PYTHON3 to its location." | fmt 
  fi
fi
