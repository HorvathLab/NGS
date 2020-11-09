#!/bin/sh
set -x
TAG=0
COMMIT=0
while true; do
    if [ "$1" = "tag" ]; then
	TAG=1
	COMMIT=1
	shift;
    elif [ "$1" = "commit" ]; then
	COMMIT=1
	shift;
    else
        break
    fi
done
PACKAGE="$1"
if [ "$PACKAGE" = "" ]; then
    echo "Usage: make.sh PACKAGE" 1>&2
    exit 1;
fi
BASE=`dirname "$0"`
BASE="$BASE/../.."
# BASE=`readlink -f "$BASE"`
cd $BASE
if [ ! -d "common" -o ! -d "common/build" ]; then
    echo "Please change directory to the base of the HorvathLabTools distribution" 1>&2
    exit 1;
fi
# PACKAGE=`readlink -f "$PACKAGE"`
PACKAGE=`basename "$PACKAGE"`
if [ ! -d "./$PACKAGE" ]; then
    echo "Valid packages: SNPlice, RNA2DNAlign, ReadCounts" 1>&2
    exit 1;
fi
OS=`uname`
AR=`uname -m`
XX="$OS-$AR"
YY="Python-3.7"

if [ "$OS" = "Darwin" ]; then
  PYTHON3=python3
  PYINST="pyinstaller -w"
  MD5SUM="md5 -r"
  XX="macOS-$AR"
else
  PYTHON3=./venv/bin/python
  PYINST=./venv/bin/pyinstaller
  MD5SUM=md5sum
fi
VER=`$PYTHON3 $PACKAGE/src/release.py VERSION | tr -d -c '0-9.'`
PROGS=`$PYTHON3 $PACKAGE/src/release.py PROGRAMS`

# Source (Python-3.7) distribution
rm -rf build/$PACKAGE-${VER}.${YY} dist/$PACKAGE-${VER}.${YY}.tgz
mkdir -p build/$PACKAGE-${VER}.${YY}
INCLUDES=`$PYTHON3 $PACKAGE/src/release.py INCLUDES`
for d in src scripts; do
 mkdir -p build/$PACKAGE-${VER}.${YY}/$d
 for p in $INCLUDES $PACKAGE; do
  if [ -d $p/$d ]; then
    rsync --copy-links -a $p/$d build/$PACKAGE-${VER}.${YY}
  fi
 done
done
for d in docs data; do
 mkdir -p build/$PACKAGE-${VER}.${YY}/$d
 for p in $PACKAGE; do
  if [ -d $p/$d ]; then
    rsync --copy-links -a $p/$d build/$PACKAGE-${VER}.${YY}
  fi
 done
done
mkdir build/$PACKAGE-${VER}.${YY}/bin
for p in $PROGS; do
  base=`basename $p .py`
  cp build/$PACKAGE-${VER}.${YY}/scripts/wrapper.sh build/$PACKAGE-${VER}.${YY}/bin/$base
done
rm build/$PACKAGE-${VER}.${YY}/scripts/wrapper.sh
find build/$PACKAGE-${VER}.${YY} -name ".svn" -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${YY} -name "*.pyc" -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${YY} -name "*~" -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${YY} -type d -empty -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${YY} -type d -name __pycache__ -exec rm -rf {} \;

# Binary distribution
rm -rf build/$PACKAGE-${VER}.${XX} dist/$PACKAGE-${VER}.${XX}.tgz
rm -rf $PACKAGE/bin
# export TCL_LIBRARY=/tools/EPD/lib/tcl8.5
mkdir -p $PACKAGE/bin
for p in $PROGS; do
  if [ -f build/$PACKAGE-${VER}.${YY}/src/$p ]; then
    base=`basename $p .py`
    rm -rf build/$base
    rm -f ${base}.spec
    rm -rf $PACKAGE/bin/$base $PACKAGE/bin/$base.app
    $PYINST --hidden-import pkg_resources.py2_warn --hidden-import pysam.libctabixproxies --hidden-import json --distpath $PACKAGE/bin build/$PACKAGE-${VER}.${YY}/src/$p
    if [ -d $PACKAGE/bin/$base.app ]; then
      rsync -av $PACKAGE/bin/$base.app/ $PACKAGE/bin/_bin
      mkdir -p $PACKAGE/bin/_bin/Contents/MacOS/tcl
      mkdir -p $PACKAGE/bin/_bin/Contents/MacOS/tk
    else
      rsync -av $PACKAGE/bin/$base/ $PACKAGE/bin/_bin
      cp ./venv/lib/libstdc++.so.6 $PACKAGE/bin/_bin
    fi
    rm -rf $PACKAGE/bin/$base $PACKAGE/bin/$base.app
    cp common/scripts/wrapper.sh $PACKAGE/bin/$base
    rm -rf build/$base
    rm -f ${base}.spec
  fi
done

mkdir -p build/$PACKAGE-${VER}.${XX}

for d in bin scripts; do
  mkdir -p build/$PACKAGE-${VER}.${XX}/$d
  for p in $INCLUDES $PACKAGE; do
  if [ -d $p/$d ]; then
    rsync --copy-links -a $p/$d build/$PACKAGE-${VER}.${XX}
  fi
  done
done
for d in docs data; do
  mkdir -p build/$PACKAGE-${VER}.${XX}/$d
  for p in $PACKAGE; do
  if [ -d $p/$d ]; then
    rsync --copy-links -a $p/$d build/$PACKAGE-${VER}.${XX}
  fi
  done
done
for p in $INCLUDES $PACKAGE; do
  rsync --copy-links -a $p/src/*.ini build/$PACKAGE-${VER}.${XX}/bin
done
# mv build/$PACKAGE-${VER}.${XX}/bin build/$PACKAGE-${VER}.${XX}/lib
# mkdir -p build/$PACKAGE-${VER}.${XX}/bin
# for s in $PROGS; do
#   base=`basename $s .py`
#   cp common/scripts/wrapper.sh build/$PACKAGE-${VER}.${XX}/bin/$base
# done
rm build/$PACKAGE-${VER}.${XX}/scripts/wrapper.sh
find build/$PACKAGE-${VER}.${XX} -name ".svn" -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${XX} -name "*.pyc" -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${XX} -name "*~" -exec rm -rf {} \;
# find build/$PACKAGE-${VER}.${XX} -type d -empty -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${YY} -type d -name __pycache__ -exec rm -rf {} \;

mkdir -p dist
tar -czf dist/$PACKAGE-${VER}.${XX}.tgz -C build $PACKAGE-${VER}.${XX}
tar -czf dist/$PACKAGE-${VER}.${YY}.tgz -C build $PACKAGE-${VER}.${YY}
( cd dist; $MD5SUM $PACKAGE-${VER}.*.tgz > $PACKAGE-${VER}.md5 )
if [ "$COMMIT" -eq 1 ]; then
  git commit -a -m "Release $PACKAGE-${VER} commit"; git push
fi
if [ "$TAG" -eq 1 ]; then
  git tag -f $PACKAGE-${VER}
fi
if [ "$COMMIT" -eq 1 ]; then
  git push -f --tags
fi
