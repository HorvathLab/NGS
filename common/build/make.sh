#!/bin/sh
# set -x
PACKAGE="$1"
if [ "$PACKAGE" = "" ]; then
    echo "Usage: make.sh PACKAGE" 1>&2
    exit 1;
fi
BASE=`dirname "$0"`
BASE="$BASE/../.."
BASE=`readlink -f "$BASE"`
cd $BASE
if [ ! -d "common" -o ! -d "common/build" ]; then
    echo "Please change directory to the base of the HorvathLabTools distribution" 1>&2
    exit 1;
fi
PACKAGE=`readlink -f "$PACKAGE"`
PACKAGE=`basename "$PACKAGE"`
if [ ! -d "./$PACKAGE" ]; then
    echo "Valid packages: SNPlice, RNA2DNAlign" 1>&2
    exit 1;
fi
VER=`python27 $PACKAGE/src/version.py VERSION | tr -d -c '0-9.'`
OS=`uname`
AR=`uname -m`
XX="$OS-$AR"
YY="Python-2.7"

# Source (Python-2.7) distribution
rm -rf build/$PACKAGE-${VER}.${YY} dist/$PACKAGE-${VER}.${YY}.tgz
mkdir -p build/$PACKAGE-${VER}.${YY}
for d in src data scripts; do
  mkdir -p build/$PACKAGE-${VER}.${YY}/$d
  if [ -d common/$d ]; then
    rsync -a common/$d build/$PACKAGE-${VER}.${YY}
  fi
  if [ -d $PACKAGE/$d ]; then
    rsync -a $PACKAGE/$d build/$PACKAGE-${VER}.${YY}
  fi
done
find build/$PACKAGE-${VER}.${YY} -name ".svn" -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${YY} -name "*.pyc" -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${YY} -name "*~" -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${YY} -type d -empty -exec rm -rf {} \;
rm build/$PACKAGE-${VER}.${YY}/scripts/wrapper.sh

# Binary distribution
rm -rf build/$PACKAGE-${VER}.${XX} dist/$PACKAGE-${VER}.${XX}.tgz
PROGS=`python27 $PACKAGE/src/version.py PROGRAMS`
rm -rf $PACKAGE/bin
export TCL_LIBRARY=/tools/EPD/lib/tcl8.5
for p in $PROGS; do
  /tools/EPD/bin/cxfreeze --include-path=common/src --include-modules=hashlib,ctypes,platform,pysam.TabProxies,numpy.core --target-dir=$PACKAGE/bin $PACKAGE/src/$p
done
mkdir -p build/$PACKAGE-${VER}.${XX}
for d in bin data scripts; do
  mkdir -p build/$PACKAGE-${VER}.${XX}/$d
  if [ -d common/$d ]; then
    rsync -a common/$d build/$PACKAGE-${VER}.${XX}
  fi
  if [ -d $PACKAGE/$d ]; then
    rsync -a $PACKAGE/$d build/$PACKAGE-${VER}.${XX}
  fi
done
rsync -a common/src/*.ini build/$PACKAGE-${VER}.${XX}/bin
rsync -a $PACKAGE/src/*.ini build/$PACKAGE-${VER}.${XX}/bin
mv build/$PACKAGE-${VER}.${XX}/bin build/$PACKAGE-${VER}.${XX}/lib
mkdir -p build/$PACKAGE-${VER}.${XX}/bin
for s in $PROGS; do
  base=`basename $s .py`
  cp common/scripts/wrapper.sh build/$PACKAGE-${VER}.${XX}/bin/$base
done
rm build/$PACKAGE-${VER}.${XX}/scripts/wrapper.sh
find build/$PACKAGE-${VER}.${XX} -name ".svn" -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${XX} -type d -empty -exec rm -rf {} \;

mkdir -p dist
tar -czf dist/$PACKAGE-${VER}.${XX}.tgz -C build $PACKAGE-${VER}.${XX}
tar -czf dist/$PACKAGE-${VER}.${YY}.tgz -C build $PACKAGE-${VER}.${YY}
( cd dist; md5sum $PACKAGE-${VER}.*.tgz > $PACKAGE-${VER}.md5 )
git tag -d $PACKAGE-${VER}
git tag $PACKAGE-${VER}
