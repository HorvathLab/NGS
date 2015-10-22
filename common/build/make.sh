#!/bin/sh
PACKAGE="$1"
BASE=`dirname "$0"`
BASE="$BASE/../.."
BASE=`readlink -f "$BASE"`
cd $BASE
if [ ! -d "common" -o ! -d "common/build" ]; then
    echo "Please change directory to the base of the HorvathLabTools distribution"
fi
if [ ! -d "$PACKAGE" ]; then
    echo "Valid packages: SNPlice, RNA2DNAlign" 1>&2
    exit 1;
fi
VER=`python27 $PACKAGE/src/version.py | tr -d -c '0-9.'`
OS=`uname`
AR=`uname -m`
XX="$OS-$AR"
rm -rf build/$PACKAGE-${VER}.${XX} dist/$PACKAGE-${VER}.${XX}.tgz
PROGS=`fgrep -l 'from version import' $PACKAGE/src/*.py`
rm -rf $PACKAGE/bin
export TCL_LIBRARY=/tools/EPD/lib/tcl8.5
for p in $PROGS; do
  /tools/EPD/bin/cxfreeze --include-path=common/src --include-modules=hashlib,ctypes,platform,pysam.TabProxies --target-dir=$PACKAGE/bin $p
done
mkdir -p build/$PACKAGE-${VER}.${XX}
for d in bin data scripts; do
  mkdir build/$PACKAGE-${VER}.${XX}/$d
  if [ -d common/$d ]; then
    rsync -a common/$d build/$PACKAGE-${VER}.${XX}
  fi
  if [ -d $PACKAGE/$d ]; then
    rsync -a $PACKAGE/$d build/$PACKAGE-${VER}.${XX}
  fi
done
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
# rm -rf build/$PACKAGE-${VER}.${XX}
