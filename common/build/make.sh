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
BASE=`readlink -f "$BASE"`
cd $BASE
if [ ! -d "common" -o ! -d "common/build" ]; then
    echo "Please change directory to the base of the HorvathLabTools distribution" 1>&2
    exit 1;
fi
PACKAGE=`readlink -f "$PACKAGE"`
PACKAGE=`basename "$PACKAGE"`
if [ ! -d "./$PACKAGE" ]; then
    echo "Valid packages: SNPlice, RNA2DNAlign, ReadCounts" 1>&2
    exit 1;
fi
VER=`apython3 $PACKAGE/src/version.py VERSION | tr -d -c '0-9.'`
OS=`uname`
AR=`uname -m`
XX="$OS-$AR"
YY="Python-3.7"
PROGS=`apython3 $PACKAGE/src/version.py PROGRAMS`

# Source (Python-3.7) distribution
rm -rf build/$PACKAGE-${VER}.${YY} dist/$PACKAGE-${VER}.${YY}.tgz
mkdir -p build/$PACKAGE-${VER}.${YY}
INCLUDES=`apython3 $PACKAGE/src/version.py INCLUDES`
for d in src docs data scripts; do
 mkdir -p build/$PACKAGE-${VER}.${YY}/$d
 for p in $INCLUDES $PACKAGE; do
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
    ./venv/bin/pyinstaller --hidden-import pkg_resources.py2_warn --hidden-import pysam.libctabixproxies --hidden-import json --distpath $PACKAGE/bin build/$PACKAGE-${VER}.${YY}/src/$p
    # unzip -lv $PACKAGE/bin/$base/base_library.zip > $PACKAGE/bin/${base}_library.toc
    rsync -av $PACKAGE/bin/$base/ $PACKAGE/bin/_bin
    rm -rf $PACKAGE/bin/$base
    cp common/scripts/wrapper.sh $PACKAGE/bin/$base
    # /tools/anaconda3/bin/pyinstaller --hidden-import pysam.libctabixproxies --hidden-import json --distpath $PACKAGE/bin build/$PACKAGE-${VER}.${YY}/src/$p
    rm -rf build/$base
    rm -f ${base}.spec
  fi
done

mkdir -p build/$PACKAGE-${VER}.${XX}

for d in bin docs data scripts; do
  mkdir -p build/$PACKAGE-${VER}.${XX}/$d
  for p in $INCLUDES $PACKAGE; do
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
find build/$PACKAGE-${VER}.${XX} -type d -empty -exec rm -rf {} \;
find build/$PACKAGE-${VER}.${YY} -type d -name __pycache__ -exec rm -rf {} \;

mkdir -p dist
tar -czf dist/$PACKAGE-${VER}.${XX}.tgz -C build $PACKAGE-${VER}.${XX}
tar -czf dist/$PACKAGE-${VER}.${YY}.tgz -C build $PACKAGE-${VER}.${YY}
( cd dist; md5sum $PACKAGE-${VER}.*.tgz > $PACKAGE-${VER}.md5 )
if [ "$COMMIT" -eq 1 ]; then
  git commit -a -m "Release $PACKAGE-${VER} commit"; git push
fi
if [ "$TAG" -eq 1 ]; then
  git tag -f $PACKAGE-${VER}
fi
if [ "$COMMIT" -eq 1 ]; then
  git push -f --tags
fi
