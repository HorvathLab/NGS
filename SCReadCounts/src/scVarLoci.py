#!/bin/env python3
import sys
import os
import os.path
import traceback
from os.path import join, dirname, realpath
try:
    sys.path.append(join(dirname(realpath(__file__)),
                         '..', '..', 'common', 'src'))
except NameError:
    pass

from collections import defaultdict

from util import ReadGroupFactory

groupFactory = ReadGroupFactory()
groupOptions = [t[1] for t in groupFactory.list(type="CellBarcode")] + [""]
groupMap = {}
for s,n,d in sorted(groupFactory.list(type="CellBarcode")):
    groupMap[n] = s

umiOptions = [""] + [t[1] for t in groupFactory.list(type="UMI")]
umiMap = {}
for s,n,d in sorted(groupFactory.list(type="UMI")):
    umiMap[n] = s

from variantloci import SCVariantLociByFetch

bamfilename = sys.argv[1]
region = None
mincells = 2
minumipercell = 3
minqual = 30
mindist = 5
cellbarcodes,umis = sys.argv[2].split(':')
assert cellbarcodes in groupMap
assert umis in umiMap
cellbarcodes = groupFactory.get(groupMap[cellbarcodes])
umis = groupFactory.get(umiMap[umis])

if len(sys.argv) > 3:
    mincells,minumipercell = map(int,sys.argv[3].split(":"))
if len(sys.argv) > 4:
    region = sys.argv[4]
    if region.strip() in ("", "-", "None"):
        region = None

scan = SCVariantLociByFetch(bamfilename,region=region,
                            minqual=minqual,mindist=mindist,
                            cellbarcodes=cellbarcodes,umis=umis,
                            mincells=mincells,minumipercell=minumipercell)
for chrom,pos,refnuc,altnuc,freq in scan.loci():
    # anf = []
    # for cb,umis in freq.get(altnuc,{}).items():
    #     anf.append(len(umis))
    # ans = sum(1 for _ in filter(lambda x: x>=minumipercell,anf))
    # print(chrom,pos,refnuc,altnuc,ans,sep='\t')
    print(chrom,pos,refnuc,altnuc,sep='\t')
