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

from variantloci import VariantLociByPileup, VariantLociByFetch

bamfilename = sys.argv[1]
region = None
mincnt = 1
if len(sys.argv) > 2:
    region = sys.argv[2]
    if region.strip() in ("", "-", "None"):
        region = None
if len(sys.argv) > 3:
    mincnt = int(sys.argv[3])

# scan = VariantLociByPileup(bamfilename,region=region,mincnt=mincnt)
scan = VariantLociByFetch(bamfilename,region=region,mincnt=mincnt)
for chrom,pos,refnuc,altnuc,freq in scan.loci():
    # print(chrom,pos,refnuc,altnuc,freq.get(refnuc,0),freq.get(altnuc,0),sep='\t')
    print(chrom,pos,refnuc,altnuc,sep='\t')
