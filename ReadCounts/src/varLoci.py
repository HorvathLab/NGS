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
minreads = 1
if len(sys.argv) > 2:
    minreads = int(sys.argv[2])
if len(sys.argv) > 3:
    region = sys.argv[3]
    if region.strip() in ("", "-", "None"):
        region = None

# scan = VariantLociByPileup(bamfilename,region=region,minreads=minreads)
scan = VariantLociByFetch(bamfilename,region=region,minreads=minreads)
for chrom,pos,refnuc,altnuc,freq in scan.loci():
    # print(chrom,pos,refnuc,altnuc,freq.get(refnuc,0),freq.get(altnuc,0),sep='\t')
    print(chrom,pos,refnuc,altnuc,sep='\t')
