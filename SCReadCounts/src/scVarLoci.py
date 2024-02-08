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

opt = defaultdict(lambda: False)
while len(sys.argv) > 1 and sys.argv[1] in ("-v","-h"):
   opt[sys.argv.pop(1)[1]] = True

if len(sys.argv) <= 1 or opt['h']:
    print("""
Usage:
  scVarLoci [ options ] <bamfile> [ parameters ]

Parameters (with defaults):
  acceptlist=
  cellbarcodes=STARsolo
  minbasequal=25
  mincells=3
  mindist=5
  minmappingqual=60
  minvarumipercell=3
  region=
  outfile=
  umibarcodes=STARsolo

Options:
  -v 	Verbose
  -h	Help
""".strip())
    sys.exit(0)

verbose = opt['v']

bamfilename = sys.argv[1]
kwargs = dict(region = None, minbasequal = 25,
              minmappingqual = 60, mindist = 5,
              mincells = 3, minvarumipercell = 3,
              cellbarcodes = 'STARsolo', acceptlist=None)

for kv in sys.argv[2:]:
    k,v = kv.split('=')
    try:
        v = eval(v)
    except:
        pass
    kwargs[k] = v

if kwargs.get('cellbarcodes',None) != None:
    if not kwargs.get('umibarcodes',None):
        kwargs['umibarcodes'] = kwargs['cellbarcodes']
    if kwargs.get('acceptlist') != None:
        if kwargs.get('acceptlist') in (None,"","None","-"):
            readgroupparam = "*:acceptlist=None"
        else:
            readgroupparam = "*:acceptlist='%s'"%(kwargs.get('acceptlist',))
    else:
        readgroupparam = ""
    kwargs['cellbarcodes'] = groupFactory.get(groupMap[kwargs['cellbarcodes']],readgroupparam)
    kwargs['umibarcodes'] = groupFactory.get(umiMap[kwargs['umibarcodes']])

if verbose:
    print("Parameters:",file=sys.stderr)
    for k,v in sorted(kwargs.items()):
        if k in ('cellbarcodes','umibarcodes'):
            print("  %s=%s"%(k,v._name),file=sys.stderr)
        else:
            if v == None:
                print("  %s="%(k,),file=sys.stderr)
            else:
                print("  %s=%s"%(k,v),file=sys.stderr)
    print(file=sys.stderr)

outfile = kwargs.get('outfile')
if 'outfile' in kwargs:
    del kwargs['outfile']

scan = SCVariantLociByFetch(bamfilename,**kwargs)
if not outfile:
    outfile = sys.stdout
else:
    outfile = open(outfile,'wt')
for chrom,pos,refnuc,altnuc,freq in scan.loci():
    # print(chrom,pos,refnuc,altnuc,freq,sep='\t')
    altcellcnt = sum(1 for _ in filter(lambda x: x>=kwargs.get('minvarumipercell'),map(lambda k: len(freq.get(altnuc,{}).get(k,{})),filter(lambda k: k != "-",freq.get(altnuc,{}).keys()))))
    print(chrom,pos,refnuc,altnuc,altcellcnt,sep='\t',file=outfile)

if outfile != sys.stdout:
    outfile.close()
