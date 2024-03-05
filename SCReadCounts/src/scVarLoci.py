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
  acceptlist=            # File of cell-barcodes to consider
  cellbarcodes=STARsolo  # Name of cell-barcode extraction strategy
  maxedits=5             # Max. edits for acceptable read alignment
  minbasequal=25         # Min. base quality score
  mincells=3             # Min. number of cells with variant
  mindist=5              # Min. number of bases between variants
  minmappingqual=60      # Min. mapping quality for reads
  minvarumipercell=3     # Min. number of UMIs with variant per cell
  minvarreadspercell=0   # Min. number of reads with variant per cell
  region=                # Restriction analysis to chrom. region(s)
  outfile=               # Place output in file, otherwise standard out
  outputcells=False      # Output the number of cells with the variant
  umibarcodes=STARsolo   # Name of UMI extraction strategy

Options:
  -v 	Verbose
  -h	Help
""".strip())
    sys.exit(0)

verbose = opt['v']

bamfilename = sys.argv[1]
kwargs = dict(region = None, minbasequal = 25, outfile = "", 
              minmappingqual = 60, mindist = 5, maxedits = 5,
              mincells = 3, minvarumipercell = 3, minvarreadspercell = 0,
              cellbarcodes = 'STARsolo', acceptlist="", outputcells=False)

for kv in sys.argv[2:]:
    badsplit=False
    try:
        k,v = map(str.strip,kv.split('='))
    except ValueError:
        badsplit = True
    if badsplit or not k or k not in kwargs:
        raise RuntimeError("Bad parameter setting: "+kv)
    try:
        v = eval(v)
    except:
        pass
    kwargs[k] = v

if not kwargs.get('umibarcodes',None):
    kwargs['umibarcodes'] = kwargs['cellbarcodes']
if kwargs.get('acceptlist') != None:
    if kwargs.get('acceptlist') in ("","None","-"):
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
outputcells = kwargs['outputcells']
del kwargs['outputcells'] 

scan = SCVariantLociByFetch(bamfilename,**kwargs)
if not outfile:
    outfile = sys.stdout
else:
    outfile = open(outfile,'wt')
if outputcells:
    print("CHROM","POS","REF","ALT","SNVCELLS",sep='\t',file=outfile)
else:
    print("CHROM","POS","REF","ALT",sep='\t',file=outfile)
for chrom,pos,refnuc,altnuc,freq in scan.loci():
    # print(chrom,pos,refnuc,altnuc,freq,sep='\t')
    if outputcells:
        ncells = 0
        for cb,umicnt in freq[altnuc].items():
            umis = len(umicnt.keys())
            if len(umicnt.keys()) < kwargs.get('minvarumipercell'):
                continue
            if sum(umicnt.values()) < kwargs.get('minvarreadspercell'): 
                continue
            ncells += 1
        print(chrom,pos,refnuc,altnuc,ncells,sep='\t',file=outfile)
    else:
        print(chrom,pos,refnuc,altnuc,sep='\t',file=outfile)

if outfile != sys.stdout:
    outfile.close()
