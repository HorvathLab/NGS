#!/bin/env python3

import sys, os, copy, traceback, re, time
from os.path import join, dirname, realpath
try:
    sys.path.append(join(dirname(realpath(__file__)),
                         '..', '..', 'common', 'src'))
except NameError:
    pass

from collections import defaultdict
from operator import itemgetter

from util import ReadGroupFactory

from rcio import scrcbinwriter, scrctxtwriter

from selectedloci import SCSelectedLociByFetch

opt = defaultdict(lambda: False)
while len(sys.argv) > 1 and sys.argv[1] in ("-v","-h"):
   opt[sys.argv.pop(1)[1]] = True

import sys, time

if len(sys.argv) <= 1 or opt['h']:
    print("""
Usage:
  scCountLoci [ options ] <bamfile> <snvlocifile> [ parameters ]

Parameters (with defaults):
  acceptlist=            # File of cell-barcodes to consider
  cellbarcodes=STARsolo  # Name of cell-barcode extraction strategy
  maxedits=3             # Max. edits for acceptable read alignment
  minbasequal=25         # Min. base quality score
  mincells=3             # Min. number of cells for cell and variant constraints
  minmappingqual=60      # Min. mapping quality for reads
  minreadspercell=3      # Min. number of reads per cell
  minumipercell=3        # Min. number of UMIs per cell
  minvarreadspercell=3   # Min. number of variant reads per cell
  minvarumipercell=3     # Min. number of variant UMIs per cell
  outfile=               # Place output in compact binary file
  outputumis=True        # Output UMI counts, not reads
  umibarcodes=STARsolo   # Name of UMI extraction strategy

Options:
  -v    Verbose
  -h    Help
""".strip())
    sys.exit(0)

verbose = opt['v']

bf = sys.argv[1]
sf = sys.argv[2]
kwargs = {'minbasequal': 25, 'minmappingqual': 60, 'gapsize': 1000, 'maxedits': 3, 
          'cellbarcodes': 'STARsolo', 'mincells': 3, 'minvarumipercell': 3, 'minumipercell': 3,
          'outfile': "", 'acceptlist': "", 'minvarreadspercell': 3, 'minreadspercell': 3, 'outputumis': True}
for kv in sys.argv[3:]:
    k,v = kv.split('=')
    try:
       v = eval(v)
    except:
       pass
    kwargs[k] = v
snvs = [ l.split()[:4] for l in open(sf) ]

groupFactory = ReadGroupFactory()
groupOptions = [t[1] for t in groupFactory.list(type="CellBarcode")] + [""]
groupMap = {}
for s,n,d in sorted(groupFactory.list(type="CellBarcode")):
    groupMap[n] = s

umiOptions = [t[1] for t in groupFactory.list(type="UMI")]
umiMap = {}
for s,n,d in sorted(groupFactory.list(type="UMI")):
    umiMap[n] = s

if not 'umibarcodes' in kwargs and not kwargs.get('umibarcodes',None):
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
outputumis = kwargs.get('outputumis',True)

if verbose:
    print("Parameters:",file=sys.stderr)
    for k,v in sorted(kwargs.items()):
        if k in ('cellbarcodes','umibarcodes') and v:
            print("  %s=%s"%(k,v._name),file=sys.stderr)
        else:
            print("  %s=%s"%(k,v),file=sys.stderr)
    print(file=sys.stderr)

outfile = kwargs.get('outfile')
if 'outfile' in kwargs:
    del kwargs['outfile']
if 'outputumis' in kwargs:
    del kwargs['outputumis']

start = time.time()
locifreq = SCSelectedLociByFetch(bf,snvs,**kwargs)
if not outfile:
    writer = scrctxtwriter()
else:
    assert(outfile.endswith('.brc'))
    writer = scrcbinwriter(outfile)
for chr,pos,refalt,freq in locifreq.loci():
    ref = refalt[0]
    for alt in refalt[1:]:
        for cb in freq:
            allumi = set()
            vals = dict()
            totalreads = 0
            vals1 = dict()
            for nuc in 'ACGT':
                allumi.update(freq[cb][nuc])
                vals[nuc] = len(freq[cb][nuc])
                vals1[nuc] =  sum(freq[cb][nuc].values())         
                totalreads += vals1[nuc]
            if len(allumi) < kwargs.get('minumipercell'):
                continue
            if totalreads < kwargs.get('minreadspercell'):
                continue
            if vals[alt] < kwargs.get('minvarumipercell'):
                continue
            if vals1[alt] < kwargs.get('minvarreadspercell'):
                continue
            if outputumis:
                writer.writecounts(chr,pos,ref,alt,None,cb,vals[ref],vals[alt])
            else:
                writer.writecounts(chr,pos,ref,alt,None,cb,vals1[ref],vals1[alt])
writer.close()
elapsed = time.time() - start
# print("Elapsed: %d:%02d.%02d"%(int(elapsed//60),int(elapsed%60),int(round(100*(elapsed-int(elapsed))))))
