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
  scCountLoci [ options ] <bamfile> [ <snvlocifile> | - ] [ parameters ]

Parameters (with defaults):
  acceptlist=            # File of cell-barcodes to consider
  cellbarcodes=STARsolo  # Name of cell-barcode extraction strategy
  maxedits=5             # Max. edits for acceptable read alignment
  minbasequal=25         # Min. base quality score
  mincells=3             # Min. number of cells for cell and variant constraints
  minmappingqual=60      # Min. mapping quality for reads
  minreadspercell=0      # Min. number of reads per cell
  minumipercell=3        # Min. number of UMIs per cell
  minvarreadspercell=0   # Min. number of variant reads per cell
  minvarumipercell=3     # Min. number of variant UMIs per cell
  outfile=               # Place output in compact binary file
  outputumis=True        # Output UMI counts, not reads
  outputcells=False      # Output cell counts, not UMI or read counts per cell
  umibarcodes=STARsolo   # Name of UMI extraction strategy

Options:
  -v    Verbose
  -h    Help
""".strip())
    sys.exit(0)

verbose = opt['v']

bf = sys.argv[1]
sf = sys.argv[2]
kwargs = {'minbasequal': 25, 'minmappingqual': 60, 'gapsize': 1000, 'maxedits': 5, 
          'cellbarcodes': 'STARsolo', 'mincells': 3, 'minvarumipercell': 3, 'minumipercell': 3,
          'acceptlist': "", 'minvarreadspercell': 0, 'minreadspercell': 0, 
          'outputumis': True, 'outputcells': False, 'outfile': "", 'umibarcodes': None, }
for kv in sys.argv[3:]:
    badsplit = False
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
if sf == "-":
    snvs = [ l.split()[:4] for l in sys.stdin if not l.startswith('CHR') ]
else:
    snvs = [ l.split()[:4] for l in open(sf) if not l.startswith('CHR') ]

groupFactory = ReadGroupFactory()
groupOptions = [t[1] for t in groupFactory.list(type="CellBarcode")] + [""]
groupMap = {}
for s,n,d in sorted(groupFactory.list(type="CellBarcode")):
    groupMap[n] = s

umiOptions = [t[1] for t in groupFactory.list(type="UMI")]
umiMap = {}
for s,n,d in sorted(groupFactory.list(type="UMI")):
    umiMap[n] = s

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
outputumis = kwargs.get('outputumis',True)
outputcells = kwargs.get('outputcells',False)

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
if 'outputcells' in kwargs:
    del kwargs['outputcells']

start = time.time()
locifreq = SCSelectedLociByFetch(bf,snvs,**kwargs)
if not outputcells:
    if not outfile:
        writer = scrctxtwriter()
    else:
        assert(outfile.endswith('.brc'))
        writer = scrcbinwriter(outfile)
else:
    if not outfile:
        writer = sys.stdout
    else:
        writer = open(outfile,'wt')
    writer.write("\t".join(["CHROM","POS","REF","ALT","SNVCELLS","CELLS"])+"\n")
for chr,pos,refalt,freq in locifreq.loci():
    ref = refalt[0]
    for alt in refalt[1:]:
        cellcnt = 0
        varcellcnt = 0
        refcellcnt = 0
        bothcellcnt = 0
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
            cellcnt += 1
            varbad = False
            if vals[alt] < kwargs.get('minvarumipercell'):
                varbad = True
            elif vals1[alt] < kwargs.get('minvarreadspercell'):
                varbad = True
            if not varbad:
                varcellcnt += 1
            refbad = False
            if vals[ref] < kwargs.get('minvarumipercell'):
                refbad = True
            elif vals1[ref] < kwargs.get('minvarreadspercell'):
                refbad = True
            if not refbad:
                refcellcnt += 1
            if not refbad and not varbad:
                bothcellcnt += 1
            if varbad:
                continue
            if not outputcells:
                if outputumis:
                    writer.writecounts(chr,pos,ref,alt,None,cb,vals[ref],vals[alt])
                else:
                    writer.writecounts(chr,pos,ref,alt,None,cb,vals1[ref],vals1[alt])
        if outputcells:
            writer.write("\t".join(map(str,[chr,pos,ref,alt,varcellcnt,cellcnt]))+"\n")
if not outputcells:
    writer.close()
else:
    if writer != sys.stdout:
        writer.close()
elapsed = time.time() - start
# print("Elapsed: %d:%02d.%02d"%(int(elapsed//60),int(elapsed%60),int(round(100*(elapsed-int(elapsed))))))
