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
from variantloci import VariantLociByFetch

opt = defaultdict(lambda: False)
while len(sys.argv) > 1 and sys.argv[1] in ("-v","-h"):
   opt[sys.argv.pop(1)[1]] = True

if len(sys.argv) <= 1 or opt['h']:
    print("""
Usage:
  varLoci [ options ] <bamfile> [ parameters ]

Parameters (with defaults):
  maxedits=5             # Max. edits for acceptable read alignment
  minbasequal=25         # Min. base quality score
  mindist=5              # Min. number of bases between variants
  minmappingqual=60      # Min. mapping quality for reads
  minvarcnt=5            # Min. number of reads with variant
  region=                # Restriction analysis to chrom. region
  outfile=               # Place output in file, otherwise standard out

Options:
  -v 	Verbose
  -h	Help
""".strip())
    sys.exit(0)

verbose = opt['v']

bamfilename = sys.argv[1]
kwargs = dict(region = None, minvarcnt = 5, maxedits = 5, minbasequal = 25, minmappingqual = 60, mindist = 5, outfile="")
for kv in sys.argv[2:]:
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

scan = VariantLociByFetch(bamfilename,**kwargs)
if not outfile:
    outfile = sys.stdout
else:
    outfile = open(outfile,'wt')
print("CHROM","POS","REF","ALT",sep='\t',file=outfile)
for chrom,pos,refnuc,altnuc,freq in scan.loci():
    print(chrom,pos,refnuc,altnuc,sep='\t')

if outfile != sys.stdout:
    outfile.close()
