#!/bin/env python3
import sys
import os
import os.path
import glob
import copy
import traceback
import re
import csv
import tempfile
import urllib.request, urllib.parse, urllib.error
import shutil
import atexit
import subprocess
import time
import math
from collections import defaultdict, Counter
from os.path import join, dirname, realpath
try:
    sys.path.append(join(dirname(realpath(__file__)),
                         '..', '..', 'common', 'src'))
except NameError:
    pass
from optparse_gui import OptionParser, OptionGroup, GUI, UserCancelledError, ProgressText
from util import ReadFilterFactory, BadRead, ReadGroupFactory
from fisher import *
from pileups import SerialPileups, ThreadedPileups, MultiprocPileups
from chromreg import ChromLabelRegistry
from operator import itemgetter

from os.path import join, dirname, realpath, split
if getattr(sys, 'frozen', False):
    scriptdir = dirname(realpath(sys.executable))
    if not scriptdir.endswith('/bin') and not scriptdir.endswith('/MacOS'):
        scriptdir = realpath(os.path.join(scriptdir,".."))
    scriptdirs = [scriptdir]
else:
    scriptdir = dirname(realpath(sys.argv[0]))
    scriptdir1 = realpath(join(scriptdir, '..', '..', 'ReadCounts', 'src'))                                                 
    scriptdirs = [scriptdir,scriptdir1]

try:
    scriptextn = "." + os.path.split(sys.argv[0])[1].rsplit('.', 1)[1]
except:
    scriptextn = ""

from execute import Execute                                                                                                 
execprog = Execute(*scriptdirs,extn=scriptextn)     

from release import RELEASE, VERSION                                                                                        
VERSION = "%s (%s:%s)"%(VERSION,RELEASE,VERSION)

def excepthook(etype, value, tb):
    traceback.print_exception(etype, value, tb)
    print("Type <Enter> to Exit...", end=' ', file=sys.stderr)
    sys.stderr.flush()
    input()

toremove = []

def cleanup():
    for d in toremove:
        shutil.rmtree(d, ignore_errors=True)
atexit.register(cleanup)

if len(sys.argv) == 2 and sys.argv[1] == '--GUI':
    from optparse_gui.needswx import *
    sys.exit(1)

if len(sys.argv) == 1:
    if not GUI():
        print("Graphical user-interface unavailable.",file=sys.stderr)
        sys.exit(1)
    from optparse_gui import OptionParserGUI
    parser = OptionParserGUI(version=VERSION)
    error_kwargs = {'exit': False}
    sys.excepthook = excepthook
else:
    parser = OptionParser(version=VERSION)
    error_kwargs = {}

filterFactory = ReadFilterFactory()
filterOptions = [t[0] for t in filterFactory.list()]
filterDesc = []
filterMap = {}
for n,s,d in sorted(filterFactory.list()):
    filterDesc.append("%s (%s)"%(n,d.strip('.')))
    filterMap[n] = s

groupFactory = ReadGroupFactory()
groupOptions = [t[1] for t in groupFactory.list(type="CellBarcode")] + [""]
groupDesc = []
groupMap = {}
for s,n,d in sorted(groupFactory.list(type="CellBarcode")):
    groupDesc.append("%s (%s)"%(n,d.strip('.')))
    groupMap[n] = s
umiOptions = [""] + [t[1] for t in groupFactory.list(type="UMI")]
umiDesc = []
umiMap = {}
for s,n,d in sorted(groupFactory.list(type="UMI")):
    umiDesc.append("%s (%s)"%(n,d.strip('.')))
    umiMap[n] = s

minreads_default = 5
maxreads_default = None
threads_default = 0
filter_default = "Basic"
readgroup_default = "UMI-tools"
umicount_default = None

advanced = OptionGroup(parser, "Advanced")
parser.add_option("-s", "--snvs", type="files", dest="snvs", default=None,
                  help="Single-Nucleotide-Variant files. Required.", name="SNV Files",
                  notNone=True, remember=True,
                  filetypes=[("SNV Files", "*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignment files in indexed BAM format. Required.", name="Read Alignment Files",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignment Files (indexed BAM)", "*.bam")])
parser.add_option("-C", "--cellbarcode", type="choice", dest="cellbarcode", default=readgroup_default, remember=True,
                    choices=groupOptions, name="Cell Barcode",
                    help="Group reads based on cell-barcodes extracted from read name/identifiers or BAM-file tags. Options: %s. Default: %s."%(", ".join(groupDesc),readgroup_default))
parser.add_option("-U", "--umicount", type="choice", dest="umicount", default=umicount_default, remember=True,
                    choices=umiOptions, name="UMI Count",
                    help="Count unique identifiers (UMI) based on read name/identifiers or BAM-file tags. Options: %s. Default: None, count reads not UMIs."%(", ".join(umiDesc),))
parser.add_option("-f", "--alignmentfilter", type="choice", dest="filter", default=filter_default, remember=True,
                  help="Alignment filtering strategy. Options: %s. Default: Basic."%(", ".join(filterDesc),), choices = filterOptions,
                  name="Alignment Filter")
advanced.add_option("-m", "--minreads", type="int", dest="minreads", default=minreads_default, remember=True,
                    help="Minimum number of good reads at SNV locus per alignment file. Default=5.", name="Min. Reads")
advanced.add_option("-M", "--maxreads", type="string", dest="maxreads", default=maxreads_default, remember=True,
                    help="Scale read counts at high-coverage loci to ensure at most this many good reads at SNV locus per alignment file. Values greater than 1 indicate absolute read counts, otherwise the value indicates the coverage distribution percentile. Default=No maximum.", name="Max. Reads")
advanced.add_option("-D", "--directional", action="store_true", dest="directional", default=False, remember=True,
                    help="Output directional (forward and reverse complement) VAF and read counts. Default: False", name="Directional")
advanced.add_option("-b","--barcode_acceptlist", type="file", dest="acceptlist", default=None,
                  help="File of white-space separated, acceptable cell-barcodes. Overrides accept list, if any, specified by Cell Barcode option. Use None to remove a default accept list.", name="Valid Cell Barcodes",
                  remember=True,
                  filetypes=[("Valid Cell Barcodes", "*.txt;*.tsv")])
advanced.add_option("-t", "--threads", type="int", dest="threads", default=threads_default, remember=True,
                    help="Worker threads. Indicate no threading/multiprocessing with 0. Default=0.", name="Threads")
# advanced.add_option("--alignmentfilterparam", type="string", dest="filterparam", default="", remember=True,
#                     help="Override parameters for selected alignment filter. Default: Do not override.", name="Alignment Filter Param.")
# advanced.add_option("--readgroupparam", type="string", dest="readgroupparam", default="", remember=True,
#                     help="Override parameters for selected read group. Default: Do not override.", name="Read Group Param.")
advanced.add_option("-F", "--force", action="store_true", dest="force", default=False, remember=True,
                    help="Force all output files to be re-computed, even if already present. Default: False.", name="Force")
advanced.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, remember=True,
                    help="Quiet.", name="Quiet")
parser.add_option("-o", "--output", type="savefile", dest="output", remember=True,
                  help="Output file. Required.", notNone=True, default=None,
                  name="Output File", filetypes=[("All output formats", "*.xlsx;*.xls;*.csv;*.tsv;*.txt"),
                                                 ("Excel", "*.xlsx"), ("Excel2003", "*.xls"),
                                                 ("CSV", "*.csv"), ("TSV", "*.tsv"), ("Text", "*.txt")])
parser.add_option_group(advanced)

opt = None
while True:
    if 'exit' in error_kwargs:
        try:
            opt, args = parser.parse_args(opts=opt)
        except UserCancelledError:
            sys.exit(0)
    else:
        opt, args = parser.parse_args()

    try:
        if opt.maxreads not in (None,""):
            opt.maxreads = float(opt.maxreads)
        else:
            opt.maxreads = None
    except ValueError:
        parser.error("Bad Max. Read option",**error_kwargs)
        continue

    break

readfilter = filterFactory.get(filterMap[opt.filter])
if opt.cellbarcode not in ("","None","-",None):
    readgroupparam = ""
    if opt.acceptlist != None:
        if opt.acceptlist in ("","None","-"):
            readgroupparam = "*:acceptlist=None"
        else:
            readgroupparam = "*:acceptlist='%s'"%(opt.acceptlist,)
    readgroup = groupFactory.get(groupMap[opt.cellbarcode],readgroupparam)
else:
    readgroup = None

if opt.umicount not in ("","None","-",None):
    umicount = groupFactory.get(umiMap[opt.umicount])
else:
    umicount = None

progress = None
if not opt.output:
    opt.quiet = True
progress = ProgressText(quiet=opt.quiet)

doublequote = lambda s: '"%s"'%(s,)
indent = lambda s,n: "\n".join([(" "*n)+l for l in s.splitlines()])

args = []
args.extend(["-s",doublequote(" ".join(opt.snvs))])
args.extend(["-r",doublequote(" ".join(opt.alignments))])
if opt.filter != filter_default:
    args.extend(["-f",doublequote(opt.filter)])
if opt.minreads != minreads_default:
    args.extend(["-m",str(opt.minreads)])
if opt.maxreads != maxreads_default:
    args.extend(["-M",str(opt.maxreads)])
if opt.directional:
    args.extend(["-D"])
if opt.cellbarcode != readgroup_default:
    args.extend(["-C",doublequote(groupMap[opt.cellbarcode] if readgroup != None else "")])
if opt.acceptlist != None and readgroup != None:
    args.extend(["-b",doublequote(opt.acceptlist)])
if opt.umicount not in ("",None,"None","-"):
    args.extend(["-U",doublequote(opt.umicount if umicount != None else "")])
if opt.threads != threads_default:
    args.extend(["-t",str(opt.threads)])
if opt.force:
    args.extend(["-F"])
if opt.quiet:
    args.extend(["-q"])
args.extend(["-o",doublequote(opt.output)])

cmdargs = " ".join(args)

execution_log = """
scReadCounts Options:
  SNV Files (-s):             %s
  Read Files (-r):            %s
  Read/Alignment Filter (-f): %s%s
  Outfile File (-o):          %s

  Advanced:
    Min. Reads (-m)           %s (applied only to VAF, FVAF, RVAF)
    Max. Reads (-M):          %s
    Directional Counts (-D):  %s
    Cell Barcode (-C):        %s%s
    Valid Cell Barcode (-b):  %s
    UMI Count (-U):           %s%s
    Threads (-t):             %s
    Quiet (-q):               %s

Command-Line: scReadCounts %s
"""%(", ".join(opt.snvs),
     ", ".join(opt.alignments),
     opt.filter,
     "" if readfilter == None else "\n"+indent(readfilter.tostr(),10),
     opt.output,
     opt.minreads,
     opt.maxreads,
     opt.directional,
     None if readgroup == None else opt.cellbarcode,
     "" if readgroup == None else "\n"+indent(readgroup.tostr(),12),
     "" if opt.acceptlist else opt.acceptlist,
     None if umicount == None else opt.umicount,
     "" if umicount == None else "\n"+indent(umicount.tostr(),12),
     opt.threads,
     opt.quiet,
     cmdargs)

progress.message(execution_log)

args = []
args.extend(["-s"," ".join(opt.snvs)])
args.extend(["-r"," ".join(opt.alignments)])
args.extend(["-f",opt.filter])
args.extend(["-o",opt.output])
args.extend(["-m",0])
args.extend(["-B",10000])
if opt.maxreads != maxreads_default:
    args.extend(["-M",opt.maxreads])
if readgroup != None:
    args.extend(["-G",groupMap[opt.cellbarcode]])
    if opt.acceptlist:
        args.extend(["-b",opt.acceptlist])
if umicount != None:
    args.extend(["-U",umiMap[opt.umicount]])
args.extend(["-t",opt.threads])
if opt.quiet:
    args.extend(["-q"])
args = [ str(x) for x in args ]

if os.path.exists(opt.output) and not opt.force:

    progress.message("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    progress.message("Skipping readCounts, output file present.")
    progress.message(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")

else:

    progress.message("\n>>>>>>>>>>>>>>>>>>>>>")
    progress.message("Execute readCounts...")
    progress.message(">>>>>>>>>>>>>>>>>>>>>\n")
    execprog.execute("readCounts",*args)
    opt.force = True

outbase,extn = opt.output.rsplit('.',1)
outmatrix1 = outbase + '.cnt.matrix.' + extn
outmatrix2 = outbase + '.vaf-m%d.matrix.'%(opt.minreads,) + extn

if opt.directional:
    matrix1 = "FRef;FVar;RRef;RVar"
    matrix2 = "FVAF;RVAF"
else:
    matrix1 = "Ref;Var"
    matrix2 = "VAF"

args = []
args.extend(["-c",opt.output])
args.extend(["-M",matrix1])
args.extend(["-m",0])
if opt.quiet:
    args.extend(["-q"])
args.extend(["-o",outmatrix1])
args = [ str(x) for x in args ]

if os.path.exists(outmatrix1) and not opt.force:

    progress.message("\n"+">"*(52+len(matrix1)))
    progress.message("Skipping readCountsMatrix for "+matrix1+", output file present.")
    progress.message(">"*(52+len(matrix1))+"\n")

else:

    progress.message("\n"+">"*(39+len(matrix1)))
    progress.message("Execute readCountsMatrix for "+matrix1+" matrix...")
    progress.message(">"*(39+len(matrix1))+"\n")
    execprog.execute("readCountsMatrix",*args)

args = []
args.extend(["-c",opt.output])
args.extend(["-M",matrix2])
args.extend(["-m",opt.minreads])
if opt.quiet:
    args.extend(["-q"])
args.extend(["-o",outmatrix2])
args = [ str(x) for x in args ]

if os.path.exists(outmatrix2) and not opt.force:

    progress.message("\n"+">"*(52+len(matrix2)))
    progress.message("Skipping readCountsMatrix for "+matrix2+", output file present.")
    progress.message(">"*(52+len(matrix2))+"\n")

else:

    progress.message("\n"+">"*(39+len(matrix2)))
    progress.message("Execute readCountsMatrix for "+matrix2+" matrix...")
    progress.message(">"*(39+len(matrix2))+"\n")
    execprog.execute("readCountsMatrix",*args)
