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

groupFactory = ReadGroupFactory()
groupOptions = [""] + [t[0] for t in groupFactory.list()]

minreads_default = 5
maxreads_default = None
tpb_default = 0
filter_default = "Basic"
readgroup_default = "UMI-tools"

advanced = OptionGroup(parser, "Advanced")
parser.add_option("-s", "--snvs", type="files", dest="snvs", default=None,
                  help="Single-Nucleotide-Variant files. Required.", name="SNV Files",
                  notNone=True, remember=True,
                  filetypes=[("SNV Files", "*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignment files in indexed BAM format. Required.", name="Read Alignment Files",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignment Files (indexed BAM)", "*.bam")])
parser.add_option("-f", "--alignmentfilter", type="choice", dest="filter", default=filter_default, remember=True,
                  help="Alignment filtering strategy. Default: Basic.", choices = filterOptions,
                  name="Alignment Filter")
advanced.add_option("-m", "--minreads", type="int", dest="minreads", default=minreads_default, remember=True,
                    help="Minimum number of good reads at SNV locus per alignment file. Default=5.", name="Min. Reads")
advanced.add_option("-M", "--maxreads", type="string", dest="maxreads", default=maxreads_default, remember=True,
                    help="Scale read counts at high-coverage loci to ensure at most this many good reads at SNV locus per alignment file. Values greater than 1 indicate absolute read counts, otherwise the value indicates the coverage distribution percentile. Default=No maximum.", name="Max. Reads")
advanced.add_option("-t", "--threadsperbam", type="int", dest="tpb", default=tpb_default, remember=True,
                    help="Worker threads per alignment file. Indicate no threading with 0. Default=0.", name="Threads/BAM")
advanced.add_option("-G", "--readgroup", type="choice", dest="readgroup", default=readgroup_default, remember=True,
                    choices=groupOptions, name="Read Group",
                    help="Additional read grouping based on read name/identifier strings or BAM-file RG. Default: None, group reads by BAM-file only.")
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

readfilter = filterFactory.get(opt.filter)
if opt.readgroup:
    readgroup = groupFactory.get(opt.readgroup)
else:
    readgroup = None

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
if opt.readgroup != readgroup_default:
    args.extend(["-G",doublequote(opt.readgroup if readgropup != None else "")])
if opt.tpb != tpb_default:
    args.extend(["-t",str(opt.tpb)])
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
    Min. Reads (-m)           %s (applied only to VAF matrix)
    Max. Reads (-M):          %s
    Read Groups (-G):         %s%s
    Threads per BAM (-t):     %s
    Quiet (-q):               %s

Command-Line: scReadCounts %s
"""%(", ".join(opt.snvs),
     ", ".join(opt.alignments),
     opt.filter,
     "" if readfilter == None else "\n"+indent(readfilter.tostr(),10),
     opt.output,
     opt.minreads,
     opt.maxreads,
     None if readgroup == None else opt.readgroup,
     "" if readgroup == None else "\n"+indent(readgroup.tostr(),12),
     opt.tpb,
     opt.quiet,
     cmdargs)

progress.message(execution_log)

args = []
args.extend(["-s"," ".join(opt.snvs)])
args.extend(["-r"," ".join(opt.alignments)])
args.extend(["-f",opt.filter])
args.extend(["-o",opt.output])
args.extend(["-m",0])
if opt.maxreads != maxreads_default:
    args.extend(["-M",opt.maxreads])
if readgroup != None:
    args.extend(["-G",opt.readgroup])
args.extend(["-t",opt.tpb])
if opt.quiet:
    args.extend(["-q"])
args = [ str(x) for x in args ]

if os.path.exists(opt.output):

    progress.message("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    progress.message("Skipping readCounts, output file present.")
    progress.message(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")

else:

    progress.message("\n>>>>>>>>>>>>>>>>>>>>>")
    progress.message("Execute readCounts...")
    progress.message(">>>>>>>>>>>>>>>>>>>>>\n")
    execprog.execute("readCounts",*args)

outbase,extn = opt.output.rsplit('.',1)
outmatrix1 = outbase + '.cnt.matrix.' + extn
outmatrix2 = outbase + '.vaf-m%d.matrix.'%(opt.minreads,) + extn

args = []
args.extend(["-c",opt.output])
args.extend(["-M","Ref;Var"])
args.extend(["-m",0])
if opt.quiet:
    args.extend(["-q"])
args.extend(["-o",outmatrix1])
args = [ str(x) for x in args ]

if os.path.exists(outmatrix1):

    progress.message("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    progress.message("Skipping readCountsMatrix for Ref;Var, output file present.")
    progress.message(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")

else:

    progress.message("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    progress.message("Execute readCountsMatrix for Ref;Var matrix...")
    progress.message(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
    execprog.execute("readCountsMatrix",*args)

args = []
args.extend(["-c",opt.output])
args.extend(["-M","VAF"])
args.extend(["-m",opt.minreads])
if opt.quiet:
    args.extend(["-q"])
args.extend(["-o",outmatrix2])
args = [ str(x) for x in args ]

if os.path.exists(outmatrix2):

    progress.message("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    progress.message("Skipping readCountsMatrix for VAF, output file present.")
    progress.message(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")

else:

    progress.message("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    progress.message("Execute readCountsMatrix for VAF matrix...")
    progress.message(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
    execprog.execute("readCountsMatrix",*args)
