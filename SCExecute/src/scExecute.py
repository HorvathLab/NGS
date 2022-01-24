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
from util import ReadGroupFactory
from split import SplitBAM
import multiprocessing

from release import RELEASE, VERSION
VERSION = "%s (%s:%s)"%(VERSION,RELEASE,VERSION)

def excepthook(etype, value, tb):
    traceback.print_exception(etype, value, tb)
    print("Type <Enter> to Exit...", end=' ', file=sys.stderr)
    sys.stderr.flush()
    input()

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

groupFactory = ReadGroupFactory()
groupOptions = [""] + [t[0] for t in groupFactory.list()]
groupDesc = []
for o,d in sorted(groupFactory.list()):
    groupDesc.append("%s (%s)"%(o,d.strip('.')))
readgroup_default = "STARsolo"
threads_default = 1
batch_default = 100

advanced = OptionGroup(parser, "Advanced")
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignment files in (indexed) BAM format. Required.", name="Read Alignment Files",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignment Files (BAM)", "*.bam")])
parser.add_option("-G", "--readgroup", type="choice", dest="readgroup", default=readgroup_default, remember=True,
                  choices=groupOptions, name="Read Group", notNone=True,
                  help="Read group / barcode extraction strategy. Options: %s. Default: %s."%(", ".join(groupDesc),readgroup_default))
parser.add_option("-C", "--command", type="string", dest="command", default=None, notNone=True, remember=True,
                  help="Command to execute for each read-group specific BAM file. The BAM filename replaces {} in the command or placed at the end of the command if no {} is present. Required.", name="Command")
advanced.add_option("-R", "--region", type="str", dest="region", default="", remember=True,
                    help="Restrict reads to those aligning to a specific region. Default: No restriction.", name="Region")
advanced.add_option("-t", "--threads", type="int", dest="threads", default=threads_default, remember=True,
                    help="Number of concurrent instances of command to run. Default: 1.", name="Threads")
advanced.add_option("-B","--batch", type="int", dest="batch", default=batch_default,
                    help="Batch-size per BAM-file pass. Default: 10", name="Batch Size",remember=True)
advanced.add_option("-i", "--index", action="store_true", dest="index", default=False, remember=True,
                    help="Index read-group specific BAM file before executing command.", name="Index")
advanced.add_option("-b","--barcode_acceptlist", type="file", dest="acceptlist", default=None,
                    help="File of white-space separated, acceptable read group values (barcode accept list). Overrides value, if any, specified by Read Group. Use None to remove a default accept list.", name="Valid Read Groups", remember=True,
                    filetypes=[("Valid Read Groups File", "*.txt;*.tsv")])
advanced.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, remember=True,
                    help="Quiet.", name="Quiet")
# advanced.add_option("-d", "--debug", action="store_true", dest="debug", default=False, remember=True,
#                     help="Debug.", name="Debug")
parser.add_option_group(advanced)

region = None
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
        if opt.region:
            if ':' in opt.region:
                contig,rest = opt.region.split(':',1)
                contig = contig.strip()
                start,stop = rest.split('-',1)
                if start.strip() == "":
                    start = None
                else:
                    start = int(start)
                if stop.strip() == "":
                    stop = None
                else:
                    stop = int(stop)
            else:
                contig = opt.region
                start = None
                stop = None
            region = (contig,start,stop)
    except ValueError:
        parser.error("Bad Region option",**error_kwargs)
        continue

    break

opt.debug = False

readgroupparam = ""
if opt.acceptlist != None:
    if opt.acceptlist in ("","None","-"):
        readgroupparam = "*:acceptlist=None"
    else:
        readgroupparam = "*:acceptlist='%s'"%(opt.acceptlist,)
readgroup = groupFactory.get(opt.readgroup,readgroupparam)

progress = ProgressText(quiet=opt.quiet)

doublequote = lambda s: '"%s"'%(s,)
indent = lambda s,n: "\n".join([(" "*n)+l for l in s.splitlines()])

args = []
args.extend(["-r",doublequote(" ".join(opt.alignments))])
if opt.readgroup != readgroup_default:
    args.extend(["-G",doublequote(opt.readgroup)])
args.extend(["-C",doublequote(opt.command)])
if opt.region != "":
    args.extend(["-R",doublequote(opt.region)])
if opt.threads != threads_default:
    args.extend(["-t",str(opt.threads)])
if opt.batch != batch_default:
    args.extend(["-B",str(opt.batch)])
if opt.acceptlist != None:
    args.extend(["-b",doublequote(opt.acceptlist)])
if opt.quiet:
    args.extend(["-q"])

cmdargs = " ".join(args)

execution_log = """
scExecute Options:
  Read Files (-r):            %s
  Command/Script (-C):        %s
  Read Groups (-G):           %s%s

  Advanced:
    Region (-R):              %s
    Threads (-t):             %s
    Batch size (-B):          %s
    Valid Read Groups (-b):   %s
    Quiet (-q):               %s

Command-Line: scExecute %s
"""%(", ".join(opt.alignments),
     opt.command,
     opt.readgroup,
     "\n"+indent(readgroup.tostr(),12),
     opt.region,
     opt.threads,
     opt.batch,
     "" if opt.acceptlist == None else opt.acceptlist,
     opt.quiet,
     cmdargs)

progress.message(execution_log)

execution_queue = multiprocessing.Queue(opt.batch)

def clean(bamfile):
    for fn in (bamfile, bamfile+".bai"):
        if os.path.exists(fn):
            os.unlink(fn)

def execution_worker(execution_queue):
    while True:
        i0,origbam,i1,rg,bamfile = execution_queue.get()
        if bamfile == None:
            break
        bambase = os.path.split(os.path.abspath(origbam))[1].rsplit('.',1)[0]
        subst = {"{}": bamfile, 
                 "{BAMFILE}": bamfile,
                 "{BAMBASE}": bambase,
                 "{BFINDEX}": str(i0),
                 "{RGINDEX}": str(i1),
                 "{CBINDEX}": str(i1),
                 "{READGRP}": rg,
                 "{BARCODE}": rg,
        } 
        command = opt.command
        for k,v in subst.items():
            if k in command:
                command = command.replace(k,v)
        if command == opt.command:
            command = opt.command + " " + bamfile
        progress.message("Executing: %s"%(command,))
        subprocess.run(command,shell=True,check=True)
        clean(bamfile)

threads = []
for i in range(opt.threads):
    t = multiprocessing.Process(target=execution_worker, args=(execution_queue,))
    # t.daemon = True
    t.start()
    threads.append(t)

tmpdir = tempfile.TemporaryDirectory(prefix=".scexec",dir='.')
tmpdirname = tmpdir.name

allrg = dict()
k1 = 0
k = 0
for bamfile in opt.alignments:
    for rg,splitfile in SplitBAM(bamfile, readgroup, opt.batch, tmpdirname, opt.index, region).iterator():
        if rg not in allrg:
            allrg[rg] = k1
            k1 += 1
        execution_queue.put((k,bamfile,allrg[rg],rg,splitfile))
        k += 1
        
for i in range(opt.threads):
    execution_queue.put((None,None,None,None,None))

for t in threads:
    t.join()
