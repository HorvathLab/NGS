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
parser.add_option("-C", "--command", type="string", dest="command", default="", remember=True,
                  help="Command to execute for each read-group specific BAM file. The BAM filename replaces {} in the command or placed at the end of the command if no {} is present. At least one of Command/--command/-C and/or File Template/--filetemplate/-F must be specified.", name="Command")
parser.add_option("-F", "--filetemplate", type="string", dest="filetempl", default="", remember=True,
                  help="Filename template for each read-group specific BAM file. Use {ORIBASE} and {BARCODE} to construct the filename, must end in \".bam\". At least one of Command/--command/-C and/or File Template/--filetemplate/-F must be specified.", name="File Template")
advanced.add_option("-D", "--directory", type="str", dest="directory", default="", remember=True,
                    help="Working directory for running command on each read-group specific BAM file. Default: Current working directory.", name="Directory")
advanced.add_option("-o", "--outtemplate", type="string", dest="stdouttempl", default="", remember=True,
                  help="Filename template for the standard output of each read-group specific command execution. Use {ORIBASE} and {BARCODE} to construct the filename. Default: Standard Output of command not captured.", name="Output Template")
advanced.add_option("-l", "--logtemplate", type="string", dest="stderrtempl", default="", remember=True,
                  help="Filename template for the standard error of each read-group specific command execution. Use {ORIBASE} and {BARCODE} to construct the filename. Default: Standard error of command not captured.", name="Error Template")
advanced.add_option("-O", "--allouttemplate", type="string", dest="allouttempl", default="", remember=True,
                  help="Filename template for the standard output and standard error of each read-group specific command execution. Use {ORIBASE} and {BARCODE} to construct the filename. Default: Output/Error of command not captured.", name="All Output Template")
advanced.add_option("-L", "--limit", type="string", dest="limit", default="", remember=True,
                    help="Generate at most this many read-group specific BAM files. Default: No limit.", name="Limit")
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
limit = None
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

    if not opt.command.strip() and not opt.filetempl.strip():
        parser.error("One of Command/File Template must be specified",**error_kwargs)
        continue

    try:
       if opt.limit.strip() not in ("",None):
           limit = int(opt.limit)
       else:
           limit = None
    except ValueError:
        parser.error("Bad Limit option",**error_kwargs)
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
if opt.filetempl:
    opt.filetempl = opt.filetempl.strip()
if opt.command:
    opt.command = opt.command.strip()
if opt.directory:
    opt.directory = opt.directory.strip()
if opt.allouttempl:
    opt.allouttempl = opt.allouttempl.strip()
if opt.stderrtempl:
    opt.stderrtempl = opt.stderrtempl.strip()
if opt.stdouttempl:
    opt.stdouttempl = opt.stdouttempl.strip()

progress = ProgressText(quiet=opt.quiet)

doublequote = lambda s: '"%s"'%(s,)
indent = lambda s,n: "\n".join([(" "*n)+l for l in s.splitlines()])

args = []
args.extend(["-r",doublequote(" ".join(opt.alignments))])
if opt.readgroup != readgroup_default:
    args.extend(["-G",doublequote(opt.readgroup)])
if opt.command != "":
    args.extend(["-C",doublequote(opt.command)])
if opt.filetempl != "":
    args.extend(["-F",doublequote(opt.filetempl)])
if opt.directory != "":
    args.extend(["-D",doublequote(opt.directory)])
if opt.stdouttempl != "":
    args.extend(["-o",doublequote(opt.stdouttempl)])
if opt.stderrtempl != "":
    args.extend(["-l",doublequote(opt.stderrtempl)])
if opt.allouttempl != "":
    args.extend(["-O",doublequote(opt.allouttempl)])
if opt.limit not in ("",None):
    args.extend(["-L",str(opt.limit)])
if opt.region != "":
    args.extend(["-R",doublequote(opt.region)])
if opt.threads != threads_default:
    args.extend(["-t",str(opt.threads)])
if opt.batch != batch_default:
    args.extend(["-B",str(opt.batch)])
if opt.index:
    args.extend(["-i",])
if opt.acceptlist != None:
    args.extend(["-b",doublequote(opt.acceptlist)])
if opt.quiet:
    args.extend(["-q"])

cmdargs = " ".join(args)

execution_log = """
scExecute Options:
  Read Files (-r):            %s
  Read Groups (-G):           %s%s
  Command/Script (-C):        %s
  File Template (-F):         %s

  Advanced:
    Directory (-D):           %s
    Output Template (-o):     %s
    Error Template (-l):      %s
    All Output Template (-O): %s
    Limit (-L):               %s
    Region (-R):              %s
    Threads (-t):             %s
    Batch size (-B):          %s
    Index BAMs (-i):          %s
    Valid Read Groups (-b):   %s
    Quiet (-q):               %s

Command-Line: scExecute %s
"""%(", ".join(opt.alignments),
     opt.readgroup,
     "\n"+indent(readgroup.tostr(),12),
     opt.command,
     opt.filetempl,
     opt.directory,
     opt.stdouttempl,
     opt.stderrtempl,
     opt.allouttempl,
     opt.limit if opt.limit not in ("",None) else "",
     opt.region,
     opt.threads,
     opt.batch,
     str(opt.index),
     "" if opt.acceptlist == None else opt.acceptlist,
     opt.quiet,
     cmdargs)

progress.message(execution_log)

execution_queue = multiprocessing.Queue(opt.batch)
output_lock = multiprocessing.Lock()

def clean(bamfile):
    for fn in (bamfile, bamfile+".bai"):
        if os.path.exists(fn):
            os.unlink(fn)

def substtempl(templ,subst):
    result = templ
    for k,v in subst.items():
        if k in result:
            result = result.replace(k,v)
    return result

def execution_worker(execution_queue,output_lock):
    while True:
        i0,origbam,i1,rg,bampath,rgind,rgtot,bamtot = execution_queue.get()
        if bampath == None:
            break
        bamfile = os.path.split(os.path.abspath(bampath))[1]
        bambase = bamfile.rsplit('.',1)[0]
        orifile = os.path.split(os.path.abspath(origbam))[1]
        oribase = orifile.rsplit('.',1)[0]
        subst = {"{}": bampath, 
                 "{CBPATH}": bampath, # Full path to generated BAM file
                 "{CBFILE}": bamfile, # Filename part of generated BAM file
                 "{CBBASE}": bambase, # Filename of generated BAM file, no .bam
                 "{RGPATH}": bampath, # Full path to generated BAM file
                 "{RGFILE}": bamfile, # Filename part of generated BAM file
                 "{RGBASE}": bambase, # Filename of gerated BAM file, no .bam
                 "{ORIPATH}": origbam, # Input/Original BAM, full path
                 "{ORIFILE}": orifile, # Input/Original BAM, filename part
                 "{ORIBASE}": oribase, # Input/Original BAM, filename part, no .bam
                 "{BFINDEX}": str(i0), # Original BAM file index
                 "{RGINDEX}": str(i1), # Read Group index
                 "{CBINDEX}": str(i1), # Cell-Barcode index
                 "{READGRP}": rg, # Read Group 
                 "{BARCODE}": rg, # Cell-Barcode
        } 
        thebamfile = bampath
        if opt.filetempl != "":
            thebamfile = substtempl(opt.filetempl,subst)
            if opt.command == "":
                output_lock.acquire()
                try:
                    progress.message("Filename: %s"%(thebamfile,))
                finally:
                    output_lock.release()
            shutil.copyfile(bampath,thebamfile)
            if os.path.exists(bampath+'.bai'):
                shutil.copyfile(bampath+'.bai',thebamfile+'.bai')
            clean(bampath)
        thestderrfile = None
        if opt.stderrtempl != "":
            thestderrfile = substtempl(opt.stderrtempl,subst)
        thestdoutfile = None
        if opt.stdouttempl != "":
            thestdoutfile = substtempl(opt.stdouttempl,subst)
        thealloutfile = None
        if opt.allouttempl != "":
            thealloutfile = substtempl(opt.allouttempl,subst)
        theworkingdir = "."
        if opt.directory != "":
            theworkingdir = substtempl(opt.directory,subst)
        if opt.command != "":
            subst['{}'] = thebamfile
            thecommand = substtempl(opt.command,subst)
            if thecommand == opt.command:
                thecommand = opt.command + " " + thebamfile
            output_lock.acquire()
            try:
                progress.message("Executing [%*d/%*d]: %s"%(int(math.log(rgtot,10)),rgind+1,int(math.log(rgtot,10)),rgtot,thecommand,))
            finally:
                output_lock.release()
            try:
                stdout=None; stderr=None
                if opt.stderrtempl != "":
                    stderr=open(thestderrfile,'w')
                if opt.stdouttempl != "":
                    stdout=open(thestdoutfile,'w')
                if opt.allouttempl != "":
                    stdout=open(thealloutfile,'w')
                    stderr=subprocess.STDOUT
                result = None
                start = time.time()
                result = subprocess.run(thecommand,shell=True,check=True,
                                        cwd=theworkingdir,stdin=subprocess.DEVNULL,
                                        stdout=stdout,stderr=stderr)
                elapsed = time.time()-start
                output_lock.acquire()
                try:
                    progress.message("Complete  [%*d/%*d]: %s:%s elapsed."%(int(math.log(rgtot,10)),rgind+1,int(math.log(rgtot,10)),rgtot,int(elapsed/60),int(elapsed%60)))
                finally:
                    output_lock.release()
            except subprocess.SubprocessError as e: 
                output_lock.acquire()
                try:
                    progress.message("Command \"%s\" failed, return code %s."%(thecommand,e.returncode,))
                finally:
                    output_lock.release()
                    break
            finally:
                if stdout != None:
                    stdout.close()
                if stderr != None and stderr != subprocess.STDOUT:
                    stderr.close()
                clean(bampath)

threads = []
for i in range(opt.threads):
    t = multiprocessing.Process(target=execution_worker, args=(execution_queue,output_lock))
    # t.daemon = True
    t.start()
    threads.append(t)

tmpdir = tempfile.TemporaryDirectory(prefix=".scexec",dir='.')
tmpdirname = tmpdir.name

allrg = dict()
k1 = 0
k = 0
for bamfile in opt.alignments:
    for rg,splitfile,rgind,rgtot in SplitBAM(bamfile, readgroup, opt.batch, tmpdirname, opt.index, region, limit).iterator():
        if rg not in allrg:
            allrg[rg] = k1
            k1 += 1
        execution_queue.put((k,bamfile,allrg[rg],rg,splitfile,rgind,rgtot,len(opt.alignments)))
    k += 1
        
for i in range(opt.threads):
    execution_queue.put((None,None,None,None,None,None,None,None))

for t in threads:
    t.join()
