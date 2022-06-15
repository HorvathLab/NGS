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
groupOptions = [t[1] for t in groupFactory.list(type="CellBarcode")]
groupDesc = []
groupMap = {}
for s,n,d in sorted(groupFactory.list(type="CellBarcode")):
    groupDesc.append("%s (%s)"%(n,d.strip('.')))
    groupMap[n] = s

readgroup_default = "STARsolo"

threads_default = 1
batch_default = 100

advanced = OptionGroup(parser, "Advanced")
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignment files in (indexed) BAM format. Required.", name="Read Alignment Files",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignment Files (BAM)", "*.bam")])
parser.add_option("-G", "--cellbarcode", type="choice", dest="readgroup", default=readgroup_default, remember=True,
                  choices=groupOptions, name="Cell Barcode", notNone=True,
                  help="Cell barcode extraction strategy. Options: %s. Default: %s."%(", ".join(groupDesc),readgroup_default))
parser.add_option("-C", "--command", type="string", dest="command", default="", remember=True,
                  help="Command to execute for each read-group specific BAM file. The BAM filename replaces {} in the command or placed at the end of the command if no {} is present. At least one of Command/--command/-C and/or File Template/--filetemplate/-F must be specified.", name="Command")
parser.add_option("-F", "--filetemplate", type="string", dest="filetempl", default="", remember=True,
                  help="Filename template for each read-group specific BAM file. Use {BAMBASE} and {BARCODE} to construct the filename. The read-group specific BAM file should end in \".bam\" and will not be deleted after the command, if specified, is run. At least one of Command/--command/-C and/or File Template/--filetemplate/-F must be specified.", name="File Template")
advanced.add_option("-D", "--directory", type="str", dest="directory", default="", remember=True,
                    help="Working directory for running command on each read-group specific BAM file. Default: Current working directory.", name="Directory Template")
advanced.add_option("-o", "--outtemplate", type="string", dest="stdouttempl", default="", remember=True,
                  help="Filename template for the standard output of each read-group specific command execution. Use {BAMBASE} and {BARCODE} to construct the filename. Default: Standard Output of command not captured.", name="Output Template")
advanced.add_option("-e", "--errtemplate", type="string", dest="stderrtempl", default="", remember=True,
                  help="Filename template for the standard error of each read-group specific command execution. Use {BAMBASE} and {BARCODE} to construct the filename. Default: Standard error of command not captured.", name="Error Template")
advanced.add_option("-O", "--allouttemplate", type="string", dest="allouttempl", default="", remember=True,
                  help="Filename template for the standard output and standard error of each read-group specific command execution. Use {BAMBASE} and {BARCODE} to construct the filename. Default: Output/Error of command not captured.", name="All Output Template")
advanced.add_option("-L", "--limit", type="string", dest="limit", default="", remember=True,
                    help="Generate at most this many read-group specific BAM files. Default: No limit.", name="Limit")
advanced.add_option("-R", "--region", type="str", dest="region", default="", remember=True,
                    help="Restrict reads to those aligning to a specific region specified as chrom:start-end. Default: No restriction.", name="Region")
advanced.add_option("--regions", type="file", dest="regions", default="", remember=True,
                    help="Restrict reads to those aligning to specific region(s). If a GTF-format genome annotation file (extn *.gtf) is provided, restrict to genic regions; otherwise, one region per line specified as chrom:start-end in a file with extn *.txt. Default: No restriction.", name="Region File", filetypes=[("Genome Annotation", "*.gtf"),("Region List", "*.txt")])
advanced.add_option("-t", "--threads", type="int", dest="threads", default=threads_default, remember=True,
                    help="Number of concurrent instances of command to run. Default: 1.", name="Threads")
advanced.add_option("--cpuaffinity", action="store_true", dest="affinity", default=False, remember=True,
                    help="Constrain each instance of a command to a single CPU. Default: False.", name="CPU Affinity")
advanced.add_option("-B","--batch", type="int", dest="batch", default=batch_default,
                    help="Batch-size per BAM-file pass. Default: 10", name="Batch Size",remember=True)
advanced.add_option("-i", "--index", action="store_true", dest="index", default=False, remember=True,
                    help="Index read-group specific BAM file before executing command.", name="Index")
advanced.add_option("-b","--barcode_acceptlist", type="file", dest="acceptlist", default=None,
                    help="File of white-space separated, acceptable cell barcode values. Overrides value, if any, specified by Cell Barcode. Use None to remove a default accept list.", name="Valid Cell Barcodes", remember=True,
                    filetypes=[("Valid Cell Barcodes File", "*.txt;*.tsv")])
advanced.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, remember=True,
                    help="Quiet.", name="Quiet")
# advanced.add_option("-d", "--debug", action="store_true", dest="debug", default=False, remember=True,
#                     help="Debug.", name="Debug")
parser.add_option_group(advanced)


def extract_region(regionstr):
    if ':' in regionstr:
        contig,rest = regionstr.split(':',1)
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
    return (contig,start,stop)

regions = None
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
            regions = [extract_region(opt.region.strip())]
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
readgroup = groupFactory.get(groupMap[opt.readgroup],readgroupparam)
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
    args.extend(["-e",doublequote(opt.stderrtempl)])
if opt.allouttempl != "":
    args.extend(["-O",doublequote(opt.allouttempl)])
if opt.limit not in ("",None):
    args.extend(["-L",str(opt.limit)])
if opt.region != "":
    args.extend(["-R",doublequote(opt.region)])
if opt.regions not in ("",None):
    args.extend(["--regions",doublequote(opt.regions)])
if opt.threads != threads_default:
    args.extend(["-t",str(opt.threads)])
if opt.affinity:
    args.extend(["--cpuaffinity",])
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
    Directory Template (-D):  %s
    Output Template (-o):     %s
    Error Template (-e):      %s
    All Output Template (-O): %s
    Limit (-L):               %s
    Region (-R):              %s
    Region List (--regions):  %s
    Threads (-t):             %s
    Affinity (--cpuaffinity): %s
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
     opt.regions if opt.regions not in ("",None) else "",
     opt.threads,
     str(opt.affinity),
     opt.batch,
     str(opt.index),
     "" if opt.acceptlist == None else opt.acceptlist,
     opt.quiet,
     cmdargs)

progress.message(execution_log)

starttime = time.time()

if opt.regions:
    progress.message("Reading regions...")
    regions = []
    try:
        if opt.regions.lower().endswith('.gtf'):
            for l in open(opt.regions):
                if l.strip() == "" or l.startswith('#'):
                    continue
                sl = l.split()
                if sl[2] == "gene":
                    chrom = sl[0]
                    start = int(sl[3])
                    stop = int(sl[4])
                    regions.append((chrom,start,stop))
        elif opt.regions.lower().endswith('.txt'):
            for l in open(opt.regions):
                if l.strip() == "" or l.startswith('#'):
                    continue
                regions.append(extract_region(l.strip()))
    except (ValueError,IndexError,IOError):
        progress.message("Bad Region List")
        sys.exit(1)
    if len(regions) == 0:
        progress.message("No regions selected!")
        sys.exit(1)
    progress.message("Restriction to %d genomic regions"%(len(regions),))

execution_queue = multiprocessing.Queue(opt.batch)
done_queue = multiprocessing.Queue()
output_lock = multiprocessing.Lock()
totaltime = multiprocessing.Value('d',0.0)
totaljobs = multiprocessing.Value('I',0)

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

def execution_worker(index,execution_queue,done_queue,output_lock):
    if opt.affinity:
        try:
            cpus = sorted(os.sched_getaffinity(0))
            os.sched_setaffinity(0, {cpus[index%len(cpus)]})
        except:
            pass
    while True:
        i0,origbam,i1,rg,bampath,rgind,rgtot,bamtot = execution_queue.get()
        if bampath == None:
            done_queue.put((None,None,None,None,None))
            break
        if opt.directory != "":
            bampath = os.path.realpath(bampath)
        bamfile = os.path.split(bampath)[1]
        bambase = bamfile.rsplit('.',1)[0]
        if opt.directory != "":
            origbam = os.path.realpath(origbam)
        orifile = os.path.split(origbam)[1]
        oribase = orifile.rsplit('.',1)[0]
        subst = {"{}": bampath, 
                 "{CBPATH}": bampath, # Full path to generated BAM file
                 "{CBFILE}": bamfile, # Filename part of generated BAM file
                 "{CBBASE}": bambase, # Filename of generated BAM file, no .bam
                 "{RGPATH}": bampath, # Full path to generated BAM file
                 "{RGFILE}": bamfile, # Filename part of generated BAM file
                 "{RGBASE}": bambase, # Filename of gerated BAM file, no .bam
                 "{BAMPATH}": origbam, # Input/Original BAM, full path
                 "{BAMFILE}": orifile, # Input/Original BAM, filename part
                 "{BAMBASE}": oribase, # Input/Original BAM, filename part, no .bam
                 "{BFINDEX}": str(i0), # Original BAM file index
                 "{RGINDEX}": str(i1), # Read Group index
                 "{CBINDEX}": str(i1), # Cell-Barcode index
                 "{READGRP}": rg, # Read Group 
                 "{BARCODE}": rg, # Cell-Barcode
                 "{SCEXEC_WORKER_INDEX}": str(index), # Index of this worker
                 "{WORKER}": str(index), # Index of this worker
        } 
        theworkingdir = "."
        if opt.directory != "":
            theworkingdir = substtempl(opt.directory,subst)
            theworkingdir = os.path.realpath(theworkingdir)
            if not os.path.exists(theworkingdir):
                os.makedirs(theworkingdir)
        thebamfile = bampath
        if opt.filetempl != "":
            thebamfile = substtempl(opt.filetempl,subst)
            if not thebamfile.startswith(os.sep) and opt.directory != "":
                thebamfile = os.path.realpath(os.path.join(theworkingdir,thebamfile))
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
            if not thestderrfile.startswith(os.sep) and opt.directory != "":
                thestderrfile = os.path.realpath(os.path.join(theworkingdir,thestderrfile))
        thestdoutfile = None
        if opt.stdouttempl != "":
            thestdoutfile = substtempl(opt.stdouttempl,subst)
            if not thestdoutfile.startswith(os.sep) and opt.directory != "":
                thestdoutfile = os.path.realpath(os.path.join(theworkingdir,thestdoutfile))
        thealloutfile = None
        if opt.allouttempl != "":
            thealloutfile = substtempl(opt.allouttempl,subst)
            if not thealloutfile.startswith(os.sep) and opt.directory != "":
                thealloutfile = os.path.realpath(os.path.join(theworkingdir,thealloutfile))                                  
        if opt.command != "":
            subst['{}'] = thebamfile
            thecommand = substtempl(opt.command,subst)
            if thecommand == opt.command:
                thecommand = opt.command + " " + thebamfile
            prog,rest = thecommand.split(None,1)
            if os.path.exists(prog) and not prog.startswith(os.sep) and opt.directory != "":
                prog = os.path.realpath(os.path.join(os.getcwd(),prog))
                thecommand = prog + " " + rest
            output_lock.acquire()
            try:
                progress.message("Executing [%d/%d]: %s"%(rgind+1,rgtot,thecommand,))
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
                                        stdout=stdout,stderr=stderr,
                                        env=dict(os.environ,SCEXEC_WORKER_INDEX=str(index)))
                elapsed = time.time()-start
                done_queue.put((i0,bamtot,rgind,rgtot,elapsed))
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

def timestr(seconds):
    secs = seconds%60
    if seconds < 60:
        return "%.2f sec"%(secs,)
    secs = int(secs)
    mins = int(seconds/60)
    if mins < 60:
        return "%d:%02d min:sec"%(mins,secs,)
    mins = int(mins%60)
    hrs = int(mins/60)
    return "%d:%02d hrs:min"%(hrs,mins,)
    
def done_worker(done_queue,output_lock,totaljobs,totaltime,starttime):
    donecount = 0
    count = 0
    sumtime = 0
    done=defaultdict(set)
    rgtots=defaultdict(int)
    while True:
        i0,bamtot,rgind,rgtot,elapsed = done_queue.get()
        if i0 == None:
            donecount += 1
            if donecount == opt.threads:
                break
            else:
                continue
        if i0 not in rgtots:
            rgtots[i0] = rgtot
            avergtot = sum(rgtots.values())/len(rgtots)
            rgtotest = int(round(sum(rgtots.values()) + avergtot*(len(rgtots)-bamtot)))
        if rgind not in done[i0]:
            done[i0].add(rgind)
            count += 1
            totaltime.value = (totaltime.value + elapsed)
            totaljobs.value = (totaljobs.value + 1)
            sumtime = (time.time()-starttime)
            avetime = sumtime/count
        remaining = (rgtotest-count)*avetime
        output_lock.acquire()
        try:
            progress.message("Completed [%d/%d]: Runtime %s, est. time to finish %s, %.2f%% complete."%(rgind+1,rgtot,timestr(elapsed),timestr(remaining),100*(count)/rgtotest))
        finally:
            output_lock.release()

threads = []
for i in range(opt.threads):
    t = multiprocessing.Process(target=execution_worker, args=(i,execution_queue,done_queue,output_lock))
    # t.daemon = True
    t.start()
    threads.append(t)

tmpdir = tempfile.TemporaryDirectory(prefix=".scexec",dir='.')
tmpdirname = tmpdir.name

allrg = dict()
k1 = 0
k = 0
made_done_worker = False
for bamfile in opt.alignments:
    for rg,splitfile,rgind,rgtot in SplitBAM(bamfile, readgroup, opt.batch, tmpdirname, opt.index, regions, limit).iterator():
        if rg not in allrg:
            allrg[rg] = k1
            k1 += 1
        execution_queue.put((k,bamfile,allrg[rg],rg,splitfile,rgind,rgtot,len(opt.alignments)))
        if not made_done_worker:
            output_lock.acquire()
            try:
                progress.message("Time to extract first batch of SCbams: %s"%(timestr(time.time()-starttime,)))
            finally:
                output_lock.release()
            t = multiprocessing.Process(target=done_worker, args=(done_queue,output_lock,totaljobs,totaltime,time.time()))
            t.start()
            threads.append(t)
            made_done_worker = True
    k += 1
        
for i in range(opt.threads):
    execution_queue.put((None,None,None,None,None,None,None,None))

for t in threads:
    t.join()

elapsed = time.time()-starttime
totalcount = totaljobs.value
avejobruntime = (totaltime.value/totalcount)
totaljobtime = math.ceil(totalcount/opt.threads)*avejobruntime
progress.message("Final: Runtime %s, %d jobs, ave. job runtime %s, scExecute overhead %.2f%%."%(timestr(elapsed),totalcount,timestr(avejobruntime),100*(elapsed-totaljobtime)/elapsed))

