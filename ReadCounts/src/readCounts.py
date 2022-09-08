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
for s,n,d in sorted(filterFactory.list()):
    filterDesc.append("%s (%s)"%(s,d.strip('.')))

groupFactory = ReadGroupFactory()
groupOptions = [""] + [t[0] for t in groupFactory.list()]
groupDesc = []
for s,n,d in sorted(groupFactory.list()):
    groupDesc.append("%s (%s)"%(s,d.strip('.')))

minreads_default = 10
maxreads_default = None
threads_default = 0
filter_default = "Basic"

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
                  help="Alignment filtering strategy. Options: %s. Default: Basic."%(", ".join(filterDesc),), choices = filterOptions,
                  name="Alignment Filter")
advanced.add_option("-m", "--minreads", type="int", dest="minreads", default=minreads_default, remember=True,
                    help="Minimum number of good reads at SNV locus per alignment file. Default=10.", name="Min. Reads")
advanced.add_option("-M", "--maxreads", type="string", dest="maxreads", default=maxreads_default, remember=True,
                    help="Scale read counts at high-coverage loci to ensure at most this many good reads at SNV locus per alignment file. Values greater than 1 indicate absolute read counts, otherwise the value indicates the coverage distribution percentile. Default=No maximum.", name="Max. Reads")
advanced.add_option("-G", "--readgroup", type="choice", dest="readgroup", default=None, remember=True,
                    choices=groupOptions, name="Read Group",
                    help="Additional read grouping based on read name/identifier strings or BAM-file tags. Options: %s. Default: None, group reads by BAM-file only."%(", ".join(groupDesc),))
advanced.add_option("-b","--barcode_acceptlist", type="file", dest="acceptlist", default=None,
                  help="File of white-space separated, acceptable read group values (barcode accept list). Overrides value, if any, specified by Read Group. Use None to remove a default accept list.", name="Valid Read Groups",
                  remember=True,
                  filetypes=[("Valid Read Groups File", "*.txt;*.tsv")])
advanced.add_option("-U", "--umicount", type="choice", dest="umigroup", default=None, remember=True,
                    choices=groupOptions, name="UMI Count",
                    help="Count unique identifiers (UMI) based on read name/identifier strings or BAM-file tags. Options: %s. Default: None, count reads not UMIs."%(", ".join(groupDesc),))
advanced.add_option("-E","--extended",type="multichoice",dest="extended",default=None, remember=True, 
               help="Generate extended output, one or more comma-separated values: Genotype likelihood, Read filtering statistics. Default: No extended ouptut.", name="Extended Output", multichoices=["Genotype likelihood","Read filtering statistics"])
# advanced.add_option("--alignmentfilterparam", type="string", dest="filterparam", default="", remember=True,
#                     help="Override parameters for selected alignment filter. Default: Do not override.", name="Alignment Filter Param.")
# advanced.add_option("--readgroupparam", type="string", dest="readgroupparam", default="", remember=True,
#                     help="Override parameters for selected read group. Default: Do not override.", name="Read Group Param.")
advanced.add_option("-B", "--snvbatchsize", type="int", dest="snvbatchsize", default=None, remember=True,
                    help="Manage memory footprint by making multiple passes through the BAM file, one for each batch of SNVs. Default=All SNVs (single pass).", name="SNV Batch Size")
advanced.add_option("-t", "--threads", type="int", dest="threads", default=threads_default, remember=True,
                    help="Worker threads. Indicate no threading/multiprocessing with 0. Default=0.", name="Threads")
advanced.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, remember=True,
                    help="Quiet.", name="Quiet")
# advanced.add_option("-d", "--debug", action="store_true", dest="debug", default=False, remember=True,
#                     help="Debug.", name="Debug")
parser.add_option("-o", "--output", type="savefile", dest="output", remember=True,
                  help="Output file. Leave empty for console ouptut.", default="",
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

opt.debug = False
readfilter = filterFactory.get(opt.filter)
if opt.readgroup:
    readgroupparam = ""
    if opt.acceptlist != None:
        if opt.acceptlist in ("","None","-"):
            readgroupparam = "*:acceptlist=None"
        else:
            readgroupparam = "*:acceptlist='%s'"%(opt.acceptlist,)
    readgroup = groupFactory.get(opt.readgroup,readgroupparam)
else:
    readgroup = None
if opt.umigroup:
    umigroup = groupFactory.get(opt.umigroup)
else:
    umigroup = None
    
if opt.extended and isinstance(opt.extended,str):
    opt.extended = [ e.strip() for e in opt.extended.split(',') ]

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
if opt.snvbatchsize != None:
    args.extend(["-B",str(opt.snvbatchsize)])
if readgroup != None:
    args.extend(["-G",doublequote(opt.readgroup)])
    if opt.acceptlist != None:
        args.extend(["-b",doublequote(opt.acceptlist)])
if umigroup != None:
    args.extend(["-U",doublequote(opt.umigroup)])
if opt.threads != threads_default:
    args.extend(["-t",str(opt.threads)])
if opt.extended:
    args.extend(["-E",doublequote(",".join(opt.extended))])
if opt.quiet:
    args.extend(["-q"])
if opt.output:
    args.extend(["-o",doublequote(opt.output)])

cmdargs = " ".join(args)

execution_log = """
readCounts Options:
  SNV Files (-s):             %s
  Read Files (-r):            %s
  Read/Alignment Filter (-f): %s%s
  Outfile File (-o):          %s

  Advanced:
    Min. Reads (-m)           %s
    Max. Reads (-M):          %s
    SNV Batch Size (-B):      %s
    Read Groups (-G):         %s%s
    Valid Read Groups (-b):   %s
    UMI Count (-U):           %s%s
    Threads (-t):             %s
    Extended Output (-E):     %s
    Quiet (-q):               %s

Command-Line: readCounts %s
"""%(", ".join(opt.snvs),
     ", ".join(opt.alignments),
     opt.filter,
     "" if readfilter == None else "\n"+indent(readfilter.tostr(),10),
     opt.output,
     opt.minreads,
     opt.maxreads,
     "" if opt.snvbatchsize == None else opt.snvbatchsize,
     None if readgroup == None else opt.readgroup,
     "" if readgroup == None else "\n"+indent(readgroup.tostr(),12),
     "" if opt.acceptlist == None else opt.acceptlist,
     None if umigroup == None else opt.umigroup,
     "" if umigroup == None else "\n"+indent(umigroup.tostr(),12),
     opt.threads,
     ", ".join(opt.extended) if opt.extended else "None",
     opt.quiet,
     cmdargs)

progress.message(execution_log)

if opt.maxreads == None:
    opt.maxreads = 1e+20

def print_memory(msg=""):
  import resource
  rusage = resource.getrusage(resource.RUSAGE_SELF)
  if msg:
      msg = "["+msg+"] "
  print("%sMemory: %.2fGb"%(msg,rusage.ru_maxrss/(1024**2)),flush=True)

from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable, BEDFile, VCFFile

progress.stage("Read SNV data", len(opt.snvs))
snvheaders = [_f for _f in """
CHROM POS REF ALT
""".split() if _f]

snvdata = {}
snvchroms = defaultdict(set)
for filename in opt.snvs:

    base, extn = filename.rsplit('.', 1)
    extn = extn.lower()
    if extn == 'csv':
        snvs = CSVFileTable(filename=filename)
    elif extn == 'vcf':
        snvs = VCFFile(filename=filename)
    elif extn == 'tsv':
        snvs = TSVFileTable(filename=filename)
    elif extn == 'xls':
        snvs = XLSFileTable(filename=filename)
    elif extn == 'xlsx':
        snvs = XLSXFileTable(filename=filename)
    elif extn == 'txt':
        snvs = TXTFileTable(filename=filename, headers=snvheaders)
    else:
        raise RuntimeError("Unexpected SNV file extension: %s" % filename)

    for h in snvheaders:
        if h not in snvs.headers():
            raise RuntimeError(
                "Required header: %s missing from SNV file %s" % (h, filename))

    for h in snvs.headers():
        if h in snvheaders:
            continue

    for r in snvs:
        try:
            chr = str(int(float(r[snvheaders[0]])))
        except ValueError:
            chr = r[snvheaders[0]].strip()
        snvchroms[filename].add(chr)
        locus = int(float(r[snvheaders[1]]))
        ref = r[snvheaders[2]].strip()
        alt = r[snvheaders[3]].strip()
        if r.get('INFO:INDEL'):
            continue
        if len(ref) != 1:
            continue
        if not re.search(r'^[ACGT](,[ACGT])*$', alt):
            continue
        # for h in r:
        #     if r.get(h):
        #         usedsnvheaders.add(h)
        snvkey = (filename, chr, locus, ref, alt)
        if snvkey not in snvdata:
            snvdata[snvkey] = r

    progress.update()
progress.done()

chrreg = ChromLabelRegistry()

for snvfile in snvchroms:
    chrreg.add_labels(snvfile,snvchroms[snvfile])

snvdata1 = {}
for (sf, chr, locus, ref, alt), r in snvdata.items():
    chrom = chrreg.label2chrom(sf,chr)
    assert(chrom)
    snvkey = (chrom,locus,ref,alt)
    if snvkey not in snvdata1:
        snvdata1[snvkey] = (chrom,locus,ref,alt,r)
    # print(snvkey,snvdata1[snvkey])

for bamfile in opt.alignments:
    chrreg.add_bamlabels(bamfile)

chrreg.determine_chrom_order()
# chrreg.default_chrom_order()

snvdata = sorted(list(snvdata1.values()),key=lambda s: (chrreg.chrom_order(s[0]),s[1],s[2],s[3]))
progress.message("SNVs: %d" % len(snvdata))

outheaders = snvheaders + [_f for _f in """
SNVCountForward
SNVCountReverse
RefCountForward
RefCountReverse
SNVCount
RefCount
GoodReads
%BadRead
VAF
""".split() if _f]

dogl = True
if opt.extended == None or 'Genotype likelihood' not in opt.extended:
    dogl = False

glheaders = [_f for _f in """
HomoVarSc
HetSc
HomoRefSc
VarDomSc
RefDomSc
NotHomoVarpV
NotHomoRefpV
NotHetpV
VarDompV
RefDompV
NotHomoVarFDR
NotHomoRefFDR
NotHetFDR
VarDomFDR
RefDomFDR
""".split() if _f]

if dogl:
    outheaders.extend(glheaders)

dodebug = True
if opt.extended == None or 'Read filtering statistics' not in opt.extended:
    dodebug = False

debugging = [_f for _f in """
OtherCountForward
OtherCountReverse
OtherCount
FilteredSNVLociReads
SNVLociReads
""".split() if _f]
debugging.extend(sorted(filter(lambda h: h != "BadRead",BadRead.allheaders)))

if dodebug:
    outheaders.extend(debugging)

pos = outheaders.index("SNVCountForward")
outheaders.insert(pos, 'ReadGroup')

outheaders1 = copy.copy(outheaders)

emptysym = None
outrows = []

if opt.snvbatchsize == None:
    snvblocksize = len(snvdata)
else:
    snvblocksize = opt.snvbatchsize

totalsnvs = 0
progress.stage("Count reads per SNV", len(snvdata))
for snvblock in range(0,len(snvdata),snvblocksize):

  if opt.threads == 0:
    pileups = SerialPileups(snvdata[snvblock:(snvblock+snvblocksize)],
                            opt.alignments, readfilter, chrreg, readgroup, umigroup).iterator()
  elif opt.threads > 0:
    pileups = MultiprocPileups(snvdata[snvblock:(snvblock+snvblocksize)], 
                               opt.alignments, readfilter, chrreg, readgroup, umigroup, procs=opt.threads).iterator()
  elif opt.threads < 0:
    pileups = ThreadedPileups(snvdata[snvblock:(snvblock+snvblocksize)], 
                              opt.alignments, readfilter, chrreg, readgroup, umigroup, threads=-opt.threads).iterator()

  for snvchr, snvpos, ref, alt, snvextra in snvdata[snvblock:(snvblock+snvblocksize)]:
    
    snvchr1, snvpos1, ref1, alt1, total, reads, badread = next(pileups)
    assert(snvchr == snvchr1 and snvpos == snvpos1)
    
    if opt.debug:
         print(snvchr,snvpos,ref,alt, \
             " ".join(map(str,[total[i] for i in range(len(opt.alignments))])), \
             " ".join(map(str,[badread[(i, 'Good')] for i in range(len(opt.alignments))])))

    allsi = set()
    goodreads = defaultdict(list)
    for al, pos, base, si, umi in reads:
        goodreads[base].append((si, al, umi))
        allsi.add(si)

    totalsnvs += 1

    altnucs = "".join(sorted(map(str.strip,alt.split(','))))
    othernucs = "".join(sorted(set('ACGT')-set([ref]+list(altnucs))))

    if not umigroup:
        counts = defaultdict(int)
        for base in goodreads:
            for si, al, umi in goodreads[base]:
                counts[(base, "R" if al.is_reverse else "F", si)] += 1
                counts[(base, si)] += 1
                counts[si] += 1
                if len(altnucs) > 1 and base in altnucs:
                    counts[(altnucs, "R" if al.is_reverse else "F", si)] += 1
                    counts[(altnucs, si)] += 1
                if len(othernucs) > 1 and base in othernucs:
                    counts[(othernucs, "R" if al.is_reverse else "F", si)] += 1
                    counts[(othernucs, si)] += 1
    else:
        umis = defaultdict(set)
        for base in goodreads:
            for si, al, umi in goodreads[base]:
                if umi == None:
                    continue
                umis[(base, "R" if al.is_reverse else "F", si)].add(umi)
                umis[(base, si)].add(umi)
                umis[si].add(umi)
                if len(altnucs) > 1 and base in altnucs:
                    umis[(altnucs, "R" if al.is_reverse else "F", si)].add(umi)
                    umis[(altnucs, si)].add(umi)
                if len(othernucs) > 1 and base in othernucs:
                    umis[(othernucs, "R" if al.is_reverse else "F", si)].add(umi)
                    umis[(othernucs, si)].add(umi)
        counts = defaultdict(int)
        for k in umis:
            counts[k] = len(umis[k])

    for si in sorted(allsi):
        nsnvf = counts[(altnucs,"F",si)]
        nsnvr = counts[(altnucs,"R",si)]
        nsnv = counts[(altnucs,si)]
        nreff = counts[(ref, "F", si)]
        nrefr = counts[(ref, "R", si)]
        nref = counts[(ref, si)]
        notherf = counts[(othernucs,"F",si)]
        notherr = counts[(othernucs,"R",si)]
        nother = counts[(othernucs,si)]
        counted = counts[si]

        if isinstance(si,int):
            alf = opt.alignments[si]
            rg = os.path.split(alf)[1].rsplit('.', 1)[0]
        elif len(opt.alignments) == 1:
            rg = si[1]
        else:
            alf = opt.alignments[si[0]]
            rg = os.path.split(alf)[1].rsplit('.', 1)[0] + ":"+si[1]

        row = [ snvchr, snvpos, ref, alt, rg ] + \
              [nsnvf, nsnvr,
               nreff, nrefr,
               nsnv, nref,
               counted,
               100.0 * (total[si] - badread[si, 'Good']) /
                           float(total[si]) if total[si] != 0 else 0.0,
               float(nsnv)/(nsnv+nref) if (nsnv+nref)>0 else None]  + \
               ( [] if not dogl else [ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ] ) + \
               ( [] if not dodebug else [ notherf, notherr, nother, badread[si, 'Good'], total[si]] )

        if dodebug:
          for s in sorted(BadRead.allheaders):
            if s == 'BadRead':
                continue
            row.append(badread[si, s])

        outrows.append(row)

    progress.update()
progress.done()
if not opt.quiet:
    print("SNVs/sec: %.2f"%(float(totalsnvs)/progress.duration(),))

# Determine the maxreads value, if percentile, otherwise let the defaultdict take care of it
coverage = defaultdict(list)
maxreads = defaultdict(lambda: int(opt.maxreads))
if 0 < opt.maxreads < 1:
    for r in [dict(list(zip(outheaders,r))) for r in outrows]:
        coverage[r['ReadGroup']].append(r['GoodReads'])
    for al in coverage:
        n = len(coverage[al])
        percind = int(round(n*opt.maxreads))
        maxreads[al] = sorted(coverage[al])[percind]

if dogl:

  for i in range(len(outrows)):

    #Exctract the counts and rescale if necessary
    r = dict(list(zip(outheaders, outrows[i])))
    al,nsnv,nref,nother,counted = list(map(r.get,["ReadGroup","SNVCount","RefCount","OtherCount","GoodReads"]))

    if counted > maxreads[al]:
        factor = float(maxreads[al])/float(counted)
        nsnv = int(round(factor*nsnv))
        nref = int(round(factor*nref))
        nother = int(round(factor*nother))

    #Compute p-values
    pcount = 0.5
    n = nsnv + nref + nother
    nprime = nsnv + nref + nother + 4 * pcount
    q = float(nother + 2 * pcount) / (2 * nprime)
    nothomoref = binom_test_high(nsnv, n, q)
    nothomovar = binom_test_high(nref, n, q)
    if nsnv > nref:
        nothet = binom_test_high(nsnv, nsnv + nref, 0.5)
        refdom = 1.0
        vardom = nothet
    elif nref > nsnv:
        nothet = binom_test_high(nref, nsnv + nref, 0.5)
        vardom = 1.0
        refdom = nothet
    else:
        nothet = 1.0
        vardom = 1.0
        refdom = 1.0

    # And store in the output rows...
    pos = outheaders.index("NotHomoRefpV")
    outrows[i][pos] = nothomoref
    pos = outheaders.index("NotHetpV")
    outrows[i][pos] = nothet
    pos = outheaders.index("NotHomoVarpV")
    outrows[i][pos] = nothomovar
    pos = outheaders.index("VarDompV")
    outrows[i][pos] = vardom
    pos = outheaders.index("RefDompV")
    outrows[i][pos] = refdom

  # Now compute FDR and scores...

  pvkeys = [h for h in outheaders if h.endswith('pV')]
  fdrkeys = [h for h in outheaders if h.endswith('FDR')]
  allpvals = []
  n = len(outrows)
  for pvk in pvkeys:
      pos = outheaders.index(pvk)
      allpvals.extend(list(map(itemgetter(pos), outrows)))
  # print allpvals
  allfdrs = fdr(allpvals)

  for j, fdrk in enumerate(fdrkeys):
    pos1 = outheaders.index(fdrk)
    for i in range(len(outrows)):
        outrows[i][pos1] = allfdrs[(j * n) + i]

  for i in range(len(outrows)):
    r = dict(list(zip(outheaders, outrows[i])))
    homovarsc = max(0.0, min(pvscore(r["NotHetFDR"]),
                             pvscore(r["NotHomoRefFDR"])) - pvscore(r["NotHomoVarFDR"]))
    homorefsc = max(0.0, min(pvscore(r["NotHetFDR"]),
                             pvscore(r["NotHomoVarFDR"])) - pvscore(r["NotHomoRefFDR"]))
    hetsc = max(0.0, min(pvscore(r["NotHomoRefFDR"]),
                         pvscore(r["NotHomoVarFDR"])) - pvscore(r["NotHetFDR"]))
    vardomsc = pvscore(r["VarDomFDR"])
    refdomsc = pvscore(r["RefDomFDR"])
    pos = outheaders.index("HomoVarSc")
    outrows[i][pos] = homovarsc
    pos = outheaders.index("HomoRefSc")
    outrows[i][pos] = homorefsc
    pos = outheaders.index("HetSc")
    outrows[i][pos] = hetsc
    pos = outheaders.index("VarDomSc")
    outrows[i][pos] = vardomsc
    pos = outheaders.index("RefDomSc")
    outrows[i][pos] = refdomsc

if opt.output:
    filename = opt.output
    base, extn = filename.rsplit('.', 1)
    extn = extn.lower()
    if extn == 'csv':
        output = CSVFileTable(filename=filename, headers=outheaders1)
    elif extn == 'tsv':
        output = TSVFileTable(filename=filename, headers=outheaders1)
    elif extn == 'xls':
        output = XLSFileTable(
            filename=filename, headers=outheaders1, sheet='Results')
    elif extn == 'xlsx':
        output = XLSXFileTable(
            filename=filename, headers=outheaders1, sheet='Results')
    elif extn == 'txt':
        output = TXTFileTable(filename=filename, headers=outheaders1)
        emptysym = "-"
    else:
        raise RuntimeError("Unexpected output file extension: %s" % filename)
else:
    output = TXTFileTable(filename=sys.stdout, headers=outheaders1)
    emptysym = "-"

progress.stage('Output results')
if opt.output:
    outdir = os.path.split(opt.output)[0]
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir,exist_ok=True)
output.from_rows(
    map(lambda r: dict(zip(outheaders, list(map(lambda v: emptysym if (v == None) else v,r)) + [emptysym] * 50)),outrows))
progress.done()

