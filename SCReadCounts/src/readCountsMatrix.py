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
    sys.path.append(join(dirname(realpath(__file__)),
                         '..', '..', 'ReadCounts', 'src'))
except NameError:
    pass
from optparse_gui import OptionParser, OptionGroup, GUI, UserCancelledError, ProgressText
from chromreg import ChromLabelRegistry

from version import VERSION
VERSION = '%s' % (VERSION,)

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

if not GUI() and len(sys.argv) == 2 and sys.argv[1] == '--GUI':
    from optparse_gui.needswx import *
    sys.exit(1)

if GUI() and len(sys.argv) == 1:
    from optparse_gui import OptionParserGUI
    parser = OptionParserGUI(version=VERSION)
    error_kwargs = {'exit': False}
    sys.excepthook = excepthook
else:
    parser = OptionParser(version=VERSION)
    error_kwargs = {}

minreads_default = 0

parser.add_option("-c", "--counts", type="files", dest="counts", default=None,
                  help="ReadCounts files. Required.", name="ReadCounts Files",
                  notNone=True, remember=True,
                  filetypes=[("ReadCounts Files", "*.csv;*.tsv;*.xls;*.xlsx;*.txt;")])
parser.add_option("-M", "--matrix", type="choice", dest="matrix", default="Ref:Var", remember=True,
                  help="Matrix output format: Ref:Var or VAF. Default: Ref:Var.", name="Matrix",
                  choices=["Ref:Var","VAF"])
parser.add_option("-m", "--minreads", type="int", dest="minreads", default=minreads_default, remember=True,
                  help="Minimum number of good reads at SNV locus. Default=no minimum.", name="Min. Reads")
parser.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, remember=True,
                  help="Quiet.", name="Quiet")
parser.add_option("-o", "--output", type="savefile", dest="output", remember=True,
                  help="Output file. Leave empty for console ouptut.", default="",
                  name="Output File", filetypes=[("All output formats", "*.xlsx;*.xls;*.csv;*.tsv;*.txt"),
                                                 ("Excel", "*.xlsx"), ("Excel2003", "*.xls"),
                                                 ("CSV", "*.csv"), ("TSV", "*.tsv"), ("Text", "*.txt")])
opt = None
while True:
    if 'exit' in error_kwargs:
        try:
            opt, args = parser.parse_args(opts=opt)
        except UserCancelledError:
            sys.exit(0)
    else:
        opt, args = parser.parse_args()

    break

matrix = None
if opt.matrix:
    if opt.matrix == "Ref:Var":
        matrix = (lambda d: "%(Ref)s:%(Var)s"%d)
    elif opt.matrix == "VAF":
        matrix = (lambda d: "%(VAF)s"%d)

progress = None
if not opt.output:
    opt.quiet = True
progress = ProgressText(quiet=opt.quiet)

doublequote = lambda s: '"%s"'%(s,)
indent = lambda s,n: "\n".join([(" "*n)+l for l in s.splitlines()])

args = []
args.extend(["-c",doublequote(" ".join(opt.counts))])
if matrix:
    args.extend(["-M",opt.matrix])
if opt.minreads != minreads_default:
    args.extend(["-m",str(opt.minreads)])
args.extend(["-o",doublequote(opt.output)])
if opt.quiet:
    args.extend(["-q"])

cmdargs = " ".join(args)

execution_log = """
readCounts Options:
  ReadCounts Files (-c): %s
  Matrix Output (-M):    %s
  Min. Reads (-m):       %s%s
  Quiet (-q):            %s
  Outfile File (-o):     %s

Command-Line: readCountsMatrix %s
"""%(", ".join(opt.counts),
     None if not matrix else opt.matrix,
     opt.minreads,
     "" if opt.matrix not in ("Ref:Var",) or opt.minreads == 0 else " (ignored)",
     opt.quiet,
     opt.output,
     cmdargs)

progress.message(execution_log)

from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable

progress.stage("Read ReadCounts input files", len(opt.counts))
headers = "CHROM POS REF ALT ReadGroup RefCount SNVCount GoodReads".split()

allrg = set()
vafmatrix = defaultdict(dict)
for filename in opt.counts:
    base, extn = filename.rsplit('.', 1)
    extn = extn.lower()
    if extn == 'csv':
        counts = CSVFileTable(filename=filename)
    elif extn == 'vcf':
        counts = VCFFile(filename=filename)
    elif extn == 'tsv':
        counts = TSVFileTable(filename=filename)
    elif extn == 'xls':
        counts = XLSFileTable(filename=filename)
    elif extn == 'xlsx':
        counts = XLSXFileTable(filename=filename)
    elif extn == 'txt':
        counts = TXTFileTable(filename=filename, headers=headers)
    else:
        raise RuntimeError("Unexpected ReadCounts file extension: %s" % filename)

    for h in headers:
        if h not in counts.headers():
            raise RuntimeError(
                "Required header: %s missing from ReadCounts file %s" % (h, filename))

    for r in counts:
        try:
             chr = int(r[headers[0]].strip())
        except ValueError:
             chr = r[headers[0]].strip()
        locus = int(r[headers[1]].strip())
        ref = r[headers[2]].strip()
        alt = r[headers[3]].strip()
        snvkey = (filename, chr, locus, ref, alt)
        rg=r[headers[4]].strip()
        nref = int(r[headers[5]].strip())
        nsnv = int(r[headers[6]].strip())
        nall = int(r[headers[7]].strip())
        allrg.add(rg)
        if opt.matrix in ("Ref:Var",) or nall >= opt.minreads:
            if (nref+nsnv) == 0:
                vaf = "NA"
            else:
                vaf = "%.6f"%(float(nsnv)/float(nref+nsnv),)
            vafmatrix[(chr,locus,ref,alt)][rg] = matrix(dict(Ref=nref,Var=nsnv,VAF=vaf))

    progress.update()
progress.done()

allrg = sorted(allrg)
chromsortkey = lambda c: (1e+20,c) if isinstance(c,str) else (c,"")
allloci = sorted(vafmatrix,key=chromsortkey)
outheaders = [ "SNV" ] + allrg
emptysym = None

if opt.output:
    filename = opt.output
    base, extn = filename.rsplit('.', 1)
    extn = extn.lower()
    if extn == 'csv':
        output = CSVFileTable(filename=filename, headers=outheaders)
    elif extn == 'tsv':
        output = TSVFileTable(filename=filename, headers=outheaders)
    elif extn == 'xls':
        output = XLSFileTable(
            filename=filename, headers=outheaders, sheet='Results')
    elif extn == 'xlsx':
        output = XLSXFileTable(
            filename=filename, headers=outheaders, sheet='Results')
    elif extn == 'txt':
        output = TXTFileTable(filename=filename, headers=outheaders)
    else:
        raise RuntimeError("Unexpected output file extension: %s" % filename)
else:
    output = TXTFileTable(filename=sys.stdout, headers=outheaders)
    emptysym = "-"

outrows = []
for key in allloci:
        snvkey = "%s:%s_%s>%s"%(key[0],key[1],key[2],key[3])
        row = [ snvkey ]
        for rg in allrg:
            row.append(vafmatrix[key].get(rg,matrix(dict(Ref=0,Var=0,VAF="NA"))))
        outrows.append(row)

progress.stage('Output results')
if opt.output:
    outdir = os.path.split(opt.output)[0]
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir,exist_ok=True)
output.from_rows(
    [dict(list(zip(outheaders, r + [emptysym] * 50))) for r in outrows])
progress.done()

