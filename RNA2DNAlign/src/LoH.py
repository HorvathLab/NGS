#!/bin/env python2.7
import sys
import os
import os.path
import glob
import copy
import traceback
import time
import re
import csv
import tempfile
import urllib.request, urllib.parse, urllib.error
import shutil
import atexit
import subprocess
import time
import math
from collections import defaultdict
from os.path import join, dirname, realpath
try:
    sys.path.append(join(dirname(realpath(__file__)),
                         '..', '..', 'common', 'src'))
except NameError:
    pass
from optparse_gui import OptionParser, OptionGroup, GUI, UserCancelledError, ProgressText
from util import *

from version import VERSION
VERSION = '%s' % (VERSION,)


def excepthook(etype, value, tb):
    traceback.print_exception(etype, value, tb)
    print("Type <Enter> to Exit...", end=' ', file=sys.stderr)
    sys.stderr.flush()
    input()
sys.excepthook = excepthook

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
else:
    parser = OptionParser(version=VERSION)
    error_kwargs = {}

advanced = OptionGroup(parser, "Advanced")
parser.add_option("-s", "--snvs", type="files", dest="snvs", default=None,
                  help="Single-Nucleotide-Variant files. Required.", name="SNVs",
                  notNone=True, remember=True,
                  filetypes=[("SNVs", "*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignments in BAM/SAM format. Required.", name="Read Alignments",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignments (BAM/SAM Format)", "*.bam;*.sam")])
advanced.add_option("-M", "--mincount", type="int", dest="mincount", default=3, remember=True,
                    help="Minimum number of reads for reference and variant allelels to apply LoH test. Default: 3.",
                    name="Min. Count")
advanced.add_option("-F", "--full", action="store_true", dest="full", default=False, remember=True,
                    help="Output extra diagnostic read count fields. Default=False.", name="All Fields")
advanced.add_option("-U", "--uniquereads", action="store_true", dest="unique", default=False, remember=True,
                    help="Consider only distinct reads.", name="Unique Reads")
advanced.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, remember=True,
                    help="Quiet.", name="Quiet")
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

    break

progress = None
if not opt.output:
    opt.quiet = True
progress = ProgressText(quiet=opt.quiet)

from pysamimport import pysam                                                   
from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable, BEDFile, VCFFile

progress.stage("Read SNV data", len(opt.snvs))
snvheaders = [_f for _f in """
CHROM POS REF ALT
""".split() if _f]

snvdata = {}
extrasnvheaders = []
usedsnvheaders = set()
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
        if h not in extrasnvheaders:
            extrasnvheaders.append(h)

    for r in snvs:
        chr = r[snvheaders[0]]
        locus = int(r[snvheaders[1]])
        ref = r[snvheaders[2]]
        alt = r[snvheaders[3]]
        if r.get('INFO:INDEL'):
            continue
        if len(ref) != 1:
            continue
        if not re.search(r'^[ACGT](,[ACGT])*$', alt):
            continue
        for h in r:
            if r.get(h):
                usedsnvheaders.add(h)
        cannonr = (",".join(["%s:%s" % t for t in sorted(r.items())]))
        snvkey = (chr, locus, ref, alt, cannonr)
        if snvkey not in snvdata:
            snvdata[snvkey] = (chr, locus, ref, alt, r)

    progress.update()
progress.done()
snvdata = sorted(snvdata.values())
extrasnvheaders = [h for h in extrasnvheaders if h in usedsnvheaders]
progress.message("SNVs: %d" % len(snvdata))

samfiles = []
for al in opt.alignments:
    if al.lower().endswith('.bam'):
        samfile = pysam.Samfile(al, "rb")
    elif al.lower().endswith('.sam'):
        samfile = pysam.Samfile(al, "r")
    else:
        raise RuntimeError("Unexpected alignments file extension: %s." % al)
    samfiles.append(samfile)

outheaders = snvheaders + [_f for _f in """
SNVCount
NoSNVCount
Prob
LogOdds
P-Value
Bonferroni
FDR
""".split() if _f]

debugging = [_f for _f in """
OtherCount
GoodReads
RemovedDuplicateReads
FilteredSNVLociReads
SNVLociReads
""".split() if _f]

outheaders.extend(debugging)

pos = outheaders.index("SNVCount")
for h in reversed(extrasnvheaders):
    outheaders.insert(pos, h)

outheaders1 = copy.copy(outheaders)
if not opt.full:
    for dh in debugging:
        if dh in outheaders1:
            outheaders1.remove(dh)

emptysym = None
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
    else:
        raise RuntimeError("Unexpected output file extension: %s" % filename)
else:
    output = TXTFileTable(filename=sys.stdout, headers=outheaders1)
    emptysym = "-"

outrows = []

from fisher import fisher_exact, bonferroni, fdr, lod, binom_test
pvalues = []

progress.stage("Count reads per SNV", len(snvdata))

filter = SNVPileupReadFilter()

for snvchr, snvpos, ref, alt, snvextra in snvdata:

    reads = []
    total = 0
    snvpos1 = snvpos - 1
    for i, samfile in enumerate(samfiles):
        for pileupcolumn in samfile.pileup(snvchr, snvpos1, snvpos1 + 1, truncate=True):
            total += pileupcolumn.n
            for pileupread in pileupcolumn.pileups:
                try:
                    al, pos, base = filter.test(pileupread)
                except BadRead:
                    continue
                reads.append((al, pos, base, i))

    goodreads = defaultdict(list)
    for al, pos, base, si in reads:
        goodreads[base].append((si, al))

    # Deduplicate the reads based on the read sequence or the
    # start and end of the alignment or ???
    duplicates_removed = 0
    if opt.unique:
        for base in goodreads:
            seen = set()
            retain = list()
            for si, al in goodreads[base]:
                if (si, al.seq) not in seen:
                    retain.append((si, al))
                    seen.add((si, al.seq))
                else:
                    duplicates_removed += 1
            goodreads[base] = retain

    # goodreads now contains the relevant read alignments.

    if len(goodreads[ref]) < opt.mincount or len(goodreads[alt]) < opt.mincount:
        progress.update()
        continue

    counts = defaultdict(int)
    for base in goodreads:
        for si, al in goodreads[base]:
            counts[base] += 1

    nsnv = sum([counts[nuc] for nuc in alt.split(',')])
    nref = counts[ref]
    nother = sum(counts.values()) - nsnv - nref
    p = emptysym
    pval = emptysym
    logprob = emptysym

    psnv = nsnv / float(nsnv + nref)
    pref = nref / float(nsnv + nref)
    p = psnv / (psnv + pref)
    logodds = math.log(float(nsnv) / float(nref), 2.0)
    pval = binom_test(nsnv, nsnv + nref, 0.5)
    pvalues.append(pval)

    row = [ snvchr, snvpos, ref, alt ] + \
          [ snvextra.get(k, emptysym) for k in extrasnvheaders ] + \
          [nsnv,
           nref,
           p,
           logodds,
           pval,
           nother,
           sum(counts.values()),
           duplicates_removed,
           len(reads),
           total]
    outrows.append(row)

    progress.update()
progress.done()

progress.stage('Multiple-test correction and FDR computation')
bonf = bonferroni(pvalues)
fdr = fdr(pvalues)

i = 0
pvalpos = outheaders.index('P-Value')
for r in outrows:
    if len(r) > pvalpos:
        if r[pvalpos] != emptysym:
            r.insert(pvalpos + 1, fdr[i])
            r.insert(pvalpos + 1, bonf[i])
            i += 1
        else:
            r.insert(pvalpos + 1, emptysym)
            r.insert(pvalpos + 1, emptysym)
progress.done()

progress.stage('Output results')
output.from_rows(
    [dict(list(zip(outheaders, r + [emptysym] * 50))) for r in outrows])
progress.done()
