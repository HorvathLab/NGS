#!/bin/env python27
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
import urllib
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
from util import *
from fisher import *
from operator import itemgetter

from version import VERSION
VERSION = '1.0.3 (%s)' % (VERSION,)

def excepthook(etype, value, tb):
    traceback.print_exception(etype, value, tb)
    print >>sys.stderr, "Type <Enter> to Exit...",
    sys.stderr.flush()
    raw_input()

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

advanced = OptionGroup(parser, "Advanced")
parser.add_option("-s", "--snps", type="files", dest="snps", default=None,
                  help="Single-Nucleotide-Polymophisms. Required.", name="SNPs",
                  notNone=True, remember=True,
                  filetypes=[("SNPs", "*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignments in BAM/SAM format. Required.", name="Read Alignments",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignments (BAM/SAM Format)", "*.bam;*.sam")])
advanced.add_option("-m", "--minreads", type="int", dest="minreads", default=10, remember=True,
                    help="Minimum number of good reads at SNP locus per alignment file. Default=10.", name="Min. Reads")
advanced.add_option("-F", "--full", action="store_true", dest="full", default=False, remember=True,
                    help="Output extra diagnostic read count fields. Default=False.", name="All Fields")
advanced.add_option("-f", "--alignmentfilter", action="store_false", dest="filter", default=True, remember=True,
                    help="(Turn off) alignment filtering by length, edits, etc.", name="Filter Alignments")
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

import pysam
from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable, BEDFile, VCFFile

progress.stage("Read SNP data", len(opt.snps))
snpheaders = filter(None, """
CHROM POS REF ALT
""".split())

snpdata = {}
extrasnpheaders = []
usedsnpheaders = set()
for filename in opt.snps:

    base, extn = filename.rsplit('.', 1)
    extn = extn.lower()
    if extn == 'csv':
        snps = CSVFileTable(filename=filename)
    elif extn == 'vcf':
        snps = VCFFile(filename=filename)
    elif extn == 'tsv':
        snps = TSVFileTable(filename=filename)
    elif extn == 'xls':
        snps = XLSFileTable(filename=filename)
    elif extn == 'xlsx':
        snps = XLSXFileTable(filename=filename)
    elif extn == 'txt':
        snps = TXTFileTable(filename=filename, headers=snpheaders)
    else:
        raise RuntimeError("Unexpected SNP file extension: %s" % filename)

    for h in snpheaders:
        if h not in snps.headers():
            raise RuntimeError(
                "Required header: %s missing from SNP file %s" % (h, filename))

    for h in snps.headers():
        if h in snpheaders:
            continue
        if h not in extrasnpheaders:
            extrasnpheaders.append(h)

    for r in snps:
        chr = r[snpheaders[0]]
        locus = int(r[snpheaders[1]])
        ref = r[snpheaders[2]]
        alt = r[snpheaders[3]]
        if r.get('INFO:INDEL'):
            continue
        if len(ref) != 1:
            continue
        if not re.search(r'^[ACGT](,[ACGT])*$', alt):
            continue
        for h in r:
            if r.get(h):
                usedsnpheaders.add(h)
        # cannonr = (",".join(map(lambda t: "%s:%s"%t,sorted(r.items()))))
        # snpkey = (chr,locus,ref,alt,cannonr)
        snpkey = (chr, locus, ref, alt)
        # if str(locus) not in "1337612851669120781870889032162692092190017123282449235487873132320332288913526053039530834572380747246127615660316838509772413593731219137348879478140329":
        #     continue
        if snpkey not in snpdata:
            snpdata[snpkey] = (chr, locus, ref, alt, r)

    progress.update()
progress.done()

snpdata = sorted(snpdata.values())
extrasnpheaders = filter(lambda h: h in usedsnpheaders, extrasnpheaders)
progress.message("SNPs: %d" % len(snpdata))

samfiles = []
for al in opt.alignments:
    if al.lower().endswith('.bam'):
        samfile = pysam.Samfile(al, "rb")
    elif al.lower().endswith('.sam'):
        samfile = pysam.Samfile(al, "r")
    else:
        raise RuntimeError("Unexpected alignments file extension: %s." % al)
    samfiles.append(samfile)

outheaders = snpheaders + filter(None, """
SNPCountForward
SNPCountReverse
RefCountForward
RefCountReverse
SNPCount
RefCount
GoodReads
%BadRead
HomoVarSc
HetSc
HomoRefSc
VarDomSc
RefDomSc
""".split())

debugging = filter(None, """
OtherCountForward
OtherCountReverse
OtherCount
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
RemovedDuplicateReads
FilteredSNPLociReads
SNPLociReads
""".split())
debugging.extend(sorted(BadRead.allheaders))

outheaders.extend(debugging)

pos = outheaders.index("SNPCountForward")
outheaders.insert(pos, 'AlignedReads')
# for h in reversed(extrasnpheaders):
#    outheaders.insert(pos,h)

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

progress.stage("Count reads per SNP", len(snpdata))

if opt.filter:
    readfilter = SNPPileupReadFilter()
else:
    readfilter = BasicFilter()
totalsnps = 0
for snpchr, snppos, ref, alt, snpextra in snpdata:
    reads = []
    total = Counter()
    badread = Counter()
    snppos1 = snppos - 1
    for i, samfile in enumerate(samfiles):
        for pileupcolumn in samfile.pileup(snpchr, snppos1, snppos1 + 1, truncate=True):
            total[i] += pileupcolumn.n
            for pileupread in pileupcolumn.pileups:
                try:
                    al, pos, base, nseg = readfilter.test(pileupread)
                except BadRead, e:
                    badread[(i, e.message)] += 1
                    continue
                reads.append((al, pos, base, i))
                badread[i, 'Good'] += 1

    goodreads = defaultdict(list)
    for al, pos, base, si in reads:
        goodreads[base].append((si, al))

    # Deduplicate the reads based on the read sequence or the
    # start and end of the alignment or ???
    duplicates_removed = Counter()
    if opt.unique:
        for base in goodreads:
            seen = set()
            retain = list()
            for si, al in goodreads[base]:
                if (si, al.seq) not in seen:
                    retain.append((si, al))
                    seen.add((si, al.seq))
                else:
                    duplicates_removed[si] += 1
            goodreads[base] = retain

    # goodreads now contains the relevant read alignments.

    totalsnps += 1

    counts = defaultdict(int)
    for base in goodreads:
        for si, al in goodreads[base]:
            counts[(base, "R" if al.is_reverse else "F", si)] += 1
    mincounted = 1e+20
    for si, alf in enumerate(opt.alignments):
        counted = sum(map(lambda t: counts[(t[0], t[1], si)], [
                      (n, d) for n in 'ACGT' for d in 'FR']))
        mincounted = min(counted, mincounted)
    if mincounted < opt.minreads:
        continue
																																				
    for si, alf in enumerate(opt.alignments):
        nsnpf = sum(map(lambda nuc: counts[(nuc, "F", si)], map(str.strip,alt.split(','))))
        nsnpr = sum(map(lambda nuc: counts[(nuc, "R", si)], map(str.strip,alt.split(','))))
        nsnp = nsnpr + nsnpf
        nreff = counts[(ref, "F", si)]
        nrefr = counts[(ref, "R", si)]
        nref = nreff + nrefr
        othernucs = set('ACGT') - set([ref] + alt.split(','))
        notherf = sum(map(lambda nuc: counts[(nuc, "F", si)], othernucs))
        notherr = sum(map(lambda nuc: counts[(nuc, "R", si)], othernucs))
        nother = notherf + notherr
        counted = sum(map(lambda t: counts[(t[0], t[1], si)], [
                      (n, d) for n in 'ACGT' for d in 'FR']))

        pcount = 0.5
        n = nsnp + nref + nother
        nprime = nsnp + nref + nother + 4 * pcount
        q = float(nother + 2 * pcount) / (2 * nprime)
        nothomoref = binom_test_high(nsnp, n, q)
        nothomovar = binom_test_high(nref, n, q)
        if nsnp > nref:
            nothet = binom_test_high(nsnp, nsnp + nref, 0.5)
            refdom = 1.0
            vardom = nothet
        elif nref > nsnp:
            nothet = binom_test_high(nref, nsnp + nref, 0.5)
            vardom = 1.0
            refdom = nothet
        else:
            nothet = 1.0
            vardom = 1.0
            refdom = 1.0

        # [ snpextra.get(k,emptysym) for k in extrasnpheaders ] + \
        row = [ snpchr, snppos, ref, alt ] + \
              [ os.path.split(alf)[1].rsplit('.', 1)[0] ] + \
              [nsnpf, nsnpr,
               nreff, nrefr,
               nsnp, nref,
               counted,
               100.0 * (total[si] - badread[si, 'Good']) /
               float(total[si]) if total[si] != 0 else 0.0,
               -1, -1, -1, -1, -1,
               notherf, notherr,
               nother,
               nothomovar, nothomoref, nothet, vardom, refdom,
               -1, -1, -1, -1, -1,
               duplicates_removed[si],
               badread[si, 'Good'],
               total[si]]

        for s in sorted(BadRead.allheaders):
            row.append(badread[si, s])
        outrows.append(row)

    progress.update()
progress.done()

pvkeys = filter(lambda h: h.endswith('pV'), outheaders)
fdrkeys = filter(lambda h: h.endswith('FDR'), outheaders)
allpvals = []
n = len(outrows)
for pvk in pvkeys:
    pos = outheaders.index(pvk)
    allpvals.extend(map(itemgetter(pos), outrows))
# print allpvals
allfdrs = fdr(allpvals)

for j, fdrk in enumerate(fdrkeys):
    pos1 = outheaders.index(fdrk)
    for i in range(len(outrows)):
        outrows[i][pos1] = allfdrs[(j * n) + i]

for i in range(len(outrows)):
    r = dict(zip(outheaders, outrows[i]))
    homovarsc = max(0.0, min(pvscore(r["NotHetFDR"]), pvscore(
        r["NotHomoRefFDR"])) - pvscore(r["NotHomoVarFDR"]))
    homorefsc = max(0.0, min(pvscore(r["NotHetFDR"]), pvscore(
        r["NotHomoVarFDR"])) - pvscore(r["NotHomoRefFDR"]))
    hetsc = max(0.0, min(pvscore(r["NotHomoRefFDR"]), pvscore(
        r["NotHomoVarFDR"])) - pvscore(r["NotHetFDR"]))
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

progress.stage('Output results')
output.from_rows(
    map(lambda r: dict(zip(outheaders, r + [emptysym] * 50)), outrows))
progress.done()
