#!/bin/env python27
import sys
import os
import os.path
import glob
import copy
import traceback
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
from pileups import SerialPileups, ThreadedPileups
from chromreg import ChromLabelRegistry
from operator import itemgetter

from version import VERSION
VERSION = '1.0.4 (%s)' % (VERSION,)

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
                  help="Single-Nucleotide-Polymophism files. Required.", name="SNP Files",
                  notNone=True, remember=True,
                  filetypes=[("SNP Files", "*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignment files in indexed BAM format. Required.", name="Read Alignment Files",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignment Files (indexed BAM)", "*.bam")])
advanced.add_option("-m", "--minreads", type="int", dest="minreads", default=10, remember=True,
                    help="Minimum number of good reads at SNP locus per alignment file. Default=10.", name="Min. Reads")
advanced.add_option("-F", "--full", action="store_true", dest="full", default=False, remember=True,
                    help="Output extra diagnostic read count fields. Default=False.", name="All Fields")
advanced.add_option("-f", "--alignmentfilter", action="store_false", dest="filter", default=True, remember=True,
                    help="(Turn off) alignment filtering by length, edits, etc.", name="Filter Alignments")
advanced.add_option("-U", "--uniquereads", action="store_true", dest="unique", default=False, remember=True,
                    help="Consider only distinct reads.", name="Unique Reads")
advanced.add_option("-t", "--threadsperbam", type="int", dest="tpb", default=1, remember=True,
                    help="Worker threads per alignment file. Indicate no threading with 0. Default=1.", name="Threads/BAM")
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

    break

progress = None
if not opt.output:
    opt.quiet = True
progress = ProgressText(quiet=opt.quiet)

from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable, BEDFile, VCFFile

progress.stage("Read SNP data", len(opt.snps))
snpheaders = filter(None, """
CHROM POS REF ALT
""".split())

snpdata = {}
# extrasnpheaders = []
# usedsnpheaders = set()
snpchroms = defaultdict(set)
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
        # if h not in extrasnpheaders:
        #     extrasnpheaders.append(h)

    for r in snps:
        chr = r[snpheaders[0]].strip()
	snpchroms[filename].add(chr)
        locus = int(r[snpheaders[1]].strip())
        ref = r[snpheaders[2]].strip()
        alt = r[snpheaders[3]].strip()
        if r.get('INFO:INDEL'):
            continue
        if len(ref) != 1:
            continue
        if not re.search(r'^[ACGT](,[ACGT])*$', alt):
            continue
        # for h in r:
        #     if r.get(h):
        #         usedsnpheaders.add(h)
        snpkey = (filename, chr, locus, ref, alt)
        if snpkey not in snpdata:
            snpdata[snpkey] = r

    progress.update()
progress.done()

chrreg = ChromLabelRegistry()

for snpfile in snpchroms:
    chrreg.add_labels(snpfile,snpchroms[snpfile])

snpdata1 = {}
for (sf, chr, locus, ref, alt), r in snpdata.iteritems():
    chrom = chrreg.label2chrom(sf,chr)
    assert(chrom)
    snpkey = (chrom,locus,ref,alt)
    if snpkey not in snpdata1:
        snpdata1[snpkey] = (chrom,locus,ref,alt,r)

for bamfile in opt.alignments:
    chrreg.add_bamlabels(bamfile)

chrreg.determine_chrom_order()

snpdata = sorted(snpdata1.values(),key=lambda s: (chrreg.chrom_order(s[0]),s[1],s[2],s[3]))
# extrasnpheaders = filter(lambda h: h in usedsnpheaders, extrasnpheaders)
progress.message("SNPs: %d" % len(snpdata))

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

# if opt.debug:
#     import random
#     random.seed(1234567)
#     snpdata = sorted(random.sample(snpdata,10000))
#     snpdata = sorted(sorted(random.sample(snpdata,200))*5)
#     snpdata = sorted(random.sample(snpdata,200))*5

if opt.filter:
    readfilter = SNPPileupReadFilter()
else:
    readfilter = BasicFilter()

if opt.tpb == 0:
    pileups = SerialPileups(snpdata, opt.alignments, readfilter, chrreg).iterator()
else:
    pileups = ThreadedPileups(snpdata, opt.alignments, readfilter, chrreg, threadsperbam=opt.tpb).iterator()

progress.stage("Count reads per SNP", len(snpdata))

totalsnps = 0
start = time.time()
# for i in range(len(snpdata)):
for snpchr, snppos, ref, alt, snpextra in snpdata:
    
##     if opt.debug:
## 	if totalsnps % 100 == 0 and totalsnps > 0:
## 	    print "SNPs/sec: %.2f"%(float(totalsnps)/(time.time()-start),)

    snpchr1, snppos1, ref1, alt1, total, reads, badread = pileups.next()
    assert(snpchr == snpchr1 and snppos == snppos1)
    
##     if opt.debug:
##         print snpchr,snppos,ref,alt, \
##               " ".join(map(str,map(lambda i: total[i],range(len(opt.alignments))))), \
##               " ".join(map(str,map(lambda i: badread[(i, 'Good')],range(len(opt.alignments)))))

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
print "SNPs/sec: %.2f"%(float(totalsnps)/(time.time()-start),)

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
