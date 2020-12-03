#!/bin/env python2.7
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
from util import *
from fisher import *
from pileups import SerialPileups, ThreadedPileups
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

advanced = OptionGroup(parser, "Advanced")
parser.add_option("-s", "--snvs", type="files", dest="snvs", default=None,
                  help="Single-Nucleotide-Variant files. Required.", name="SNV Files",
                  notNone=True, remember=True,
                  filetypes=[("SNV Files", "*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignment files in indexed BAM format. Required.", name="Read Alignment Files",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignment Files (indexed BAM)", "*.bam")])
advanced.add_option("-m", "--minreads", type="int", dest="minreads", default=10, remember=True,
                    help="Minimum number of good reads at SNV locus per alignment file. Default=10.", name="Min. Reads")
advanced.add_option("-M", "--minphased", type="int", dest="minphased", default=10, remember=True,
                    help="Minimum number of phased reads at SNV locus summed over all alignment files. Default=10.", name="Min. Phased")
advanced.add_option("-F", "--full", action="store_true", dest="full", default=False, remember=True,
                    help="Output extra diagnostic read count fields. Default=False.", name="All Fields")
# advanced.add_option("-f", "--alignmentfilter", action="store_false", dest="filter", default=True, remember=True,
#                     help="(Turn off) alignment filtering by length, edits, etc.", name="Filter Alignments")
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

opt.minreads = max(opt.minreads,1)

opt.minphased = max(opt.minphased,0)

progress = None
if not opt.output:
    opt.quiet = True
progress = ProgressText(quiet=opt.quiet)

from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable, BEDFile, VCFFile

progress.stage("Read SNV data", len(opt.snvs))
snvheaders = [_f for _f in """
CHROM POS REF ALT
""".split() if _f]

snvdata = {}
extrasnvheaders = []
usedsnvheaders = set()
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
        if h not in extrasnvheaders:
            extrasnvheaders.append(h)

    for r in snvs:
        chr = r[snvheaders[0]].strip()
        snvchroms[filename].add(chr)
        locus = int(r[snvheaders[1]].strip())
        ref = r[snvheaders[2]].strip()
        alt = r[snvheaders[3]].strip()
        # if r.get('varType') != 'SNP':
        #     continue
        if r.get('INFO:INDEL'):
            continue
        if len(ref) != 1:
            continue
        if not re.search(r'^[ACGT](,[ACGT])*$', alt):
            continue
        for h in r:
            if r.get(h):
                usedsnvheaders.add(h)
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

for bamfile in opt.alignments:
    chrreg.add_bamlabels(bamfile)

chrreg.determine_chrom_order()

snvdata = sorted(list(snvdata1.values()),key=lambda s: (chrreg.chrom_order(s[0]),s[1],s[2],s[3]))
extrasnvheaders = [h for h in extrasnvheaders if h in usedsnvheaders]
progress.message("SNVs: %d\n" % len(snvdata))

outheaders = snvheaders + [_f for _f in """
SNVCountPhase0
SNVCountPhase1
SNVCountPhase2
RefCountPhase0
RefCountPhase1
RefCountPhase2
SNVCountPhased
RefCountPhased
SNVCount
RefCount
PhasedCount
GoodReads
%BadRead
""".split() if _f]

debugging = [_f for _f in """
OtherCount
RemovedDuplicateReads
FilteredSNVLociReads
SNVLociReads
""".split() if _f]
debugging.extend(sorted(BadRead.allheaders))

outheaders.extend(debugging)

pos = outheaders.index("SNVCountPhase0")
# outheaders.insert(pos, 'AlignedReads')
for h in reversed(extrasnvheaders):
   outheaders.insert(pos,h)

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
#     snvdata = sorted(random.sample(snvdata,10000))
#     snvdata = sorted(sorted(random.sample(snvdata,200))*5)
#     snvdata = sorted(random.sample(snvdata,200))*5

# if opt.filter:
#     readfilter = SNVPileupReadFilter()
# else:
readfilter = BasicFilter()

if opt.tpb == 0:
    pileups = SerialPileups(snvdata, opt.alignments, readfilter, chrreg).iterator()
else:
    pileups = ThreadedPileups(snvdata, opt.alignments, readfilter, chrreg, threadsperbam=opt.tpb).iterator()

progress.stage("Count reads per SNV", len(snvdata))

totalsnvs = 0
start = time.time()
for snvchr, snvpos, ref, alt, snvextra in snvdata:
    
    snvchr1, snvpos1, ref1, alt1, total, reads, badread = next(pileups)
    assert(snvchr == snvchr1 and snvpos == snvpos1)
    
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

    totalsnvs += 1

    counts = defaultdict(int)
    counts1 = defaultdict(int)
    for base in goodreads:
        for si, al in goodreads[base]:
            try:
                hp = al.opt('HP')
            except KeyError:
                hp = 0
            counts[(base, "R" if al.is_reverse else "F", si)] += 1
            counts1[(base, hp, si)] += 1
    mincounted = 1e+20
    for si, alf in enumerate(opt.alignments):
        counted = sum([counts[(t[0], t[1], si)] for t in [
                      (n, d) for n in 'ACGT' for d in 'FR']])
        mincounted = min(counted, mincounted)
    if mincounted < opt.minreads:
        continue
    phased = 0
    for si, alf in enumerate(opt.alignments):
        phased += sum([counts1[(t[0], t[1], si)] for t in [
                          (n, hp) for n in 'ACGT' for ph in (1,2)]])
    if phased < opt.minphased:
        continue    
    for si, alf in enumerate(opt.alignments):
        nsnvf = sum([counts[(nuc, "F", si)] for nuc in list(map(str.strip,alt.split(',')))])
        nsnvr = sum([counts[(nuc, "R", si)] for nuc in list(map(str.strip,alt.split(',')))])
        nsnv = nsnvr + nsnvf
        nreff = counts[(ref, "F", si)]
        nrefr = counts[(ref, "R", si)]
        nref = nreff + nrefr
        othernucs = set('ACGT') - set([ref] + alt.split(','))
        notherf = sum([counts[(nuc, "F", si)] for nuc in othernucs])
        notherr = sum([counts[(nuc, "R", si)] for nuc in othernucs])
        nother = notherf + notherr
        counted = sum([counts[(t[0], t[1], si)] for t in [
                      (n, d) for n in 'ACGT' for d in 'FR']])
        nrefph0 = counts1[(ref, 0, si)]
        nrefph1 = counts1[(ref, 1, si)]
        nrefph2 = counts1[(ref, 2, si)]
        nrefphased = nrefph1 + nrefph2
        nsnvph0 = sum([counts1[(nuc, 0, si)] for nuc in list(map(str.strip,alt.split(',')))])
        nsnvph1 = sum([counts1[(nuc, 1, si)] for nuc in list(map(str.strip,alt.split(',')))])
        nsnvph2 = sum([counts1[(nuc, 2, si)] for nuc in list(map(str.strip,alt.split(',')))])
        nsnvphased = nsnvph1 + nsnvph2
        nphased = nrefphased + nsnvphased
        # [ os.path.split(alf)[1].rsplit('.', 1)[0] ] + \
        row = [ snvchr, snvpos, ref, alt ] + \
              list(map(snvextra.get,extrasnvheaders)) + \
              [nsnvph0, nsnvph1, nsnvph2,
               nrefph0, nrefph1, nrefph2,
               nsnvphased,nrefphased,
               nsnv, nref,
               nphased,
               counted,
               100.0 * (total[si] - badread[si, 'Good']) /
               float(total[si]) if total[si] != 0 else 0.0,
               nother,
               duplicates_removed[si],
               badread[si, 'Good'],
               total[si]]

        for s in sorted(BadRead.allheaders):
            row.append(badread[si, s])
        outrows.append(row)

    progress.update()
progress.done()
if not opt.quiet:
    print("SNVs/sec: %.2f"%(float(totalsnvs)/(time.time()-start),))

progress.stage('Output results')
output.from_rows(
    [dict(list(zip(outheaders, r + [emptysym] * 50))) for r in outrows])
progress.done()
