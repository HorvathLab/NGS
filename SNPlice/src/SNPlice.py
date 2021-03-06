#!/bin/env python3
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
from collections import defaultdict, Counter
from os.path import join, dirname, realpath, split
try:
    scriptdir = dirname(realpath(__file__))
except NameError:
    scriptdir = dirname(realpath(sys.argv[0]))
scriptdir1 = realpath(join(scriptdir, '..', '..', 'common', 'src'))
sys.path.append(scriptdir1)
try:
    scriptextn = "." + os.path.split(sys.argv[0])[1].rsplit('.', 1)[1]
except:
    scriptextn = ""

from execute import Execute
execprog = Execute(scriptdir,scriptdir1,extn=scriptextn)

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
parser.add_option("-s", "--snps", type="files", dest="snps", default=None,
                  help="Single-Nucleotide-Polymophisms. Required.", name="SNPs",
                  notNone=True, remember=True,
                  filetypes=[("SNPs", "*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("-j", "--junctions", type="files", dest="junctions", default=None,
                  help="Splice junctions. Required.", name="Splice Junctions",
                  notNone=True, remember=True,
                  filetypes=[("Splice Junctions (BED Format)", "*.bed")])
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignments in BAM/SAM format. Required.", name="Read Alignments",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignments (BAM/SAM Format)", "*.bam;*.sam")])
advanced.add_option("-e", "--exoncoords", type="file", dest="exoncoords", default=None,
                  help="Exon coordinates for SNV filtering. Optional.", name="Exon Coords.",
                  remember=True,
                  filetypes=[("Exonic Coordinates", "*.txt")])
advanced.add_option("-d", "--distance", type="int", dest="dist", default=50, remember=True,
                    help="Upper bound on the distance between SNP locus and splice junction. Default: 50.",
                    name="Distance Bound")
advanced.add_option("-R", "--readthrough", type="int", dest="readthrough", default=5, remember=True,
                    help="Number of bases aligning into intron. Default: 5bp.",
                    name="Read-Through")
# advanced.add_option("-M", "--matepairs", action="store_true", dest="mates", default=False, remember=True,
#                     help="Consider the mate-pair reads for the detection of splicing. Default=False.", name="Mates")
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

opts.mates = False

progress = None
if not opt.output:
    opt.quiet = True
progress = ProgressText(quiet=opt.quiet)

import pysam
from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable, BEDFile, VCFFile

progress.stage("Read SNP data", len(opt.snps))
snpheaders = [_f for _f in """
CHROM POS REF ALT
""".split() if _f]

snvdata = {}
snvchroms = defaultdict(set)
extrasnpheaders = []
usedsnpheaders = set()
for filename in opt.snps:
    filename0 = filename
    base, extn = filename.rsplit('.', 1)
    extn = extn.lower()
    tempfilename = None
    if opt.exoncoords:
        if extn != 'vcf':
            extn = 'tsv'
        tempfileno,tempfilename = tempfile.mkstemp(suffix="."+extn)
        execprog.execute("exonicFilter", "--exons", opt.exoncoords,
                         "--input", filename, "--output", tempfilename)
        filename = tempfilename
    
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
        snvchroms[filename0].add(chr)
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
        # cannonr = (",".join(map(lambda t: "%s:%s" % t, sorted(r.items()))))
        snpkey = (filename0, chr, locus, ref, alt)
        if snpkey not in snpdata:
            snpdata[snpkey] = r

    if tempfilename != None:
        os.unlink(tempfilename)
    
    progress.update()
progress.done()

chrreg = ChromLabelRegistry()

for snvfile in snvchroms:
    chrreg.add_labels(snvfile,snvchroms[snvfile])

snpdata1 = {}
for (sf, chr, locus, ref, alt), r in snvdata.items():
    chrom = chrreg.label2chrom(sf,chr)
    assert(chrom)
    snvkey = (chrom,locus,ref,alt)
    if snvkey not in snvdata1:
        snvdata1[snvkey] = (chrom,locus,ref,alt,r)

for bamfile in opt.alignments:
    chrreg.add_bamlabels(bamfile)

progress.stage("Read splice junction data", len(opt.junctions))
juncdata = set()
juncchrom = defaultdict(set)
for juncfile in opt.junctions:
    junc = BEDFile(filename=juncfile)
    for r in junc:
        chr = r['chrom']
        juncchrom[juncfile].add(chr)
        st = int(r['chromStart'])
        ed = int(r['chromEnd'])
        bs = list(map(int, r['blockSizes'].split(',')))
        assert(len(bs) == 2)
        gap = (st + bs[0], ed - bs[1])
        key = (juncfile, chr, gap)
        juncdata.add(key)
    progress.update()
progress.done()

for juncfile in juncchroms:
    chrreg.add_labels(juncfile,juncchroms[juncfile])

juncdata1 = set()
for jf, chr, gap in juncdata:
    chrom = chrreg.label2chrom(jf,chr)
    juncdata1.add((chrom, gap))

chrreg.determine_chrom_order()

snpdata = sorted(list(snvdata1.values()),key=lambda s: (chrreg.chrom_order(s[0]),s[1],s[2],s[3]))

extrasnpheaders = [h for h in extrasnpheaders if h in usedsnpheaders]
progress.message("SNPs: %d" % len(snpdata))

juncdata = sorted(((chr, gap[i], gap) for i in (0, 1) for chr, gap in juncdata1), key=lambda j: (chrreg.chrom_order(j[0]),j[1],j[2]))
progress.message("Exon/Intron junctions: %d" % len(juncdata))

outheaders = snpheaders + [_f for _f in """
NumofJuncs
Distance
Junctions
SNPJuncIntronCount
SNPJuncNoIntronCount
NoSNPJuncIntronCount
NoSNPJuncNoIntronCount
Probability
LOD
P-Value
Bonferroni
FDR
%BadRead
""".split() if _f]

debugging = [_f for _f in """
SNPMateCount
NoSNPMateCount
SNPCount
NoSNPCount
MatesCount
NotMatesCount
IntronCount
NoIntronCount
SpanningReads
RemovedDuplicateReads
FilteredReads
SNPLociReads
""".split() if _f]
debugging.extend(sorted(BadRead.allheaders))

outheaders.extend(debugging)

pos = outheaders.index("NumofJuncs")
for h in reversed(extrasnpheaders):
    outheaders.insert(pos, h)

outheaders1 = copy.copy(outheaders)
if not opt.full:
    for dh in debugging:
        if dh in outheaders1:
            outheaders1.remove(dh)

if not opt.mates:
    for mh in "SNPMateCount NoSNPMateCount MatesCount NotMatesCount".split():
        if mh in outheaders1:
            outheaders1.remove(mh)

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

from fisher import fisher_exact, bonferroni, fdr, lod
pvalues = []
if opt.filter:
    alifilter = SNPPileupReadFilter(maxsegments=2)
    matefilter = ReadFilter(maxsegments=2)
else:
    alifilter = BasicFilter()
    matefilter = BasicFilter()

if opt.tpb == 0:
    pileups = SerialPileups(snpdata, opt.alignments, alifilter, chrreg, matedist=opt.dist).iterator()
else:
    pileups = MultiprocPileups(snpdata, opt.alignments, alifilter, chrreg, procsperbam=opt.tpb, matedist=opt.dist).iterator()

juncindst = 0
progress.stage("Count reads per SNP and splice junction", len(snpdata))
for snpchr, snppos, ref, alt, snpextra in snpdata:

    if juncindst < len(juncdata):
        juncchr, juncpos, (inst, ined) = juncdata[juncindst]
        while (juncchr < snpchr) or ((juncchr == snpchr) and (juncpos < (snppos - opt.dist))):
            juncindst += 1
            if juncindst >= len(juncdata):
                break
            juncchr, juncpos, (inst, ined) = juncdata[juncindst]

    juncinded = juncindst

    if juncinded < len(juncdata):
        juncchr, juncpos, (inst, ined) = juncdata[juncinded]
        while (juncchr < snpchr) or ((juncchr == snpchr) and (juncpos < (snppos + opt.dist))):
            juncinded += 1
            if juncinded >= len(juncdata):
                break
            juncchr, juncpos, (inst, ined) = juncdata[juncinded]

    njunc = (juncinded - juncindst)

    snvchr1, snvpos1, ref1, alt1, total, reads, badread = next(pileups)
    assert(snvchr == snvchr1 and snvpos == snvpos1)

    reads = []
    total = 0
    badread = Counter()
    if njunc > 0:
        snppos1 = snppos - 1
        for i, samfile in enumerate(samfiles):
            for pileupcolumn in samfile.pileup(snpchr, snppos1, snppos1 + 1, truncate=True):
                total += pileupcolumn.n
                for pileupread in pileupcolumn.pileups:
                    try:
                        al, pos, base, nseg = alifilter.test(pileupread)
                    except BadRead as e:
                        badread[e.message] += 1
                        continue
                    reads.append((al, pos, base, i, nseg, None))
                    badread['Good'] += 1

    if opt.mates:
        matedreads = []
        for r in reads:
            read, qpos, base, sind, nseg, dummy = r
            if read.is_paired and not read.mate_is_unmapped and nseg == 1 and \
                    read.tid == read.rnext and abs(read.pos - read.pnext) < opt.dist - read.qlen:
                mate = samfiles[sind].mate(read)
                try:
                    matefilter.test(mate)
                    matedreads.append((read, qpos, base, sind, nseg, mate))
                except BadRead as e:
                    matedreads.append((read, qpos, base, sind, nseg, None))
            else:
                matedreads.append((read, qpos, base, sind, nseg, None))
        reads = matedreads

    any = False
    for juncind in range(juncindst, juncinded):
        juncchr, juncpos, (inst, ined) = juncdata[juncind]

        # If intron is to the left of the SNP, we require the query
        # to read into (or over) the intron by some number of bases
        goodreads = defaultdict(list)
        if juncpos == ined:
            for al, pos, base, si, ns, mal in reads:
                if (ined + opt.readthrough) < al.aend and (al.aend - al.qlen) < (ined - opt.readthrough):
                    goodreads[base].append((si, al, False))
                elif mal and (ined + opt.readthrough) < mal.aend and (mal.aend - mal.qlen) < (ined - opt.readthrough):
                    goodreads[base].append((si, mal, True))
        elif juncpos == inst:
            for al, pos, base, si, ns, mal in reads:
                if al.pos < (inst - opt.readthrough) and (inst + opt.readthrough) < (al.pos + al.qlen):
                    goodreads[base].append((si, al, False))
                elif mal and mal.pos < (inst - opt.readthrough) and (inst + opt.readthrough) < (mal.pos + mal.qlen):
                    goodreads[base].append((si, mal, True))

        if len(goodreads) == 0:
            continue

        any = True

        # Deduplicate the reads based on the read sequence or the
        # start and end of the alignment or ???
        duplicates_removed = 0
        if opt.unique:
            for base in goodreads:
                seen = set()
                retain = list()
                for si, al, mate in goodreads[base]:
                    if (si, al.seq) not in seen:
                        retain.append((si, al, mate))
                        seen.add((si, al.seq))
                    else:
                        duplicates_removed += 1
                goodreads[base] = retain

        # goodreads now contains the relevant read alignments.
        # detect read-through of splice site by checking the
        # alignment span on the reference.
        # if len(goodreads[ref]) < 5 or len(goodreads[alt]) < 5:
        #     continue
        spliced = 'SPLICED'
        notspliced = 'NOTSPLICED'
        other = 'OTHER'
        all = 'ALL'
        mates = 'MATES'
        notmates = 'NOTMATES'
        counts = defaultdict(int)
        for base in goodreads:
            for si, al, mate in goodreads[base]:
                if abs(al.alen - al.qlen) < 2:
                    counts[(base, notspliced)] += 1
                elif abs(al.alen - al.qlen - (ined - inst)) < 2:
                    counts[(base, spliced)] += 1
                else:
                    counts[(base, other)] += 1
                if mate:
                    counts[(base, mates)] += 1
                else:
                    counts[(base, notmates)] += 1
            counts[(base, all)] = counts[(base, notspliced)] + \
                counts[(base, spliced)] + \
                counts[(base, other)]
            counts[notspliced] += counts[(base, notspliced)]
            counts[spliced] += counts[(base, spliced)]
            counts[other] += counts[(base, other)]
            counts[all] += counts[(base, all)]
            counts[mates] += counts[(base, mates)]
            counts[notmates] += counts[(base, notmates)]

        nsnpi = sum([counts[(nuc, notspliced)] for nuc in alt.split(',')])
        nsnpex = sum([counts[(nuc, spliced)] for nuc in alt.split(',')])
        nwti = counts[(ref, notspliced)]
        nwtex = counts[(ref, spliced)]
        p = emptysym
        pval = emptysym
        lodval = emptysym
        if (nsnpi + nsnpex) > 0 and \
           (nwti + nwtex) > 0 and \
           (nsnpi + nwti) > 0:

            psnp = nsnpi / float(nsnpi + nsnpex)
            pwt = nwti / float(nwti + nwtex)
            p = psnp / (psnp + pwt)

            pval = fisher_exact(nsnpi,
                                (nsnpi + nsnpex),
                                (nsnpi + nwti),
                                (nsnpi + nsnpex + nwti + nwtex))

            pvalues.append(pval)

            l = lod(nsnpi, (nsnpi + nsnpex), (nsnpi + nwti),
                    (nsnpi + nsnpex + nwti + nwtex))
            if l != None:
                lodval = l

        row = [ snpchr, snppos, ref, alt ] + \
              [ snpextra.get(k, emptysym) for k in extrasnpheaders ] + \
              [njunc,
               juncpos - snppos + 1 * (juncpos > snppos),
               "%s:%d-%d" % (juncchr, inst, ined),
               counts[(alt, notspliced)],
               counts[(alt, spliced)],
               counts[(ref, notspliced)],
               counts[(ref, spliced)],
               p,
               lodval,
               pval,
               100.0 * (total - len(reads)) / float(total),
               counts[(alt, mates)],
               counts[(ref, mates)],
               counts[(alt, all)],
               counts[(ref, all)],
               counts[mates],
               counts[notmates],
               counts[notspliced],
               counts[spliced],
               counts[all],
               duplicates_removed,
               len(reads),
               total]
        for s in sorted(BadRead.allheaders):
            row.append(badread[s])
        outrows.append(row)

    if not any:

        # ["%s:%d_%s/%s"%(snpchr,snppos,ref,alt)] \
        outrows.append([snpchr, snppos, ref, alt]
                       + [snpextra.get(k, emptysym) for k in extrasnpheaders]
                       + [0])

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
