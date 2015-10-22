#!/bin/env python27
import sys, os, os.path, glob, copy, traceback, time, re, csv, tempfile, urllib, shutil, atexit, subprocess, time, math
from collections import defaultdict
from os.path import join, dirname, realpath
try:
    sys.path.append(join(dirname(realpath(__file__)),'..','..','common','src'))
except NameError:
    pass
from optparse_gui import OptionParser, OptionGroup, GUI, UserCancelledError, ProgressText
from util import *

from version import VERSION
VERSION='1.0.0 (%s)'%(VERSION,)

def excepthook(etype,value,tb):
    traceback.print_exception(etype,value,tb)
    print >>sys.stderr, "Type <Enter> to Exit...",
    sys.stderr.flush()
    raw_input()
sys.excepthook = excepthook

toremove = []
def cleanup():
    for d in toremove:
    	shutil.rmtree(d,ignore_errors=True)
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
parser.add_option("-s","--snps",type="files",dest="snps",default=None,
	          help="Single-Nucleotide-Polymophisms. Required.", name="SNPs",
                  notNone=True, remember=True,
		  filetypes=[("SNPs","*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("-r","--readalignments",type="files",dest="alignments",default=None,
	          help="Read alignments in BAM/SAM format. Required.", name="Read Alignments",
                  notNone=True, remember=True,
		  filetypes=[("Read Alignments (BAM/SAM Format)","*.bam;*.sam")])
advanced.add_option("-M","--mincount",type="int",dest="mincount",default=3,remember=True,
                  help="Minimum number of reads for reference and variant allelels to apply LoH test. Default: 3.",
                  name="Min. Count")
advanced.add_option("-F","--full",action="store_true",dest="full",default=False,remember=True,
	          help="Output extra diagnostic read count fields. Default=False.",name="All Fields")
advanced.add_option("-U","--uniquereads",action="store_true",dest="unique",default=False,remember=True,
                  help="Consider only distinct reads.",name="Unique Reads")
advanced.add_option("-q","--quiet",action="store_true",dest="quiet",default=False,remember=True,
                  help="Quiet.",name="Quiet")
parser.add_option("-o","--output",type="savefile",dest="output",remember=True,
                  help="Output file. Leave empty for console ouptut.",default="",
                  name="Output File", filetypes=[("All output formats","*.xlsx;*.xls;*.csv;*.tsv;*.txt"),
                                                 ("Excel","*.xlsx"),("Excel2003","*.xls"),
                                                 ("CSV","*.csv"),("TSV","*.tsv"),("Text","*.txt")])
parser.add_option_group(advanced)

opt = None
while True:
    if 'exit' in error_kwargs:
        try:
            opt,args = parser.parse_args(opts=opt)
        except UserCancelledError:
            sys.exit(0);
    else:
        opt,args = parser.parse_args()
        
    break

progress = None
if not opt.output:
    opt.quiet = True
progress = ProgressText(quiet=opt.quiet)

import pysam
from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable, BEDFile, VCFFile

progress.stage("Read SNP data",len(opt.snps))
snpheaders = filter(None,"""
CHROM POS REF ALT
""".split())

snpdata = {}
extrasnpheaders = []
usedsnpheaders = set()
for filename in opt.snps:
    
    base,extn = filename.rsplit('.',1)
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
        snps = TXTFileTable(filename=filename,headers=snpheaders)
    else:
        raise RuntimeError("Unexpected SNP file extension: %s"%filename)

    for h in snpheaders:
        if h not in snps.headers():
            raise RuntimeError("Required header: %s missing from SNP file %s"%(h,filename))

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
	if not re.search(r'^[ACGT](,[ACGT])*$',alt):
	    continue
	for h in r:
	    if r.get(h):
		usedsnpheaders.add(h)
	cannonr = (",".join(map(lambda t: "%s:%s"%t,sorted(r.items()))))
        snpkey = (chr,locus,ref,alt,cannonr)
        if snpkey not in snpdata:
            snpdata[snpkey] = (chr,locus,ref,alt,r)

    progress.update()
progress.done()
snpdata = sorted(snpdata.values())
extrasnpheaders = filter(lambda h: h in usedsnpheaders,extrasnpheaders)
progress.message("SNPs: %d"%len(snpdata))

samfiles=[]
for al in opt.alignments:
  if al.lower().endswith('.bam'):
    samfile = pysam.Samfile(al,"rb")
  elif al.lower().endswith('.sam'):
    samfile = pysam.Samfile(al,"r")
  else:
    raise RuntimeError("Unexpected alignments file extension: %s."%al)
  samfiles.append(samfile)

outheaders = snpheaders + filter(None,"""
SNPCount
NoSNPCount
Prob
LogOdds
P-Value
Bonferroni
FDR
""".split())

debugging = filter(None,"""
OtherCount
GoodReads
RemovedDuplicateReads
FilteredSNPLociReads
SNPLociReads
""".split())

outheaders.extend(debugging)

pos = outheaders.index("SNPCount")
for h in reversed(extrasnpheaders):
    outheaders.insert(pos,h)
    
outheaders1 = copy.copy(outheaders)
if not opt.full:
    for dh in debugging:
	if dh in outheaders1:
	    outheaders1.remove(dh)

emptysym = None
if opt.output:
    filename = opt.output
    base,extn = filename.rsplit('.',1)
    extn = extn.lower()
    if extn == 'csv':
        output = CSVFileTable(filename=filename,headers=outheaders1)
    elif extn == 'tsv':
        output = TSVFileTable(filename=filename,headers=outheaders1)
    elif extn == 'xls':
        output = XLSFileTable(filename=filename,headers=outheaders1,sheet='Results')
    elif extn == 'xlsx':
        output = XLSXFileTable(filename=filename,headers=outheaders1,sheet='Results')         
    elif extn == 'txt':
        output = TXTFileTable(filename=filename,headers=outheaders1)
    else:
        raise RuntimeError("Unexpected output file extension: %s"%filename)
else:
    output = TXTFileTable(filename=sys.stdout,headers=outheaders1)
    emptysym = "-"

outrows = []

from fisher import fisher_exact, bonferroni, fdr, lod, binom_test
pvalues = []

progress.stage("Count reads per SNP",len(snpdata))

filter = SNPPileupReadFilter()

for snpchr,snppos,ref,alt,snpextra in snpdata:

    reads=[]; total = 0
    snppos1 = snppos - 1
    for i,samfile in enumerate(samfiles):
        for pileupcolumn in samfile.pileup( snpchr, snppos1, snppos1+1, truncate=True):
            total += pileupcolumn.n
            for pileupread in pileupcolumn.pileups:
                try:
                    al,pos,base = filter.test(pileupread)
                except BadRead:
                    continue
                reads.append((al,pos,base,i))

    goodreads = defaultdict(list)
    for al,pos,base,si in reads:
        goodreads[base].append((si,al))

    # Deduplicate the reads based on the read sequence or the
    # start and end of the alignment or ???
    duplicates_removed = 0
    if opt.unique:
        for base in goodreads:
            seen = set(); retain = list()
            for si,al in goodreads[base]:
                if (si,al.seq) not in seen:
                    retain.append((si,al))
                    seen.add((si,al.seq))
                else:
                    duplicates_removed += 1
            goodreads[base] = retain

    # goodreads now contains the relevant read alignments.

    if len(goodreads[ref]) < opt.mincount or len(goodreads[alt]) < opt.mincount:
        progress.update()
        continue

    counts = defaultdict(int)
    for base in goodreads:
        for si,al in goodreads[base]:
            counts[base] += 1

    nsnp = sum(map(lambda nuc: counts[nuc], alt.split(',')))
    nref = counts[ref]
    nother = sum(counts.values())-nsnp-nref
    p = emptysym; pval = emptysym; logprob = emptysym;

    psnp =  nsnp/float(nsnp+nref)
    pref =  nref/float(nsnp+nref)
    p = psnp/(psnp + pref)
    logodds = math.log(float(nsnp)/float(nref),2.0)
    pval = binom_test(nsnp,nsnp+nref,0.5)
    pvalues.append(pval)

    row = [ snpchr,snppos,ref,alt ] + \
          [ snpextra.get(k,emptysym) for k in extrasnpheaders ] + \
          [ nsnp,
            nref,
            p,
            logodds,
            pval,
            nother,
            sum(counts.values()),
            duplicates_removed,
            len(reads),
	    total ]
    outrows.append(row)
    
    progress.update()
progress.done()

progress.stage('Multiple-test correction and FDR computation')
bonf = bonferroni(pvalues)
fdr =  fdr(pvalues)

i = 0
pvalpos = outheaders.index('P-Value')
for r in outrows:
    if len(r) > pvalpos:
        if r[pvalpos] != emptysym:
            r.insert(pvalpos+1,fdr[i])
            r.insert(pvalpos+1,bonf[i])
            i += 1
        else:
            r.insert(pvalpos+1,emptysym)
            r.insert(pvalpos+1,emptysym)
progress.done()

progress.stage('Output results')
output.from_rows(map(lambda r: dict(zip(outheaders,r+[emptysym]*50)),outrows))
progress.done()

