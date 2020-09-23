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
from collections import defaultdict
sys.path.append(os.path.join(os.path.split(
    os.path.abspath(__file__))[0], '..', '..', 'common', 'src'))
from optparse_gui import OptionParser, OptionGroup, GUI, \
    UserCancelledError, ProgressText

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

if GUI() and len(sys.argv) == 1:
    from optparse_gui import OptionParserGUI
    parser = OptionParserGUI(version=VERSION)
    error_kwargs = {'exit': False}
else:
    parser = OptionParser(version=VERSION)
    error_kwargs = {}

parser.add_option("-c", "--counts", type="files", dest="counts", default=None,
                  help="Read counts per SNP/Junction", name="SNP/Junction Counts",
                  remember=True, notNone=True,
                  filetypes=[("Read counts", "*.xlsx;*.xls;*.csv;*.tsv;*.txt"),
                             ("Excel", "*.xlsx"), ("Excel2003", "*.xls"),
                             ("CSV", "*.csv"), ("TSV", "*.tsv")])
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

progress = None
if not opt.output:
    opt.quiet = True
progress = ProgressText(quiet=opt.quiet)

sumkeys = [_f for _f in map(str.strip, """
SNPJuncIntronCount SNPJuncNoIntronCount NoSNPJuncIntronCount NoSNPJuncNoIntronCount SNPMateCount NoSNPMateCount SNPCount NoSNPCount MatesCount NotMatesCount IntronCount NoIntronCount SpanningReads RemovedDuplicateReads SNPLociReads""".split()) if _f]
countdata = defaultdict(dict)
progress.stage("Read SNP/Junction counts")
from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable
countheaders = None
for filename in opt.counts:
    base, extn = filename.rsplit('.', 1)
    path, base = os.path.split(base)
    extn = extn.lower()
    if extn == 'csv':
        counts = CSVFileTable(filename=filename)
    elif extn == 'tsv':
        counts = TSVFileTable(filename=filename)
    elif extn == 'xls':
        counts = XLSFileTable(filename=filename)
    elif extn == 'xlsx':
        counts = XLSXFileTable(filename=filename)
    else:
        raise RuntimeError("Unexpected count file extension: %s" % filename)

    if countheaders == None:
        countheaders = counts.headers()
    else:
        assert countheaders == counts.headers()
    assert 'CHROM' in countheaders
    assert 'POS' in countheaders
    assert 'REF' in countheaders
    assert 'ALT' in countheaders
    assert 'Junctions' in countheaders

    for r in counts:
        for k in list(r.keys()):
            if r.get(k) in ("", None):
                del r[k]
        chr = r['CHROM']
        pos = r['POS']
        ref = r['REF']
        alt = r['ALT']
        # m = re.search(r'(.*):(\d+)_(.)/(.*)$',r['SNP'])
        # assert m
        try:
            chr = int(m.group(1))
        except ValueError:
            chr = m.group(1)
        # snp = (chr,int(m.group(2)),m.group(3),m.group(4))
        snp = (chr, int(pos), ref, alt)
        m1 = re.search(r'^(.*):(\d+)-(\d+)$', r.get('Junctions', ""))
        if m1:
            try:
                chr = int(m1.group(1))
            except:
                chr = m1.group(1)
            junc = (chr, int(m1.group(2)), int(m1.group(3)))
        else:
            junc = None
        key = (snp, junc)
        if key not in countdata:
            countdata[key] = dict(iter(r.items()))
            countdata[key]['Samples'] = base
            for k in sumkeys:
                if k in r:
                    countdata[key][k] = int(r[k])
        else:
            countdata[key]['Samples'] += ',%s' % (base,)
            for k in sumkeys:
                if k in r:
                    countdata[key][k] += int(r[k])

emptysym = ""
outheaders = countheaders
outheaders.insert(0, 'Samples')
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
pvalues = []
from fisher import fisher_exact, bonferroni, fdr, lod
progress.stage("Compute statistics")
for (snpstr, junc), r in sorted(countdata.items()):
    nsnpi = r.get('SNPJuncIntronCount', 0)
    nsnpex = r.get('SNPJuncNoIntronCount', 0)
    nwti = r.get('NoSNPJuncIntronCount', 0)
    nwtex = r.get('NoSNPJuncNoIntronCount', 0)

    p = emptysym
    pval = emptysym
    lodval = emptysym
    if junc and \
       (nsnpi + nsnpex) > 0 and \
       (nwti + nwtex) > 0 and \
       (nsnpi + nwti) > 0:

        psnp = nsnpi / float(nsnpi + nsnpex)
        pwt = nwti / float(nwti + nwtex)
        p = psnp / float(psnp + pwt)

        pval = fisher_exact(nsnpi,
                            (nsnpi + nsnpex),
                            (nsnpi + nwti),
                            (nsnpi + nsnpex + nwti + nwtex))

        pvalues.append(pval)

        l = lod(nsnpi, (nsnpi + nsnpex), (nsnpi + nwti),
                (nsnpi + nsnpex + nwti + nwtex))
        if l != None:
            lodval = l

    r['Probability'] = p
    r['P-Value'] = pval
    r['LOD'] = lodval

    row = dict(iter(r.items()))
    outrows.append(row)

bonf = bonferroni(pvalues)
fdr = fdr(pvalues)

i = 0
for r in outrows:
    if r['P-Value'] != emptysym:
        r['Bonferroni'] = bonf[i]
        r['FDR'] = fdr[i]
        i += 1

    for k in outheaders:
        if k not in r:
            r[k] = emptysym

progress.stage("Output results")
output.from_rows(outrows)
