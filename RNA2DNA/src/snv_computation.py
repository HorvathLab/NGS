#!/bin/env python27
import sys
import os
import csv
import os.path
import gzip
import re
from collections import defaultdict

from version import VERSION
VERSION = '1.0.6 (%s)' % (VERSION,)

from os.path import join, dirname, realpath, split
try:
    scriptdir = dirname(realpath(__file__))
except NameError:
    scriptdir = dirname(realpath(sys.argv[0]))
sys.path.append(join(scriptdir, '..', '..', 'common', 'src'))

from optparse_gui import OptionParser, OptionGroup, ProgressText
parser = OptionParser(version=VERSION)
regexs = OptionGroup(parser, "Filename Matching")

parser.add_option("--counts", type="file", dest="counts", default=None,
                  help="Output file from readCounts. Required.", notNone=True,
                  filetypes=[("readCount Output", "*.tsv")])
parser.add_option("--cosmic", type="file", dest="cosmic", default=None,
                  help="COSMIC Mutants.",
                  filetypes=[("COSMIC Annotations", "*.tsv;*.tsv.gz")])
parser.add_option("--darned", type="file", dest="darned", default=None,
                  help="DARNED annotations.",
                  filetypes=[("DARNED Annotations", "*.txt")])
regexs.add_option("--normaldnare", type="str", dest="normaldnare", default=r'GDNA',
                  help="Germline/Normal DNA filename regular expression. Default: GDNA.",
                  remember=True, name="Germline DNA RE")
regexs.add_option("--normaltransre", type="str", dest="normaltransre", default=r'NRNA',
                  help="Normal transcriptome filename regular expression. Default: NRNA.",
                  remember=True, name="Normal Transcr. RE")
regexs.add_option("--tumordnare", type="str", dest="tumordnare", default=r'SDNA',
                  help="Somatic/Tumor DNA filename regular expression. Default: SDNA.",
                  remember=True, name="Somatic DNA RE")
regexs.add_option("--tumortransre", type="str", dest="tumortransre", default=r'TRNA',
                  help="Tumor transcriptome filename regular expression. Default: TRNA.",
                  remember=True, name="Tumor Transcr. RE")
parser.add_option_group(regexs)                 

opt, args = parser.parse_args()
regex = {}
regex["GDNA"] = opt.normaldnare
regex["NRNA"] = opt.normaltransre
regex["SDNA"] = opt.tumordnare
regex["TRNA"] = opt.tumortransre

progress = ProgressText()

base = os.path.split(os.path.abspath(opt.counts))[0]

TRNA = {}; NRNA = {}; GDNA = {}; SDNA = {}

from chromreg import ChromLabelRegistry
chrreg = ChromLabelRegistry()
labels = map(str,range(1,100)) + ["X","Y","MT"]
chrreg.add_labels(opt.counts,labels)
chrreg.default_chrom_order()
chrorder = lambda l: chrreg.chrom_order(chrreg.label2chrom(opt.counts,l))

progress.stage("Parsing read-counts")
f = open(opt.counts, 'r')
reader = csv.DictReader(f, delimiter='\t')
types2files = defaultdict(set)
files2types = defaultdict(set)
for row in reader:
    key = (row['CHROM'],row['POS'])
    filename = row['AlignedReads']
    for k in row:
        if k.endswith('Count') and row[k] != "":
            row[k] = int(row[k])
        if k.endswith('Sc') and row[k] != "":
            row[k] = float(row[k])
    if re.search(regex["GDNA"],filename) and key not in GDNA:
        GDNA[key] = row; types2files["GDNA"].add(filename); files2types[filename].add("GDNA")
    if re.search(regex["NRNA"],filename) and key not in NRNA:
        NRNA[key] = row; types2files["NRNA"].add(filename); files2types[filename].add("NRNA")
    if re.search(regex["SDNA"],filename) and key not in SDNA:
        SDNA[key] = row; types2files["SDNA"].add(filename); files2types[filename].add("SDNA")
    if re.search(regex["TRNA"],filename) and key not in TRNA:
        TRNA[key] = row; types2files["TRNA"].add(filename); files2types[filename].add("TRNA")
f.close()
progress.done()

if sum(map(len,[GDNA,SDNA,NRNA,TRNA])) == 0:
    print >>sys.stderr, "No read counts available for testing"
    sys.exit(0)

fatal = False
for f in files2types:
    if len(files2types[f]) < 1:
        print >>sys.stderr, "Filename %s does not match any read type regular expression."%(f,)
	fatal = True
    elif len(files2types[f]) > 1:
        print >>sys.stderr, "Filename %s matches more than one read type regular expression."%(f,)
	print >>sys.stderr, "Matching regular expressions: %s."%(", ".join(map(regex.get,files2types[f])))
	fatal = True
for t in "GDNA NRNA SDNA TRNA".split():
    if len(types2files[t]) < 1:
	print >>sys.stderr, "Counts from %s %s reads are missing."%("normal" if t[0] in "NG" else "tumor",
								    "DNA" if t[1]=="D" else "RNA")
    elif len(types2files[t]) > 1:
	print >>sys.stderr, "Counts from %s %s reads are found in more than one file."%("normal" if t[0] in "NG" else "tumor",
								                        "DNA" if t[1]=="D" else "RNA")
	fatal = True

if fatal:
    sys.exit(1)

from event import *
sampsig = "".join(map(str,map(lambda t: 1*(len(types2files[t])>0),"GDNA SDNA NRNA TRNA".split())))
if sampsig == "1111":
    events = AllSamplesEvent
elif sampsig == "1100":
    events = DNAOnlyEvent
elif sampsig == "1010":
    events = NormalOnlyEvent
elif sampsig == "0111":
    events = NoGDNAEvent
else:
    raise RuntimeError("Bad combination of sample files")
progress.message("Testing for events: %s."%(", ".join(events.listall()),))

events.setCounts(GDNA,SDNA,NRNA,TRNA)

cosmic_headers = []
if opt.cosmic:
    progress.stage("Parsing COSMIC annotation file")
    if opt.cosmic.endswith('.gz'):
        f = gzip.open(opt.cosmic, 'r')
    else:
        f = open(opt.cosmic, 'r')
    reader = csv.DictReader(f, delimiter='\t')
    for cos in reader:
        if cos['Mutation genome position']:
	     chr,locus = cos['Mutation genome position'].split(':',1)
	     pos_st,pos_ed = locus.split('-',1)
	     if pos_st != pos_ed:
	         continue
	     key = (chr,pos_st)
             if key in events.keys:
                 for table in (GDNA,SDNA,NRNA,TRNA):
                     if key in table:
                         for k,h in zip(('Gene name','Primary site','Site subtype 1','Primary histology'),
                                        ('Gene','Site','Sub_Site','Cancer_Type')):
                             table[key]['COSMIC '+h] = cos[k]
    f.close()
    cosmic_headers = ['COSMIC Gene','COSMIC Site','COSMIC Sub_Site','COSMIC Cancer_Type']
    progress.done()

darned_headers = []
if opt.darned:
    progress.stage("Parsing DARNED annotation file")
    with open(opt.darned, 'r') as f:
        reader = csv.DictReader((f), delimiter='\t')
        for darn in reader:
	    if darn['chrom'] and darn['coordinate']:
                key = (darn['chrom'],darn['coordinate'])
                for table in (GDNA,SDNA,NRNA,TRNA):
                    if key in table:
                        for k,h in zip(('source',),
                                       ('Cancer_Type',)):
                            table[key]['DARNED '+h] = darn[k]
    darned_headers = ['DARNED Cancer_Type']
    progress.done()
    
progress.stage("Filtering SNVs")
events.testall()
progress.done()

headers = """
AlignedReads CHROM POS REF ALT
SNVCountForward SNVCountReverse RefCountForward RefCountReverse
SNVCount RefCount
HomoVarSc HetSc HomoRefSc
VarDomSc RefDomSc
NotHomoVarpV NotHomoRefpV NotHetpV VarDompV RefDompV
NotHomoVarFDR NotHomoRefFDR NotHetFDR VarDomFDR RefDomFDR
""".split()

progress.stage("Generating event reports")
for ev in events.events:
    wh = open(os.path.join(base,'Events_%s.tsv'%ev.abbrev),'w')
    event_headers = []
    event_headers.extend(headers)
    if ev.cosmic:
        event_headers.extend(cosmic_headers)
    if ev.darned:
        event_headers.extend(darned_headers)
    writer = csv.DictWriter(wh,fieldnames=event_headers,
                            extrasaction='ignore',dialect='excel-tab')
    writer.writeheader()
    for k in sorted(ev.goodkeys,key=lambda t: (chrorder(t[0]),int(t[1]))):
        for table in (GDNA,SDNA,NRNA,TRNA):
            if k in table:
                writer.writerow(table[k])
    wh.close()
progress.done()
    
