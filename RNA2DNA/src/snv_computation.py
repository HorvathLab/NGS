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

from optparse_gui import OptionParser, OptionGroup
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
regexs.add_option("--normalexomere", type="str", dest="normalexomere", default=r'GDNA',
                  help="Germline/Normal exome filename regular expression. Default: GDNA.",
                  remember=True, name="Germline Exome RE")
regexs.add_option("--normaltransre", type="str", dest="normaltransre", default=r'NRNA',
                  help="Normal transcriptome filename regular expression. Default: NRNA.",
                  remember=True, name="Normal Transcr. RE")
regexs.add_option("--tumorexomere", type="str", dest="tumorexomere", default=r'SDNA',
                  help="Somatic/Tumor exome filename regular expression. Default: SDNA.",
                  remember=True, name="Somatic Exome RE")
regexs.add_option("--tumortransre", type="str", dest="tumortransre", default=r'TRNA',
                  help="Tumor transcriptome filename regular expression. Default: TRNA.",
                  remember=True, name="Tumor Transcr. RE")
parser.add_option_group(regexs)                 

opt, args = parser.parse_args()
opt.normalexomere = re.compile(opt.normalexomere)
opt.normaltransre = re.compile(opt.normaltransre)
opt.tumorexomere = re.compile(opt.tumorexomere)
opt.tumortransre = re.compile(opt.tumortransre)

base = os.path.split(os.path.abspath(opt.counts))[0]

TRNA = {}; NRNA = {}; GDNA = {}; SDNA = {}

from chromreg import ChromLabelRegistry
chrreg = ChromLabelRegistry()
labels = map(str,range(1,100)) + ["X","Y","MT"]
chrreg.add_labels(opt.counts,labels)
chrreg.default_chrom_order()
chrorder = chrreg.chrom_order

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
    if opt.normalexomere.search(filename) and key not in GDNA:
        GDNA[key] = row; types2files["GDNA"].add(filename); files2types[filename].add("GDNA")
    if opt.normaltransre.search(filename) and key not in NRNA:
        NRNA[key] = row; types2files["NRNA"].add(filename); files2types[filename].add("NRNA")
    if opt.tumorexomere.search(filename) and key not in SDNA:
        SDNA[key] = row; types2files["SDNA"].add(filename); files2types[filename].add("SDNA")
    if opt.tumortransre.search(filename) and key not in TRNA:
        TRNA[key] = row; types2files["TRNA"].add(filename); files2types[filename].add("TRNA")
f.close()
fatal = False
for f in files2types:
    if len(files2types[f]) < 1:
        print >>sys.stderr, "Filename %s does not match any read type regular expression."%(f,)
	fatal = True
    elif len(files2types[f]) > 1:
        print >>sys.stderr, "Filename %s matches more than one read type regular expression."%(f,)
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
elif sampsig == "0011":
    events = RNAOnlyEvent
elif sampsig == "1010":
    events = NormalOnlyEvent
elif sampsig == "0101":
    events = TumorOnlyEvent
elif sampsig == "0111":
    events = NoGDNAEvent
else:
    raise RuntimeError("Bad combination of sample files")

events.setCounts(GDNA,SDNA,NRNA,TRNA)

cosmic_headers = []
if opt.cosmic:
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
                         for k,h in zip(('Gene name','Primary site','Site subtype1','Primary histology'),
                                        ('Gene','Site','Sub_Site','Cancer_Type')):
                             table[key]['COSMIC '+h] = cos[k]
    f.close()
    cosmic_headers = ['COSMIC Gene','COSMIC Site','COSMIC Sub_Site','COSMIC Cancer_Type']

darned_headers = []
if opt.darned:
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
    
events.testall()

headers = """
AlignedReads CHROM POS REF ALT
SNVCountForward SNVCountReverse RefCountForward RefCountReverse
SNVCount RefCount
HomoVarSc HetSc HomoRefSc
VarDomSc RefDomSc
NotHomoVarpV NotHomoRefpV NotHetpV VarDompV RefDompV
NotHomoVarFDR NotHomoRefFDR NotHetFDR VarDomFDR RefDomFDR
""".split()

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
    
sys.exit(0)    





def events(TRNA, NRNA, GDNA, SDNA):
    RNA_edit(TRNA, NRNA, GDNA, SDNA)
    Tum_RNA_edit(TRNA, NRNA, GDNA, SDNA)
    VSE(TRNA, NRNA, GDNA, SDNA)
    Tum_VSE(TRNA, NRNA, GDNA, SDNA)
    VSL(TRNA, NRNA, GDNA, SDNA)
    Tum_VSL(TRNA, NRNA, GDNA, SDNA)
    TSS_event(TRNA, NRNA, GDNA, SDNA)
    LOH(TRNA, NRNA, GDNA, SDNA)

def Som(GDNA, SDNA):
    
    with open(base + 'Events_SOM.tsv', 'w') as outpt_TSS:
        outpt_TSS.write('\t'.join(header_up) + '\t' + 'Gene' + '\t' +
                        'Site' + '\t' + 'Sub_Site' + '\t' + 'Cancer_Type' + '\n')
        keys = GDNA.keys()
        for k in sorted(keys):
            chr_sample = GDNA[k]['CHROM']
            pos_sample = GDNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(GDNA[k]['HomoRefSc']) >= 10:
                if float(SDNA[k]['HomoVarSc']) >= 26 or float(SDNA[k]['HetSc']) >= 26:
                    colum_GDNA = [GDNA[k][h] for h in header_up]
                    colum_SDNA = [SDNA[k][h] for h in header_up]
                    if cosmic.has_key(key_sample) == True:
                        Gene = cosmic[key_sample][0][0]
                        site = cosmic[key_sample][0][1]
                        sub_site = cosmic[key_sample][0][2]
                        cancer_ty = cosmic[key_sample][0][3]
                    else:
                        Gene = "NA"
                        site = "NA"
                        sub_site = "NA"
                        cancer_ty = "NA"
                    outpt_TSS.write('\t'.join(colum_GDNA) + '\t' + Gene +
                                    '\t' + site + '\t' + sub_site + '\t' + cancer_ty + '\n')
                    outpt_TSS.write('\t'.join(colum_SDNA) + '\t' + Gene +
                                    '\t' + site + '\t' + sub_site + '\t' + cancer_ty + '\n')


def TSS_event(TRNA, NRNA, GDNA, SDNA):
    with open(base + 'Events_SOM.tsv', 'w') as outpt_TSS:
        outpt_TSS.write('\t'.join(header_up) + '\t' + 'Gene' + '\t' +
                        'Site' + '\t' + 'Sub_Site' + '\t' + 'Cancer_Type' + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            GDNA_SNV_count = int(GDNA[k]['SNVCount'])
            GDNA_Ref_count = int(GDNA[k]['RefCount'])
            SDNA_SNV_count = int(SDNA[k]['SNVCount'])
            SDNA_Ref_count = int(SDNA[k]['RefCount'])
            NRNA_SNV_count = int(NRNA[k]['SNVCount'])
            NRNA_Ref_count = int(NRNA[k]['RefCount'])
            TRNA_SNV_count = int(TRNA[k]['SNVCount'])
            TRNA_Ref_count = int(TRNA[k]['RefCount'])
            if float(GDNA[k]['HomoRefSc']) >= 10:
                if float(SDNA[k]['HomoVarSc']) >= 26 or float(SDNA[k]['HetSc']) >= 26:
                    if float(NRNA[k]['HomoRefSc']) >= 5:
                        if float(TRNA[k]['HomoVarSc']) >= 5 or float(TRNA[k]['HetSc']) >= 5:
                            colum_GDNA = [GDNA[k][h] for h in header_up]
                            colum_SDNA = [SDNA[k][h] for h in header_up]
                            colum_NRNA = [NRNA[k][h] for h in header_up]
                            colum_TRNA = [TRNA[k][h] for h in header_up]
                            if cosmic.has_key(key_sample) == True:
                                Gene = cosmic[key_sample][0][0]
                                site = cosmic[key_sample][0][1]
                                sub_site = cosmic[key_sample][0][2]
                                cancer_ty = cosmic[key_sample][0][3]
                            else:
                                Gene = "NA"
                                site = "NA"
                                sub_site = "NA"
                                cancer_ty = "NA"
                            writefile_TSS(colum_GDNA, colum_SDNA, colum_NRNA,
                                          colum_TRNA, Gene, site, sub_site, cancer_ty, outpt_TSS)

def RNA_Ed_Sug(SDNA, TRNA):
      with open(base + 'Events_RNAed.tsv', 'w') as outpt_RNA:
        outpt_RNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            SDNA_SNV_count = int(SDNA[k]['SNVCount'])
            SDNA_Ref_count = int(SDNA[k]['RefCount'])
            TRNA_SNV_count = int(TRNA[k]['SNVCount'])
            TRNA_Ref_count = int(TRNA[k]['RefCount'])
            if (SDNA_Ref_count >= 3 and SDNA_SNV_count == 0 and TRNA_SNV_count >= 3):
                if float(SDNA[k]['HomoRefSc']) >= 85 and float(TRNA[k]['HetSc']) >= 25:
                    colum_SDNA = [SDNA[k][h] for h in header_up]
                    colum_TRNA = [TRNA[k][h] for h in header_up]
                    if darned.has_key(key_sample) == True:
                        cancer_ty = ','.join(darned[key_sample])
                    else:
                        cancer_ty = "NA"
                    outpt_RNA.write('\t'.join(colum_SDNA) +
                                '\t' + cancer_ty + '\n')
                    outpt_RNA.write('\t'.join(colum_TRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt_RNA.write('\n')


def RNA_Two_Edit(GDNA, NRNA):

    with open(base + 'Events_RNAed.tsv', 'w') as outpt_RNA:
        outpt_RNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = NRNA.keys()
        for k in sorted(keys):
            chr_sample = NRNA[k]['CHROM']
            pos_sample = NRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            GDNA_SNV_count = int(GDNA[k]['SNVCount'])
            GDNA_Ref_count = int(GDNA[k]['RefCount'])
            NRNA_SNV_count = int(NRNA[k]['SNVCount'])
            NRNA_Ref_count = int(NRNA[k]['RefCount'])
            if (GDNA_Ref_count >= 3 and GDNA_SNV_count == 0 and NRNA_SNV_count >= 3):
                if float(GDNA[k]['HomoRefSc']) >= 85 and float(NRNA[k]['HetSc']) >= 25:
                    colum_GDNA = [GDNA[k][h] for h in header_up]
                    colum_NRNA = [NRNA[k][h] for h in header_up]
                    if darned.has_key(key_sample) == True:
                        cancer_ty = ','.join(darned[key_sample])
                    else:
                        cancer_ty = "NA"
                    outpt_RNA.write('\t'.join(colum_GDNA) +
                                '\t' + cancer_ty + '\n')
                    outpt_RNA.write('\t'.join(colum_NRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt_RNA.write('\n')


def RNA_edit(TRNA, NRNA, GDNA, SDNA):
    with open(base + 'Events_RNAed.tsv', 'w') as outpt_RNA:
        outpt_RNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            GDNA_SNV_count = int(GDNA[k]['SNVCount'])
            GDNA_Ref_count = int(GDNA[k]['RefCount'])
            SDNA_SNV_count = int(SDNA[k]['SNVCount'])
            SDNA_Ref_count = int(SDNA[k]['RefCount'])
            NRNA_SNV_count = int(NRNA[k]['SNVCount'])
            NRNA_Ref_count = int(NRNA[k]['RefCount'])
            TRNA_SNV_count = int(TRNA[k]['SNVCount'])
            TRNA_Ref_count = int(TRNA[k]['RefCount'])
    #        if (GDNA_Ref_count >= 3 and GDNA_SNV_count == 0 and SDNA_Ref_count >= 3 and SDNA_SNV_count == 0):
     #           if (TRNA_SNV_count >= 3 and NRNA_SNV_count >= 3):
            if float(GDNA[k]['HomoRefSc']) >= 85 and float(SDNA[k]['HomoRefSc']) >= 85:
																											
                        if float(TRNA[k]['HetSc']) >= 25 or float(TRNA[k]['HomoVarSc']) >= 25:
                            if float(NRNA[k]['HetSc']) >= 20 or float(NRNA[k]['HomoVarSc']) >= 20:
                                colum_GDNA = [GDNA[k][h] for h in header_up]
                                colum_SDNA = [SDNA[k][h] for h in header_up]
                                colum_NRNA = [NRNA[k][h] for h in header_up]
                                colum_TRNA = [TRNA[k][h] for h in header_up]
                                if darned.has_key(key_sample) == True:
                                    cancer_ty = ','.join(darned[key_sample])
                                else:
                                    cancer_ty = "NA"
                                writefile_RNA_Edit(
                                    colum_GDNA, colum_SDNA, colum_NRNA, colum_TRNA, cancer_ty, outpt_RNA)


def T_RNA_Ed(NRNA, TRNA):
    with open(base + 'Events_T-RNAed.tsv', 'w') as outpt_tRNA:
        outpt_tRNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HetSc']) >= 20:
                if float(NRNA[k]['HomoRefSc']) >= 20:
                    colum_NRNA = [NRNA[k][h] for h in header_up]
                    colum_TRNA = [TRNA[k][h] for h in header_up]
                    if darned.has_key(key_sample) == True:
                        cancer_ty = ','.join(darned[key_sample])
                    else:
                        cancer_ty = "NA"
                    outpt.write('\t'.join(colum_NRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\t'.join(colum_TRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\n')


def Tum_RNA_edit(TRNA, NRNA, GDNA, SDNA):
    with open(base + 'Events_T-RNAed.tsv', 'w') as outpt_tRNA:
        outpt_tRNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
	    #print "k", k
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            #GDNA_SNV_count = int(GDNA[k]['SNVCount']) ;GDNA_Ref_count = int(GDNA[k]['RefCount'])
            #SDNA_SNV_count = int(SDNA[k]['SNVCount']) ;SDNA_Ref_count = int(SDNA[k]['RefCount'])
            #NRNA_SNV_count = int(NRNA[k]['SNVCount']) ;NRNA_Ref_count = int(NRNA[k]['RefCount'])
            #TRNA_SNV_count = int(TRNA[k]['SNVCount']) ;TRNA_Ref_count = int(TRNA[k]['RefCount'])
            # if (GDNA_Ref_count >= 3 and GDNA_SNV_count == 0 and SDNA_Ref_count >= 3 and SDNA_SNV_count==0):
            # if (NRNA_Ref_count >= 3 and NRNA_SNV_count == 0):
            # if (TRNA_SNV_count>= 3):
	    #print "key_sample",key_sample
	    #print TRNA[k]['HetSc']
																		
            if float(TRNA[k]['HomoVarSc']) >= 20 or float(TRNA[k]['HetSc']) >= 20:
	#	print NRNA[k]['HomoRefSc']
                if float(NRNA[k]['HomoRefSc']) >= 9:
                    if float(GDNA[k]['HomoRefSc']) >= 20:
                        if float(SDNA[k]['HomoRefSc']) >= 20:
                            colum_GDNA = [GDNA[k][h] for h in header_up]
                            colum_SDNA = [SDNA[k][h] for h in header_up]
                            colum_NRNA = [NRNA[k][h] for h in header_up]
                            colum_TRNA = [TRNA[k][h] for h in header_up]
                            if darned.has_key(key_sample) == True:
                                cancer_ty = ','.join(darned[key_sample])
                            else:
                                cancer_ty = "NA"
                            writefile_RNA_Edit(
                                colum_GDNA, colum_SDNA, colum_NRNA, colum_TRNA, cancer_ty, outpt_tRNA)

def T_RNA_Ed_3(SDNA,NRNA, TRNA):
      with open(base + 'Events_T-RNAed.tsv', 'w') as outpt_tRNA:
        outpt_tRNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HomoVarSc']) >= 20 or float(TRNA[k]['HetSc']) >= 20:
                if float(NRNA[k]['HomoRefSc']) >= 9:
                    if float(GDNA[k]['HomoRefSc']) >= 20: #<-!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if float(SDNA[k]['HomoRefSc']) >= 20:
                            colum_SDNA = [SDNA[k][h] for h in header_up]
                            colum_NRNA = [NRNA[k][h] for h in header_up]
                            colum_TRNA = [TRNA[k][h] for h in header_up]
                            if darned.has_key(key_sample) == True:
                                cancer_ty = ','.join(darned[key_sample])
                            else:
                                cancer_ty = "NA"
			    outpt.write('\t'.join(colum_SDNA) + '\n')
			    outpt.write('\t'.join(colum_NRNA) + '\n')
                	    outpt.write('\t'.join(colum_TRNA) + '\n')
                	    outpt.write('\n')


def VSE_Sug(SDNA, TRNA):

    with open(base + 'Events_VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = str(chr_sample) + ":" + pos_sample
            #GDNA_SNV_count = int(GDNA[k]['SNVCount']) ;GDNA_Ref_count = int(GDNA[k]['RefCount'])
            #NRNA_SNV_count = int(NRNA[k]['SNVCount']) ;NRNA_Ref_count = int(NRNA[k]['RefCount'])

            if float(TRNA[k]['HomoVarSc']) >= 27 and float(SDNA[k]['HetSc']) >= 27:
                colum_SDNA = [SDNA[k][h] for h in header_up]
                colum_TRNA = [TRNA[k][h] for h in header_up]
                outpt.write('\t'.join(colum_SDNA) + '\n')
                outpt.write('\t'.join(colum_TRNA) + '\n')
                outpt.write('\n')


def VSE_Two(GDNA, NRNA):
    with open(base + 'Events_VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = NRNA.keys()
        for k in sorted(keys):
            chr_sample = NRNA[k]['CHROM']
            pos_sample = NRNA[k]['POS']
            key_sample = str(chr_sample) + ":" + pos_sample
            #GDNA_SNV_count = int(GDNA[k]['SNVCount']) ;GDNA_Ref_count = int(GDNA[k]['RefCount'])
            #NRNA_SNV_count = int(NRNA[k]['SNVCount']) ;NRNA_Ref_count = int(NRNA[k]['RefCount'])

            if float(NRNA[k]['HomoVarSc']) >= 27 and float(GDNA[k]['HetSc']) >= 27:
                colum_GDNA = [GDNA[k][h] for h in header_up]
                colum_NRNA = [NRNA[k][h] for h in header_up]
                outpt.write('\t'.join(colum_GDNA) + '\n')
                outpt.write('\t'.join(colum_NRNA) + '\n')
                outpt.write('\n')


def VSE(TRNA, NRNA, GDNA, SDNA):
    with open(base + 'Events_VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = str(chr_sample) + ":" + pos_sample
            GDNA_SNV_count = int(GDNA[k]['SNVCount'])
            GDNA_Ref_count = int(GDNA[k]['RefCount'])
            SDNA_SNV_count = int(SDNA[k]['SNVCount'])
            SDNA_Ref_count = int(SDNA[k]['RefCount'])
            NRNA_SNV_count = int(NRNA[k]['SNVCount'])
            NRNA_Ref_count = int(NRNA[k]['RefCount'])
            TRNA_SNV_count = int(TRNA[k]['SNVCount'])
            TRNA_Ref_count = int(TRNA[k]['RefCount'])
            if (float(GDNA_SNV_count) + float(GDNA_Ref_count)) != 0.0 and (float(SDNA_SNV_count) + float(SDNA_Ref_count)) != 0.0:
                if (float(NRNA_SNV_count) + float(NRNA_Ref_count)) != 0.0 and float(TRNA_SNV_count) + float(TRNA_Ref_count) != 0.0:
                    ratioGDNA = float(
                        GDNA_SNV_count) / (float(GDNA_SNV_count) + float(GDNA_Ref_count))
                    ratioSDNA = float(
                        SDNA_SNV_count) / (float(SDNA_SNV_count) + float(SDNA_Ref_count))
                    ratioNRNA = float(
                        NRNA_SNV_count) / (float(NRNA_SNV_count) + float(NRNA_Ref_count))
                    ratioTRNA = float(
                        TRNA_SNV_count) / (float(TRNA_SNV_count) + float(TRNA_Ref_count))
            # if 0.4<=ratioGDNA <0.6 and 0.4<=ratioSDNA <0.6 and 0.8<=ratioNRNA<1 and 0.8<ratioTRNA <1:
                # print "i am in"
            if float(TRNA[k]['HomoVarSc']) >= 27 and float(NRNA[k]['HomoVarSc']) >= 27:
                # if float(GDNA[k]['VarDomSc']) <= 15 and float(SDNA[k]['VarDomSc']) <=15:
                if float(GDNA[k]['HetSc']) >= 27 and float(SDNA[k]['HetSc']) >= 27:
                    # and float(GDNA[k]['RefDomSc']) <= 15 and
                    # float(SDNA[k]['RefDomSc']) <=15:
                    colum_GDNA = [GDNA[k][h] for h in header_up]
                    colum_SDNA = [SDNA[k][h] for h in header_up]
                    colum_NRNA = [NRNA[k][h] for h in header_up]
                    colum_TRNA = [TRNA[k][h] for h in header_up]
                    writefile(colum_GDNA, colum_SDNA,
                              colum_NRNA, colum_TRNA, outpt)


def TS_VSE(NRNA, TRNA):
				
    with open(base + 'Events_T-VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HomoVarSc']) >= 27 and float(NRNA[k]['HetSc']) >= 27:
                colum_NRNA = [NRNA[k][h] for h in header_up]
                colum_TRNA = [TRNA[k][h] for h in header_up]
                outpt.write('\t'.join(colum_NRNA) + '\t' + Gene + '\t' +
                            site + '\t' + sub_site + '\t' + cancer_ty + '\n')
                outpt.write('\t'.join(colum_TRNA) + '\t' + Gene + '\t' +
                            site + '\t' + sub_site + '\t' + cancer_ty + '\n')
                outpt.write('\n')

def TS_VSE_3(SDNA,NRNA, TRNA):
    with open(base + 'Events_T-VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HomoVarSc']) >= 27 and float(NRNA[k]['HetSc']) >= 27:
                if float(GDNA[k]['HetSc']) >= 27 and float(SDNA[k]['HetSc']) >= 27: #<!--------
                    colum_SDNA = [SDNA[k][h] for h in header_up]
                    colum_NRNA = [NRNA[k][h] for h in header_up]
                    colum_TRNA = [TRNA[k][h] for h in header_up]
		    outpt.write('\t'.join(colum_SDNA) + '\n')
		    outpt.write('\t'.join(colum_NRNA) + '\n')
                    outpt.write('\t'.join(colum_TRNA) + '\n')
                    outpt.write('\n')

def Tum_VSE(TRNA, NRNA, GDNA, SDNA):
    with open(base + 'Events_T-VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HomoVarSc']) >= 27 and float(NRNA[k]['HetSc']) >= 27:
                if float(GDNA[k]['HetSc']) >= 27 and float(SDNA[k]['HetSc']) >= 27:
                    colum_GDNA = [GDNA[k][h] for h in header_up]
                    colum_SDNA = [SDNA[k][h] for h in header_up]
                    colum_NRNA = [NRNA[k][h] for h in header_up]
                    colum_TRNA = [TRNA[k][h] for h in header_up]
                    writefile(colum_GDNA, colum_SDNA,
                              colum_NRNA, colum_TRNA, outpt)

def VSL_Sug(TRNA, SDNA):
    with open(base + 'Events_VSL.tsv', 'w') as outpt:
	#outpt.write( '\n')
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HomoRefSc']) >= 27 and float(SDNA[k]['HetSc']) >= 28:
                colum_SDNA = [SDNA[k][h] for h in header_up]
                colum_TRNA = [SDNA[k][h] for h in header_up]
                outpt.write('\t'.join(colum_SDNA) + '\n')
                outpt.write('\t'.join(colum_TRNA) + '\n')
                outpt.write('\n')

def VSL_Two(GDNA, NRNA):
    with open(base + 'Events_VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = NRNA.keys()
        for k in sorted(keys):
            chr_sample = NRNA[k]['CHROM']
            pos_sample = NRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(NRNA[k]['HomoRefSc']) >= 27 and float(GDNA[k]['HetSc']) >= 28:
                colum_GDNA = [GDNA[k][h] for h in header_up]
                colum_NRNA = [NRNA[k][h] for h in header_up]
                outpt.write('\t'.join(colum_NRNA) + '\n')
                outpt.write('\t'.join(colum_GDNA) + '\n')
                outpt.write('\n')


def VSL(TRNA, NRNA, GDNA, SDNA):
    with open(base + 'Events_VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HomoRefSc']) >= 27 and float(NRNA[k]['HomoRefSc']) >= 27:
                if float(GDNA[k]['HetSc']) >= 28 and float(SDNA[k]['HetSc']) >= 27:
                    colum_GDNA = [GDNA[k][h] for h in header_up]
                    colum_SDNA = [SDNA[k][h] for h in header_up]
                    colum_NRNA = [NRNA[k][h] for h in header_up]
                    colum_TRNA = [TRNA[k][h] for h in header_up]
                    writefile(colum_GDNA, colum_SDNA,
                              colum_NRNA, colum_TRNA, outpt)


def TS_VSL(NRNA, TRNA):
    with open(base + 'Events_T-VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HomoRefSc']) >= 27 and float(NRNA[k]['HetSc']) >= 27:
                colum_NRNA = [NRNA[k][h] for h in header_up]
                colum_TRNA = [TRNA[k][h] for h in header_up]
                outpt.write('\t'.join(colum_NRNA) + '\n')
                outpt.write('\t'.join(colum_TRNA) + '\n')
                outpt.write('\n')


def Tum_VSL(TRNA, NRNA, GDNA, SDNA):
    with open(base + 'Events_T-VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HomoRefSc']) >= 27 and float(NRNA[k]['HetSc']) >= 27:
                if float(GDNA[k]['HetSc']) >= 27 and float(SDNA[k]['HetSc']) >= 27:
                    colum_GDNA = [GDNA[k][h] for h in header_up]
                    colum_SDNA = [SDNA[k][h] for h in header_up]
                    colum_NRNA = [NRNA[k][h] for h in header_up]
                    colum_TRNA = [TRNA[k][h] for h in header_up]
                    writefile(colum_GDNA, colum_SDNA,
                              colum_NRNA, colum_TRNA, outpt)

def TS_VSL_3(SDNA,NRNA, TRNA):
      with open(base + 'Events_T-VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(TRNA[k]['HomoRefSc']) >= 27 and float(NRNA[k]['HetSc']) >= 27:
                if float(SDNA[k]['HetSc']) >= 27 and float(SDNA[k]['HetSc']) >= 27:
                    #colum_GDNA = [GDNA[k][h] for h in header_up]
                    colum_SDNA = [SDNA[k][h] for h in header_up]
                    colum_NRNA = [NRNA[k][h] for h in header_up]
                    colum_TRNA = [TRNA[k][h] for h in header_up]
		    outpt.write('\t'.join(colum_SDNA) + '\n')
		    outpt.write('\t'.join(colum_NRNA) + '\n')
                    outpt.write('\t'.join(colum_TRNA) + '\n')
                    outpt.write('\n')


def LOH_exome(GDNA, SDNA):
    with open(base + 'Events_LOH.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SDNA.keys()
        for k in sorted(keys):
            chr_sample = SDNA[k]['CHROM']
            pos_sample = SDNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(GDNA[k]['HetSc']) >= 23 and float(SDNA[k]['HomoVarSc']) >= 23:
                colum_GDNA = [GDNA[k][h] for h in header_up]
                colum_SDNA = [SDNA[k][h] for h in header_up]
                outpt.write('\t'.join(colum_GDNA) + '\n')
                outpt.write('\t'.join(colum_SDNA) + '\n')


def LOH(TRNA, NRNA, GDNA, SDNA):
    with open(base + 'Events_LOH.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = TRNA.keys()
        for k in sorted(keys):
            chr_sample = TRNA[k]['CHROM']
            pos_sample = TRNA[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            GDNA_SNV_count = int(GDNA[k]['SNVCount'])
            GDNA_Ref_count = int(GDNA[k]['RefCount'])
            SDNA_SNV_count = int(SDNA[k]['SNVCount'])
            SDNA_Ref_count = int(SDNA[k]['RefCount'])
            NRNA_SNV_count = int(NRNA[k]['SNVCount'])
            NRNA_Ref_count = int(NRNA[k]['RefCount'])
            TRNA_SNV_count = int(TRNA[k]['SNVCount'])
            TRNA_Ref_count = int(TRNA[k]['RefCount'])
            # if (GDNA_Ref_count >= 1 and GDNA_SNV_count >= 1 and SDNA_Ref_count == 0 and SDNA_SNV_count>=1):
            # if (NRNA_Ref_count >= 3 or NRNA_SNV_count >= 3 or NRNA_Ref_count == 0 or NRNA_SNV_count == 0 ):
            # if (TRNA_SNV_count>= 1 and TRNA_Ref_count== 0):
            if float(GDNA[k]['HetSc']) >= 23 and float(SDNA[k]['HomoVarSc']) >= 23:
                if (float(NRNA[k]['HetSc']) >= 5 or float(NRNA[k]['HomoVarSc']) >= 5 or float(NRNA[k]['HomoRefSc']) >= 5):
                    if float(TRNA[k]['HomoVarSc']) >= 5:
                        colum_GDNA = [GDNA[k][h] for h in header_up]
                        colum_SDNA = [SDNA[k][h] for h in header_up]
                        colum_NRNA = [NRNA[k][h] for h in header_up]
                        colum_TRNA = [TRNA[k][h] for h in header_up]
                        writefile(colum_GDNA, colum_SDNA,
                                  colum_NRNA, colum_TRNA, outpt)


def writefile_TSS(colum_GDNA, colum_SDNA, colum_NRNA, colum_TRNA, Gene, site, sub_site, cancer_ty, outpt):
    outpt.write('\t'.join(colum_GDNA) + '\t' + Gene + '\t' +
                site + '\t' + sub_site + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_SDNA) + '\t' + Gene + '\t' +
                site + '\t' + sub_site + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_NRNA) + '\t' + Gene + '\t' +
                site + '\t' + sub_site + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_TRNA) + '\t' + Gene + '\t' +
                site + '\t' + sub_site + '\t' + cancer_ty + '\n')
    outpt.write('\n')


def writefile_RNA_Edit(colum_GDNA, colum_SDNA, colum_NRNA, colum_TRNA, cancer_ty, outpt):
    outpt.write('\t'.join(colum_GDNA) + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_SDNA) + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_NRNA) + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_TRNA) + '\t' + cancer_ty + '\n')
    outpt.write('\n')


def writefile(colum_GDNA, colum_SDNA, colum_NRNA, colum_TRNA, outpt):
    outpt.write('\t'.join(colum_GDNA) + '\n')
    outpt.write('\t'.join(colum_SDNA) + '\n')
    outpt.write('\t'.join(colum_NRNA) + '\n')
    outpt.write('\t'.join(colum_TRNA) + '\n')
    outpt.write('\n')
#Def
if GDNA and SDNA and TRNA and NRNA:
    print "four samples def"
    def commons_dict(keys):
        sorted_keys = sorted(keys)
        for k in sorted_keys:
            if k in GDNA and k in SDNA:
                if k in TRNA and k in NRNA:
                    return (GDNA, SDNA, NRNA, TRNA)

    GDNA = commons_dict(keys)[0]
    SDNA = commons_dict(keys)[1]
    NRNA = commons_dict(keys)[2]
    TRNA = commons_dict(keys)[3]
    events(TRNA, NRNA, GDNA, SDNA)
#Def
if GDNA and SDNA and not TRNA and not NRNA :
# or (not TRNA or TRNA) and not NRNA:
    print "two with LOH and SOM def"
    def commons_dict(keys_exomes):
        sorted_keys = sorted(keys_exomes)
        for k in sorted_keys:
            if k in GDNA and k in SDNA:
                return (GDNA, SDNA)

    GDNA = commons_dict(keys)[0]
    SDNA = commons_dict(keys)[1]
    Som(GDNA, SDNA)
    LOH_exome(GDNA, SDNA)

# DEf
if GDNA and NRNA and not TRNA and not SDNA:
    print "Two samples with VSE,VSL, RNA_Ed"
    def commons_dict(keys):
        sorted_keys = sorted(keys)
        for k in sorted_keys:
            if k in GDNA and k in NRNA:
                return (GDNA, NRNA)
    GDNA = commons_dict(keys)[0]
    NRNA = commons_dict(keys)[1]
    VSL_Two(GDNA, NRNA)
    VSE_Two(GDNA, NRNA)
    RNA_Two_Edit(GDNA, NRNA)

#Sugested
if SDNA and TRNA and not GDNA and ( not NRNA or NRNA) :
    print "Two samples with VSE,VSL, RNA_ed. The Suggestions"
    def commons_dict(keys):
        sorted_keys = sorted(keys)
        for k in sorted_keys:
            if k in TRNA and k in SDNA:
                return (SDNA, TRNA)

    SDNA = commons_dict(keys)[0]
    TRNA = commons_dict(keys)[1]
    VSL_Sug(SDNA, TRNA)
    VSE_Sug(SDNA, TRNA)
    RNA_Ed_Sug(SDNA, TRNA)

#sugested
if SDNA and TRNA and GDNA and not NRNA:
    print "three samples with VSE,VSL,RNA_ed and tumors of such"
    def commons_dict(keys):
        sorted_keys = sorted(keys)
        for k in sorted_keys:
            if k in GDNA and k in TRNA and k in SDNA:
                return (GDNA, TRNA, SDNA)
    SDNA = commons_dict(keys)[2]
    TRNA = commons_dict(keys)[1]
    GDNA = commons_dict(keys)[0]
    VSL_Sug(SDNA, TRNA)
    VSE_Sug(SDNA, TRNA)
    RNA_Ed_Sug(SDNA, TRNA)
    TS_VSL_3(SDNA,TRNA, GDNA)
    TS_VSE_3(SDNA,TRNA, GDNA)
    T_RNA_Ed_3(SDNA,TRNA, GDNA)
