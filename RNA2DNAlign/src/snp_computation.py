#!/bin/env python27
import sys
import os
import csv
import os.path
import gzip
import re
from collections import defaultdict

SRNA_d = defaultdict(list)
NRNA_d = defaultdict(list)
NDNA_d = defaultdict(list)
SDNA_d = defaultdict(list)

def chrorder(chr):
    the_chr_split = chr.split(':')
    chr = the_chr_split[0]
    if chr != "X":
        if chr != "Y":
            return the_chr_split
    if chr == "X":
        return 23
    if chr == "Y":
        return 24
    return 1e+20

from version import VERSION
VERSION = '1.0.5 (%s)' % (VERSION,)

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

base = os.path.split(os.path.abspath(opt.counts))[0] + os.sep
header = []

f = open(opt.counts, 'r')
reader = csv.DictReader(f, delimiter='\t')
types2files = defaultdict(set)
files2types = {}
for row in reader:
    key = row['CHROM'] + ":" + row['POS']
    filename = row['AlignedReads']
    files2types[filename] = set()
    if opt.normalexomere.search(filename) and key not in NDNA_d:
        NDNA_d[key] = row; types2files["NDNA"].add(filename); files2types[filename].add("NDNA")
    if opt.normaltransre.search(filename) and key not in NRNA_d:
        NRNA_d[key] = row; types2files["NRNA"].add(filename); files2types[filename].add("NRNA")
    if opt.tumorexomere.search(filename) and key not in SDNA_d:
        SDNA_d[key] = row; types2files["SDNA"].add(filename); files2types[filename].add("SDNA")
    if opt.tumortransre.search(filename) and key not in SRNA_d:
        SRNA_d[key] = row; types2files["SRNA"].add(filename); files2types[filename].add("SRNA")
f.close()
fatal = False
for f in files2types:
    if len(files2types[f]) < 1:
        print >>sys.stderr, "Filename %s does not match any read type regular expression."%(f,)
	fatal = True
    elif len(files2types[f]) > 1:
        print >>sys.stderr, "Filename %s matches more than one read type regular expression."%(f,)
	fatal = True
for t in "NDNA NRNA SDNA SRNA".split():
    if len(types2files[t]) < 1:
	print >>sys.stderr, "Counts from %s %s reads are missing."%("normal" if t[0]=="N" else "tumor",
								    "DNA" if t[1]=="D" else "RNA")
    elif len(types2files[t]) > 1:
	print >>sys.stderr, "Counts from %s %s reads are found in more than one file."%("normal" if t[0]=="N" else "tumor",
								                        "DNA" if t[1]=="D" else "RNA")
	fatal = True

if fatal:
    sys.exit(1)

header_up = "AlignedReads CHROM   POS   REF   ALT   SNPCountForward   SNPCountReverse   RefCountForward   RefCountReverse   SNPCount   RefCount HomoVarSc   HetSc   HomoRefSc   VarDomSc   RefDomSc NotHomoVarpV   NotHomoRefpV   NotHetpV   VarDompV   RefDompV   NotHomoVarFDR   NotHomoRefFDR   NotHetFDR   VarDomFDR   RefDomFDR".split()

keys = set()
keys.update(NDNA_d); keys.update(SDNA_d); 
keys.update(NRNA_d); keys.update(SRNA_d);

keys_exome = set()
keys_exome.update(NDNA_d); keys_exome.update(SDNA_d);

cosmic_dic = defaultdict(list)
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
	     key = chr + ":" + pos_st
	     if key in keys:
	         cosmic_dic[key].append((cos['Gene name'], cos['Primary site'], cos['Site subtype1'], cos['Primary histology']))
    f.close()

darn_dict = defaultdict(list)
if opt.darned:
    with open(opt.darned, 'r') as f:
        reader = csv.DictReader((f), delimiter='\t')
        for darn in reader:
	    if darn['chrom']:
                key = darn['chrom'] + ":" + darn['coordinate']
		if key in keys:
                    darn_dict[key].append(darn['source'])

def events(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    RNA_edit(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    Tum_RNA_edit(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    VSE(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    Tum_VSE(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    VSL(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    Tum_VSL(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    TSS_event(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    LOH(SRNA_d, NRNA_d, NDNA_d, SDNA_d)


def Som(NDNA_d, SDNA_d):
    cosmic_dic = defaultdict(list)
				
    with open(base + 'Events_SOM.tsv', 'w') as outpt_TSS:
        outpt_TSS.write('\t'.join(header_up) + '\t' + 'Gene' + '\t' +
                        'Site' + '\t' + 'Sub_Site' + '\t' + 'Cancer_Type' + '\n')
        keys = NDNA_d.keys()
        for k in sorted(keys):
            chr_sample = NDNA_d[k]['CHROM']
            pos_sample = NDNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(NDNA_d[k]['HomoRefSc']) >= 10:
                if float(SDNA_d[k]['HomoVarSc']) >= 26 or float(SDNA_d[k]['HetSc']) >= 26:
                    colum_NDNA = [NDNA_d[k][h] for h in header_up]
                    colum_SDNA = [SDNA_d[k][h] for h in header_up]
                    if cosmic_dic.has_key(key_sample) == True:
                        Gene = cosmic_dic[key_sample][0][0]
                        site = cosmic_dic[key_sample][0][1]
                        sub_site = cosmic_dic[key_sample][0][2]
                        cancer_ty = cosmic_dic[key_sample][0][3]
                    else:
                        Gene = "NA"
                        site = "NA"
                        sub_site = "NA"
                        cancer_ty = "NA"
                    outpt_TSS.write('\t'.join(colum_NDNA) + '\t' + Gene +
                                    '\t' + site + '\t' + sub_site + '\t' + cancer_ty + '\n')
                    outpt_TSS.write('\t'.join(colum_SDNA) + '\t' + Gene +
                                    '\t' + site + '\t' + sub_site + '\t' + cancer_ty + '\n')


def TSS_event(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_SOM.tsv', 'w') as outpt_TSS:
        outpt_TSS.write('\t'.join(header_up) + '\t' + 'Gene' + '\t' +
                        'Site' + '\t' + 'Sub_Site' + '\t' + 'Cancer_Type' + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k]['SNPCount'])
            NDNA_Ref_count = int(NDNA_d[k]['RefCount'])
            SDNA_SNP_count = int(SDNA_d[k]['SNPCount'])
            SDNA_Ref_count = int(SDNA_d[k]['RefCount'])
            NRNA_SNP_count = int(NRNA_d[k]['SNPCount'])
            NRNA_Ref_count = int(NRNA_d[k]['RefCount'])
            SRNA_SNP_count = int(SRNA_d[k]['SNPCount'])
            SRNA_Ref_count = int(SRNA_d[k]['RefCount'])
            if float(NDNA_d[k]['HomoRefSc']) >= 10:
                if float(SDNA_d[k]['HomoVarSc']) >= 26 or float(SDNA_d[k]['HetSc']) >= 26:
                    if float(NRNA_d[k]['HomoRefSc']) >= 5:
                        if float(SRNA_d[k]['HomoVarSc']) >= 5 or float(SRNA_d[k]['HetSc']) >= 5:
                            colum_NDNA = [NDNA_d[k][h] for h in header_up]
                            colum_SDNA = [SDNA_d[k][h] for h in header_up]
                            colum_NRNA = [NRNA_d[k][h] for h in header_up]
                            colum_SRNA = [SRNA_d[k][h] for h in header_up]
                            if cosmic_dic.has_key(key_sample) == True:
                                Gene = cosmic_dic[key_sample][0][0]
                                site = cosmic_dic[key_sample][0][1]
                                sub_site = cosmic_dic[key_sample][0][2]
                                cancer_ty = cosmic_dic[key_sample][0][3]
                            else:
                                Gene = "NA"
                                site = "NA"
                                sub_site = "NA"
                                cancer_ty = "NA"
                            writefile_TSS(colum_NDNA, colum_SDNA, colum_NRNA,
                                          colum_SRNA, Gene, site, sub_site, cancer_ty, outpt_TSS)


def RNA_Ed(NDNA_d, NRNA_d):
    with open(base + 'Events_RNAed.tsv', 'w') as outpt_RNA:
        outpt_RNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = NRNA_d.keys()
        for k in sorted(keys):
            chr_sample = NRNA_d[k]['CHROM']
            pos_sample = NRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k]['SNPCount'])
            NDNA_Ref_count = int(NDNA_d[k]['RefCount'])
            NRNA_SNP_count = int(NRNA_d[k]['SNPCount'])
            NRNA_Ref_count = int(NRNA_d[k]['RefCount'])
            if (NDNA_Ref_count >= 3 and NDNA_SNP_count == 0 and NRNA_SNP_count >= 3):
                if float(NDNA_d[k]['HomoRefSc']) >= 85 and float(NRNA_d[k]['HetSc']) >= 25:
                    colum_NDNA = [NDNA_d[k][h] for h in header_up]
                    colum_NRNA = [NRNA_d[k][h] for h in header_up]
                    if darn_dict.has_key(key_sample) == True:
                        cancer_ty = ','.join(darn_dict[key_sample])
                    else:
                        cancer_ty = "NA"
                    outpt.write('\t'.join(colum_NDNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\t'.join(colum_NRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\n')


def RNA_edit(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_RNAed.tsv', 'w') as outpt_RNA:
        outpt_RNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k]['SNPCount'])
            NDNA_Ref_count = int(NDNA_d[k]['RefCount'])
            SDNA_SNP_count = int(SDNA_d[k]['SNPCount'])
            SDNA_Ref_count = int(SDNA_d[k]['RefCount'])
            NRNA_SNP_count = int(NRNA_d[k]['SNPCount'])
            NRNA_Ref_count = int(NRNA_d[k]['RefCount'])
            SRNA_SNP_count = int(SRNA_d[k]['SNPCount'])
            SRNA_Ref_count = int(SRNA_d[k]['RefCount'])
    #        if (NDNA_Ref_count >= 3 and NDNA_SNP_count == 0 and SDNA_Ref_count >= 3 and SDNA_SNP_count == 0):
     #           if (SRNA_SNP_count >= 3 and NRNA_SNP_count >= 3):
            if float(NDNA_d[k]['HomoRefSc']) >= 85 and float(SDNA_d[k]['HomoRefSc']) >= 85:
																											
                        if float(SRNA_d[k]['HetSc']) >= 25 or float(SRNA_d[k]['HomoVarSc']) >= 25:
                            if float(NRNA_d[k]['HetSc']) >= 20 or float(NRNA_d[k]['HomoVarSc']) >= 20:
                                colum_NDNA = [NDNA_d[k][h] for h in header_up]
                                colum_SDNA = [SDNA_d[k][h] for h in header_up]
                                colum_NRNA = [NRNA_d[k][h] for h in header_up]
                                colum_SRNA = [SRNA_d[k][h] for h in header_up]
                                if darn_dict.has_key(key_sample) == True:
                                    cancer_ty = ','.join(darn_dict[key_sample])
                                else:
                                    cancer_ty = "NA"
                                writefile_RNA_Edit(
                                    colum_NDNA, colum_SDNA, colum_NRNA, colum_SRNA, cancer_ty, outpt_RNA)


def T_RNA_Ed(NRNA_d, SRNA_d):
    with open(base + 'Events_T-RNAed.tsv', 'w') as outpt_tRNA:
        outpt_tRNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k]['HetSc']) >= 20:
                if float(NRNA_d[k]['HomoRefSc']) >= 20:
                    colum_NRNA = [NRNA_d[k][h] for h in header_up]
                    colum_SRNA = [SRNA_d[k][h] for h in header_up]
                    if darn_dict.has_key(key_sample) == True:
                        cancer_ty = ','.join(darn_dict[key_sample])
                    else:
                        cancer_ty = "NA"
                    outpt.write('\t'.join(colum_NRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\t'.join(colum_SRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\n')


def Tum_RNA_edit(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_T-RNAed.tsv', 'w') as outpt_tRNA:
        outpt_tRNA.write('\t'.join(header_up) + '\t' + 'Cancer_Type' + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
	    #print "k", k
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            #NDNA_SNP_count = int(NDNA_d[k]['SNPCount']) ;NDNA_Ref_count = int(NDNA_d[k]['RefCount'])
            #SDNA_SNP_count = int(SDNA_d[k]['SNPCount']) ;SDNA_Ref_count = int(SDNA_d[k]['RefCount'])
            #NRNA_SNP_count = int(NRNA_d[k]['SNPCount']) ;NRNA_Ref_count = int(NRNA_d[k]['RefCount'])
            #SRNA_SNP_count = int(SRNA_d[k]['SNPCount']) ;SRNA_Ref_count = int(SRNA_d[k]['RefCount'])
            # if (NDNA_Ref_count >= 3 and NDNA_SNP_count == 0 and SDNA_Ref_count >= 3 and SDNA_SNP_count==0):
            # if (NRNA_Ref_count >= 3 and NRNA_SNP_count == 0):
            # if (SRNA_SNP_count>= 3):
	    #print "key_sample",key_sample
	    #print SRNA_d[k]['HetSc']
																		
            if float(SRNA_d[k]['HomoVarSc']) >= 20 or float(SRNA_d[k]['HetSc']) >= 20:
	#	print NRNA_d[k]['HomoRefSc']
                if float(NRNA_d[k]['HomoRefSc']) >= 9:
                    if float(NDNA_d[k]['HomoRefSc']) >= 20:
                        if float(SDNA_d[k]['HomoRefSc']) >= 20:
                            colum_NDNA = [NDNA_d[k][h] for h in header_up]
                            colum_SDNA = [SDNA_d[k][h] for h in header_up]
                            colum_NRNA = [NRNA_d[k][h] for h in header_up]
                            colum_SRNA = [SRNA_d[k][h] for h in header_up]
                            if darn_dict.has_key(key_sample) == True:
                                cancer_ty = ','.join(darn_dict[key_sample])
                            else:
                                cancer_ty = "NA"
                            writefile_RNA_Edit(
                                colum_NDNA, colum_SDNA, colum_NRNA, colum_SRNA, cancer_ty, outpt_tRNA)


def VSE(NDNA_d, NRNA_d):
    with open(base + 'Events_VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = NRNA_d[k]['CHROM']
            pos_sample = NRNA_d[k]['POS']
            key_sample = str(chr_sample) + ":" + pos_sample
            #NDNA_SNP_count = int(NDNA_d[k]['SNPCount']) ;NDNA_Ref_count = int(NDNA_d[k]['RefCount'])
            #NRNA_SNP_count = int(NRNA_d[k]['SNPCount']) ;NRNA_Ref_count = int(NRNA_d[k]['RefCount'])

            if float(NRNA_d[k]['HomoVarSc']) >= 27 and float(NDNA_d[k]['HetSc']) >= 27:
                colum_NDNA = [NDNA_d[k][h] for h in header_up]
                colum_NRNA = [NRNA_d[k][h] for h in header_up]
                outpt.write('\t'.join(colum_NDNA) + '\n')
                outpt.write('\t'.join(colum_NRNA) + '\n')
                outpt.write('\n')


def VSE(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = str(chr_sample) + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k]['SNPCount'])
            NDNA_Ref_count = int(NDNA_d[k]['RefCount'])
            SDNA_SNP_count = int(SDNA_d[k]['SNPCount'])
            SDNA_Ref_count = int(SDNA_d[k]['RefCount'])
            NRNA_SNP_count = int(NRNA_d[k]['SNPCount'])
            NRNA_Ref_count = int(NRNA_d[k]['RefCount'])
            SRNA_SNP_count = int(SRNA_d[k]['SNPCount'])
            SRNA_Ref_count = int(SRNA_d[k]['RefCount'])
            if (float(NDNA_SNP_count) + float(NDNA_Ref_count)) != 0.0 and (float(SDNA_SNP_count) + float(SDNA_Ref_count)) != 0.0:
                if (float(NRNA_SNP_count) + float(NRNA_Ref_count)) != 0.0 and float(SRNA_SNP_count) + float(SRNA_Ref_count) != 0.0:
                    ratioNDNA = float(
                        NDNA_SNP_count) / (float(NDNA_SNP_count) + float(NDNA_Ref_count))
                    ratioSDNA = float(
                        SDNA_SNP_count) / (float(SDNA_SNP_count) + float(SDNA_Ref_count))
                    ratioNRNA = float(
                        NRNA_SNP_count) / (float(NRNA_SNP_count) + float(NRNA_Ref_count))
                    ratioSRNA = float(
                        SRNA_SNP_count) / (float(SRNA_SNP_count) + float(SRNA_Ref_count))
            # if 0.4<=ratioNDNA <0.6 and 0.4<=ratioSDNA <0.6 and 0.8<=ratioNRNA<1 and 0.8<ratioSRNA <1:
                # print "i am in"
            if float(SRNA_d[k]['HomoVarSc']) >= 27 and float(NRNA_d[k]['HomoVarSc']) >= 27:
                # if float(NDNA_d[k]['VarDomSc']) <= 15 and float(SDNA_d[k]['VarDomSc']) <=15:
                if float(NDNA_d[k]['HetSc']) >= 27 and float(SDNA_d[k]['HetSc']) >= 27:
                    # and float(NDNA_d[k]['RefDomSc']) <= 15 and
                    # float(SDNA_d[k]['RefDomSc']) <=15:
                    colum_NDNA = [NDNA_d[k][h] for h in header_up]
                    colum_SDNA = [SDNA_d[k][h] for h in header_up]
                    colum_NRNA = [NRNA_d[k][h] for h in header_up]
                    colum_SRNA = [SRNA_d[k][h] for h in header_up]
                    writefile(colum_NDNA, colum_SDNA,
                              colum_NRNA, colum_SRNA, outpt)


def TS_VSE(NRNA_d, SRNA_d):
				
    with open(base + 'Events_T-VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k]['HomoVarSc']) >= 27 and float(NRNA_d[k]['HetSc']) >= 27:
                colum_NRNA = [NRNA_d[k][h] for h in header_up]
                colum_SRNA = [SRNA_d[k][h] for h in header_up]
                outpt.write('\t'.join(colum_NRNA) + '\t' + Gene + '\t' +
                            site + '\t' + sub_site + '\t' + cancer_ty + '\n')
                outpt.write('\t'.join(colum_SRNA) + '\t' + Gene + '\t' +
                            site + '\t' + sub_site + '\t' + cancer_ty + '\n')
                outpt.write('\n')


def Tum_VSE(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_T-VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k]['HomoVarSc']) >= 27 and float(NRNA_d[k]['HetSc']) >= 27:
                if float(NDNA_d[k]['HetSc']) >= 27 and float(SDNA_d[k]['HetSc']) >= 27:
                    colum_NDNA = [NDNA_d[k][h] for h in header_up]
                    colum_SDNA = [SDNA_d[k][h] for h in header_up]
                    colum_NRNA = [NRNA_d[k][h] for h in header_up]
                    colum_SRNA = [SRNA_d[k][h] for h in header_up]
                    writefile(colum_NDNA, colum_SDNA,
                              colum_NRNA, colum_SRNA, outpt)


def VSL(NDNA_d, NRNA_d):
    with open(base + 'Events_VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = NRNA_d.keys()
        for k in sorted(keys):
            chr_sample = NRNA_d[k]['CHROM']
            pos_sample = NRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(NRNA_d[k]['HomoRefSc']) >= 27 and float(NDNA_d[k]['HetSc']) >= 28:
                colum_NDNA = [NDNA_d[k][h] for h in header_up]
                colum_NRNA = [NRNA_d[k][h] for h in header_up]
                outpt.write('\t'.join(colum_NRNA) + '\n')
                outpt.write('\t'.join(colum_NDNA) + '\n')
                outpt.write('\n')


def VSL(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k]['HomoRefSc']) >= 27 and float(NRNA_d[k]['HomoRefSc']) >= 27:
                if float(NDNA_d[k]['HetSc']) >= 28 and float(SDNA_d[k]['HetSc']) >= 27:
                    colum_NDNA = [NDNA_d[k][h] for h in header_up]
                    colum_SDNA = [SDNA_d[k][h] for h in header_up]
                    colum_NRNA = [NRNA_d[k][h] for h in header_up]
                    colum_SRNA = [SRNA_d[k][h] for h in header_up]
                    writefile(colum_NDNA, colum_SDNA,
                              colum_NRNA, colum_SRNA, outpt)


def TS_VSL(NRNA_d, SRNA_d):
    with open(base + 'Events_T-VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k]['HomoRefSc']) >= 27 and float(NRNA_d[k]['HetSc']) >= 27:
                colum_NRNA = [NRNA_d[k][h] for h in header_up]
                colum_SRNA = [SRNA_d[k][h] for h in header_up]
                outpt.write('\t'.join(colum_NRNA) + '\n')
                outpt.write('\t'.join(colum_SRNA) + '\n')
                outpt.write('\n')


def Tum_VSL(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_T-VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k]['HomoRefSc']) >= 27 and float(NRNA_d[k]['HetSc']) >= 27:
                if float(NDNA_d[k]['HetSc']) >= 27 and float(SDNA_d[k]['HetSc']) >= 27:
                    colum_NDNA = [NDNA_d[k][h] for h in header_up]
                    colum_SDNA = [SDNA_d[k][h] for h in header_up]
                    colum_NRNA = [NRNA_d[k][h] for h in header_up]
                    colum_SRNA = [SRNA_d[k][h] for h in header_up]
                    writefile(colum_NDNA, colum_SDNA,
                              colum_NRNA, colum_SRNA, outpt)


def LOH_exome(NDNA_d, SDNA_d):
    with open(base + 'Events_LOH.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SDNA_d.keys()
        for k in sorted(keys):
            chr_sample = SDNA_d[k]['CHROM']
            pos_sample = SDNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            if float(NDNA_d[k]['HetSc']) >= 23 and float(SDNA_d[k]['HomoVarSc']) >= 23:
                colum_NDNA = [NDNA_d[k][h] for h in header_up]
                colum_SDNA = [SDNA_d[k][h] for h in header_up]
                outpt.write('\t'.join(colum_NDNA) + '\n')
                outpt.write('\t'.join(colum_SDNA) + '\n')


def LOH(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_LOH.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k]['CHROM']
            pos_sample = SRNA_d[k]['POS']
            key_sample = chr_sample + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k]['SNPCount'])
            NDNA_Ref_count = int(NDNA_d[k]['RefCount'])
            SDNA_SNP_count = int(SDNA_d[k]['SNPCount'])
            SDNA_Ref_count = int(SDNA_d[k]['RefCount'])
            NRNA_SNP_count = int(NRNA_d[k]['SNPCount'])
            NRNA_Ref_count = int(NRNA_d[k]['RefCount'])
            SRNA_SNP_count = int(SRNA_d[k]['SNPCount'])
            SRNA_Ref_count = int(SRNA_d[k]['RefCount'])
            # if (NDNA_Ref_count >= 1 and NDNA_SNP_count >= 1 and SDNA_Ref_count == 0 and SDNA_SNP_count>=1):
            # if (NRNA_Ref_count >= 3 or NRNA_SNP_count >= 3 or NRNA_Ref_count == 0 or NRNA_SNP_count == 0 ):
            # if (SRNA_SNP_count>= 1 and SRNA_Ref_count== 0):
            if float(NDNA_d[k]['HetSc']) >= 23 and float(SDNA_d[k]['HomoVarSc']) >= 23:
                if (float(NRNA_d[k]['HetSc']) >= 5 or float(NRNA_d[k]['HomoVarSc']) >= 5 or float(NRNA_d[k]['HomoRefSc']) >= 5):
                    if float(SRNA_d[k]['HomoVarSc']) >= 5:
                        colum_NDNA = [NDNA_d[k][h] for h in header_up]
                        colum_SDNA = [SDNA_d[k][h] for h in header_up]
                        colum_NRNA = [NRNA_d[k][h] for h in header_up]
                        colum_SRNA = [SRNA_d[k][h] for h in header_up]
                        writefile(colum_NDNA, colum_SDNA,
                                  colum_NRNA, colum_SRNA, outpt)


def writefile_TSS(colum_NDNA, colum_SDNA, colum_NRNA, colum_SRNA, Gene, site, sub_site, cancer_ty, outpt):
    outpt.write('\t'.join(colum_NDNA) + '\t' + Gene + '\t' +
                site + '\t' + sub_site + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_SDNA) + '\t' + Gene + '\t' +
                site + '\t' + sub_site + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_NRNA) + '\t' + Gene + '\t' +
                site + '\t' + sub_site + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_SRNA) + '\t' + Gene + '\t' +
                site + '\t' + sub_site + '\t' + cancer_ty + '\n')
    outpt.write('\n')


def writefile_RNA_Edit(colum_NDNA, colum_SDNA, colum_NRNA, colum_SRNA, cancer_ty, outpt):
    outpt.write('\t'.join(colum_NDNA) + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_SDNA) + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_NRNA) + '\t' + cancer_ty + '\n')
    outpt.write('\t'.join(colum_SRNA) + '\t' + cancer_ty + '\n')
    outpt.write('\n')


def writefile(colum_NDNA, colum_SDNA, colum_NRNA, colum_SRNA, outpt):
    outpt.write('\t'.join(colum_NDNA) + '\n')
    outpt.write('\t'.join(colum_SDNA) + '\n')
    outpt.write('\t'.join(colum_NRNA) + '\n')
    outpt.write('\t'.join(colum_SRNA) + '\n')
    outpt.write('\n')

if NDNA_d and SDNA_d and SRNA_d and NRNA_d:
    def commons_dict(keys):
        sorted_keys = sorted(keys)
        for k in sorted_keys:
            if k in NDNA_d and k in SDNA_d:
                if k in SRNA_d and k in NRNA_d:
                    return (NDNA_d, SDNA_d, NRNA_d, SRNA_d)

    NDNA_d = commons_dict(keys)[0]
    SDNA_d = commons_dict(keys)[1]
    NRNA_d = commons_dict(keys)[2]
    SRNA_d = commons_dict(keys)[3]
    events(SRNA_d, NRNA_d, NDNA_d, SDNA_d)

if NDNA_d and SDNA_d and not SRNA_d and not NRNA_d:
    def commons_dict(keys_exomes):
        sorted_keys = sorted(keys_exomes)
        for k in sorted_keys:
            if k in NDNA_d and k in SDNA_d:
                return (NDNA_d, SDNA_d)

    NDNA_d = commons_dict(keys)[0]
    SDNA_d = commons_dict(keys)[1]
    Som(NDNA_d, SDNA_d)
    LOH_exome(NDNA_d, SDNA_d)


if NRNA_d and SRNA_d and not SDNA_d and not NDNA_d:
    def commons_dict(keys):
        sorted_keys = sorted(keys)
        for k in sorted_keys:
            if k in SRNA_d and k in NRNA_d:
                return (NRNA_d, SRNA_d)

    NRNA_d = commons_dict(keys)[0]
    SRNA_d = commons_dict(keys)[1]
    TS_VSL(NRNA_d, SRNA_d)
    TS_VSE(NRNA_d, SRNA_d)
    T_RNA_Ed(NRNA_d, SRNA_d)

if NDNA_d and NRNA_d and not SRNA_d and not SDNA_d:
    def commons_dict(keys):
        sorted_keys = sorted(keys)
        for k in sorted_keys:
            if k in NDNA_d and k in NRNA_d:
                return (NDNA_d, NRNA_d)
    NDNA_d = commons_dict(keys)[0]
    NRNA_d = commons_dict(keys)[1]
    VSL(NDNA_d, NRNA_d)
    VSE(NDNA_d, NRNA_d)
    RNA_Ed(NDNA_d, NRNA_d)
