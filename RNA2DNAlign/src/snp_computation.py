#!/bin/env python27
import sys
import os
import csv
import os.path
import gzip
import collections
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

from version import VERSION
VERSION = '1.0.0 (%s)' % (VERSION,)

from optparse_gui import OptionParser
parser = OptionParser(version=VERSION)

parser.add_option("--counts", type="file", dest="counts", default=None,
                  help="Output file from readCounts. Required.", notNone=True,
                  filetypes=[("readCount Output", "*.tsv")])
parser.add_option("--cosmic", type="file", dest="cosmic", default=None,
                  help="COSMIC Mutants.",
                  filetypes=[("COSMIC Annotations", "*.tsv;*.tsv.gz")])
parser.add_option("--darned", type="file", dest="darned", default=None,
                  help="DARNED annotations.",
                  filetypes=[("DARNED Annotations", "*.txt")])

opt, args = parser.parse_args()

base = os.path.split(os.path.abspath(opt.counts))[0] + os.sep
header = []

f = open(opt.counts, 'r')
reader = csv.reader(f, delimiter='\t')
for i, row in enumerate(reader):
    if "CHROM" in row[0]:
        header.append(row)
    else:
        if 'SRNA'in row[4] or 'TPtr'in row[4]:
            key1 = str(row[0]) + ":" + row[1]
            SRNA_d[key1].append(row)
        if 'NRNA' in row[4] or 'NTtr' in row[4]:
            key2 = str(row[0]) + ":" + row[1]
            NRNA_d[key2].append(row)
        if 'NDNA'in row[4] or 'NTex'in row[4]:
            key3 = str(row[0]) + ":" + row[1]
            NDNA_d[key3].append(row)
        if 'SDNA'in row[4] or 'TPex' in row[4]:
            key4 = str(row[0]) + ":" + row[1]
            SDNA_d[key4].append(row)
f.close()

range_column = [4, 0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 13, 14,
                15, 16, 17, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]


header_up = list(header[0][i] for i in range_column)
keys = set(SRNA_d.keys()) | set(NRNA_d.keys()) | set(
    NDNA_d.keys()) | set(SDNA_d.keys())

keys_exome = set(NDNA_d.keys()) | set(SDNA_d.keys())
cosmic_Mut = []
if opt.cosmic:
    if opt.cosmic.endswith('.gz'):
        f = gzip.open(opt.cosmic, 'r')
    else:
        f = open(opt.cosmic, 'r')
    reader = csv.reader(f, delimiter='\t')
    for i in reader:
        if 'Gene' not in i[0] and i[17] != '':
            cosmic_Mut.append(i)
    f.close()


def events(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    VSE(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    Tum_VSE(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    VSL(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    Tum_VSL(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    LOH(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    TSS_event(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    RNA_edit(SRNA_d, NRNA_d, NDNA_d, SDNA_d)
    Tum_RNA_edit(SRNA_d, NRNA_d, NDNA_d, SDNA_d)


def Som(NDNA_d, SDNA_d):
    cosmic_dic = defaultdict(list)

    with open(base + 'Events_TSS.tsv', 'w') as outpt_TSS:
        outpt_TSS.write('\t'.join(header_up) + '\t' + 'Gene' + '\t' +
                        'Site' + '\t' + 'Sub_Site' + '\t' + 'Cancer_Type' + '\n')
        for cos in cosmic_Mut:
            coord = cos[17].split(':')
            if coord:
                chr = coord[0]
                pos_un = coord[1].split('-')
                pos = pos_un[0]
                key_cs = chr + ":" + pos
                cosmic_dic[key_cs].append((cos[0], cos[7], cos[8], cos[9]))
        keys = NDNA_d.keys()
        for k in sorted(keys):
            chr_sample = NDNA_d[k][0][0]
            pos_sample = NDNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            if float(NDNA_d[k][0][15]) >= 10:
                if float(SDNA_d[k][0][13]) >= 26 or float(SDNA_d[k][0][14]) >= 26:
                    colum_NDNA = list(NDNA_d[k][0][i] for i in range_column)
                    colum_SDNA = list(SDNA_d[k][0][i] for i in range_column)
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
    cosmic_dic = defaultdict(list)
    with open(base + 'Events_TSS.tsv', 'w') as outpt_TSS:
        outpt_TSS.write('\t'.join(header_up) + '\t' + 'Gene' + '\t' +
                        'Site' + '\t' + 'Sub_Site' + '\t' + 'Cancer_Type' + '\n')
        for cos in cosmic_Mut:
            coord = cos[17].split(':')
            if coord:
                chr = coord[0]
                pos_un = coord[1].split('-')
                pos = pos_un[0]
                key_cs = chr + ":" + pos
                cosmic_dic[key_cs].append((cos[0], cos[7], cos[8], cos[9]))
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k][0][9])
            NDNA_Ref_count = int(NDNA_d[k][0][10])
            SDNA_SNP_count = int(SDNA_d[k][0][9])
            SDNA_Ref_count = int(SDNA_d[k][0][10])
            NRNA_SNP_count = int(NRNA_d[k][0][9])
            NRNA_Ref_count = int(NRNA_d[k][0][10])
            SRNA_SNP_count = int(SRNA_d[k][0][9])
            SRNA_Ref_count = int(SRNA_d[k][0][10])
            if float(NDNA_d[k][0][15]) >= 10:
                if float(SDNA_d[k][0][13]) >= 26 or float(SDNA_d[k][0][14]) >= 26:
                    if float(NRNA_d[k][0][15]) >= 5:
                        if float(SRNA_d[k][0][13]) >= 5 or float(SRNA_d[k][0][14]) >= 5:
                            colum_SRNA = list(SRNA_d[k][0][i]
                                              for i in range_column)
                            colum_NRNA = list(NRNA_d[k][0][i]
                                              for i in range_column)
                            colum_NDNA = list(NDNA_d[k][0][i]
                                              for i in range_column)
                            colum_SDNA = list(SDNA_d[k][0][i]
                                              for i in range_column)
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

darn_l = []
if opt.darned:
    with open(opt.darned, 'r') as f:
        reader = csv.reader((f), delimiter='\t')
        for i in reader:
            if 'chrom' not in i[0]:
                darn_l.append(i)


def RNA_Ed(NDNA_d, NRNA_d):
    darn_dict = defaultdict(list)
    with open(base + 'Events_RNA_Edit.tsv', 'w') as outpt_RNA:
        outpt_RNA.write('\t'.join(header_up) + 'Cancer_Type' + '\n')
        for darn in darn_l:
            key_darn = darn[0] + ":" + darn[1]
            darn_dict[key_darn].append(darn[8])
        keys = NRNA_d.keys()
        for k in sorted(keys):
            chr_sample = NRNA_d[k][0][0]
            pos_sample = NRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k][0][9])
            NDNA_Ref_count = int(NDNA_d[k][0][10])
            NRNA_SNP_count = int(NRNA_d[k][0][9])
            NRNA_Ref_count = int(NRNA_d[k][0][10])
            if (NDNA_Ref_count >= 3 and NDNA_SNP_count == 0 and NRNA_SNP_count >= 3):
                if float(NDNA_d[k][0][15]) >= 85 and float(NRNA_d[k][0][14]) >= 25:
                    colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                    colum_NDNA = list(NDNA_d[k][0][i] for i in range_column)
                    if darn_dict.has_key(key_sample) == True:
                        cancer_ty = ''.join(darn_dict[key_sample])
                    else:
                        cancer_ty = "NA"
                    outpt.write('\t'.join(colum_NRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\t'.join(colum_SRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\n')


def RNA_edit(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    darn_dict = defaultdict(list)
    with open(base + 'Events_RNA_Edit.tsv', 'w') as outpt_RNA:
        outpt_RNA.write('\t'.join(header_up) + 'Cancer_Type' + '\n')
        for darn in darn_l:
            key_darn = darn[0] + ":" + darn[1]
            darn_dict[key_darn].append(darn[8])
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k][0][9])
            NDNA_Ref_count = int(NDNA_d[k][0][10])
            SDNA_SNP_count = int(SDNA_d[k][0][9])
            SDNA_Ref_count = int(SDNA_d[k][0][10])
            NRNA_SNP_count = int(NRNA_d[k][0][9])
            NRNA_Ref_count = int(NRNA_d[k][0][10])
            SRNA_SNP_count = int(SRNA_d[k][0][9])
            SRNA_Ref_count = int(SRNA_d[k][0][10])
            if (NDNA_Ref_count >= 3 and NDNA_SNP_count == 0 and SDNA_Ref_count >= 3 and SDNA_SNP_count == 0):
                if (SRNA_SNP_count >= 3 and NRNA_SNP_count >= 3):
                    if float(NDNA_d[k][0][15]) >= 85 and float(SDNA_d[k][0][15]) >= 85:
                        if float(SRNA_d[k][0][14]) >= 25:
                            if float(NRNA_d[k][0][14]) >= 25:
                                colum_SRNA = list(SRNA_d[k][0][i]
                                                  for i in range_column)
                                colum_NRNA = list(NRNA_d[k][0][i]
                                                  for i in range_column)
                                colum_NDNA = list(NDNA_d[k][0][i]
                                                  for i in range_column)
                                colum_SDNA = list(SDNA_d[k][0][i]
                                                  for i in range_column)
                                if darn_dict.has_key(key_sample) == True:
                                    cancer_ty = ''.join(darn_dict[key_sample])
                                else:
                                    cancer_ty = "NA"
                                writefile_RNA_Edit(
                                    colum_NDNA, colum_SDNA, colum_NRNA, colum_SRNA, cancer_ty, outpt_RNA)


def T_RNA_Ed(NRNA_d, SRNA_d):
    darn_dict = defaultdict(list)
    with open(base + 'Events_Tum_RNA_Edit.tsv', 'w') as outpt_tRNA:
        outpt_tRNA.write('\t'.join(header_up) + 'Cancer_Type' + '\n')
        for darn in darn_l:
            key_darn = darn[0] + ":" + darn[1]
            darn_dict[key_darn].append(darn[8])
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k][0][14]) >= 20:
                if float(NRNA_d[k][0][15]) >= 20:
                    colum_SRNA = list(SRNA_d[k][0][i] for i in range_column)
                    colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                    if darn_dict.has_key(key_sample) == True:
                        cancer_ty = ''.join(darn_dict[key_sample])
                    else:
                        cancer_ty = "NA"
                    outpt.write('\t'.join(colum_NRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\t'.join(colum_SRNA) +
                                '\t' + cancer_ty + '\n')
                    outpt.write('\n')


def Tum_RNA_edit(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    darn_dict = defaultdict(list)
    with open(base + 'Events_Tum_RNA_Edit.tsv', 'w') as outpt_tRNA:
        outpt_tRNA.write('\t'.join(header_up) + 'Cancer_Type' + '\n')
        for darn in darn_l:
            key_darn = darn[0] + ":" + darn[1]
            darn_dict[key_darn].append(darn[8])
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            #NDNA_SNP_count = int(NDNA_d[k][0][9]) ;NDNA_Ref_count = int(NDNA_d[k][0][10])
            #SDNA_SNP_count = int(SDNA_d[k][0][9]) ;SDNA_Ref_count = int(SDNA_d[k][0][10])
            #NRNA_SNP_count = int(NRNA_d[k][0][9]) ;NRNA_Ref_count = int(NRNA_d[k][0][10])
            #SRNA_SNP_count = int(SRNA_d[k][0][9]) ;SRNA_Ref_count = int(SRNA_d[k][0][10])
            # if (NDNA_Ref_count >= 3 and NDNA_SNP_count == 0 and SDNA_Ref_count >= 3 and SDNA_SNP_count==0):
            # if (NRNA_Ref_count >= 3 and NRNA_SNP_count == 0):
            # if (SRNA_SNP_count>= 3):
            if float(SRNA_d[k][0][14]) >= 20:
                if float(NRNA_d[k][0][15]) >= 20:
                    if float(NDNA_d[k][0][15]) >= 20:
                        if float(SDNA_d[k][0][15]) >= 20:
                            colum_SRNA = list(SRNA_d[k][0][i]
                                              for i in range_column)
                            colum_NRNA = list(NRNA_d[k][0][i]
                                              for i in range_column)
                            colum_NDNA = list(NDNA_d[k][0][i]
                                              for i in range_column)
                            colum_SDNA = list(SDNA_d[k][0][i]
                                              for i in range_column)
                            if darn_dict.has_key(key_sample) == True:
                                cancer_ty = ''.join(darn_dict[key_sample])
                            else:
                                cancer_ty = "NA"
                            writefile_RNA_Edit(
                                colum_NDNA, colum_SDNA, colum_NRNA, colum_SRNA, cancer_ty, outpt_tRNA)


def VSE(NDNA_d, NRNA_d):
    with open(base + 'Events_VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = NRNA_d[k][0][0]
            pos_sample = NRNA_d[k][0][1]
            key_sample = str(chr_sample) + ":" + pos_sample
            #NDNA_SNP_count = int(NDNA_d[k][0][9]) ;NDNA_Ref_count = int(NDNA_d[k][0][10])
            #NRNA_SNP_count = int(NRNA_d[k][0][9]) ;NRNA_Ref_count = int(NRNA_d[k][0][10])

            if float(NRNA_d[k][0][13]) >= 27 and float(NDNA_d[k][0][14]) >= 27:
                colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                colum_NDNA = list(NDNA_d[k][0][i] for i in range_column)
                outpt.write('\t'.join(colum_NDNA) + '\n')
                outpt.write('\t'.join(colum_NRNA) + '\n')
                outpt.write('\n')


def VSE(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + '_Events_VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = str(chr_sample) + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k][0][9])
            NDNA_Ref_count = int(NDNA_d[k][0][10])
            SDNA_SNP_count = int(SDNA_d[k][0][9])
            SDNA_Ref_count = int(SDNA_d[k][0][10])
            NRNA_SNP_count = int(NRNA_d[k][0][9])
            NRNA_Ref_count = int(NRNA_d[k][0][10])
            SRNA_SNP_count = int(SRNA_d[k][0][9])
            SRNA_Ref_count = int(SRNA_d[k][0][10])
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
            if float(SRNA_d[k][0][13]) >= 27 and float(NRNA_d[k][0][13]) >= 27:
                # if float(NDNA_d[k][0][16]) <= 15 and float(SDNA_d[k][0][16])
                # <=15:
                if float(NDNA_d[k][0][14]) >= 27 and float(SDNA_d[k][0][14]) >= 27:
                    # and float(NDNA_d[k][0][17]) <= 15 and
                    # float(SDNA_d[k][0][17]) <=15:
                    colum_SRNA = list(SRNA_d[k][0][i] for i in range_column)
                    colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                    colum_NDNA = list(NDNA_d[k][0][i] for i in range_column)
                    colum_SDNA = list(SDNA_d[k][0][i] for i in range_column)
                    writefile(colum_NDNA, colum_SDNA,
                              colum_NRNA, colum_SRNA, outpt)


def TS_VSE(NRNA_d, SRNA_d):

    with open(base + 'Events_Tum_VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k][0][13]) >= 27 and float(NRNA_d[k][0][14]) >= 27:
                colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                colum_SRNA = list(SRNA_d[k][0][i] for i in range_column)
                outpt.write('\t'.join(colum_NRNA) + '\t' + Gene + '\t' +
                            site + '\t' + sub_site + '\t' + cancer_ty + '\n')
                outpt.write('\t'.join(colum_SRNA) + '\t' + Gene + '\t' +
                            site + '\t' + sub_site + '\t' + cancer_ty + '\n')
                outpt.write('\n')


def Tum_VSE(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_Tum_VSE.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k][0][13]) >= 27 and float(NRNA_d[k][0][14]) >= 27:
                if float(NDNA_d[k][0][14]) >= 27 and float(SDNA_d[k][0][14]) >= 27:
                    colum_SRNA = list(SRNA_d[k][0][i] for i in range_column)
                    #colum_SRNA = sorted(colum_SRNA, key = lambda  l: chrorder(l[1]))
                    colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                    colum_NDNA = list(NDNA_d[k][0][i] for i in range_column)
                    colum_SDNA = list(SDNA_d[k][0][i] for i in range_column)
                    writefile(colum_NDNA, colum_SDNA,
                              colum_NRNA, colum_SRNA, outpt)


def VSL(NDNA_d, NRNA_d):
    with open(base + 'Events_VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = NRNA_d.keys()
        for k in sorted(keys):
            chr_sample = NRNA_d[k][0][0]
            pos_sample = NRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            if float(NRNA_d[k][0][15]) >= 27 and float(NDNA_d[k][0][14]) >= 28:
                colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                colum_NDNA = list(NDNA_d[k][0][i] for i in range_column)
                outpt.write('\t'.join(colum_NRNA) + '\n')
                outpt.write('\t'.join(colum_NDNA) + '\n')
                outpt.write('\n')


def VSL(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k][0][15]) >= 27 and float(NRNA_d[k][0][15]) >= 27:
                if float(NDNA_d[k][0][14]) >= 28 and float(SDNA_d[k][0][14]) >= 27:
                    colum_SRNA = list(SRNA_d[k][0][i] for i in range_column)
                    colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                    colum_NDNA = list(NDNA_d[k][0][i] for i in range_column)
                    colum_SDNA = list(SDNA_d[k][0][i] for i in range_column)
                    writefile(colum_NDNA, colum_SDNA,
                              colum_NRNA, colum_SRNA, outpt)


def TS_VSL(NRNA_d, SRNA_d):
    with open(base + 'Events_Tum_VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k][0][15]) >= 27 and float(NRNA_d[k][0][14]) >= 27:
                colum_SRNA = list(SRNA_d[k][0][i] for i in range_column)
                colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                outpt.write('\t'.join(colum_NRNA) + '\n')
                outpt.write('\t'.join(colum_SRNA) + '\n')
                outpt.write('\n')


def Tum_VSL(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_Tum_VSL.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            if float(SRNA_d[k][0][15]) >= 27 and float(NRNA_d[k][0][14]) >= 27:
                if float(NDNA_d[k][0][14]) >= 27 and float(SDNA_d[k][0][14]) >= 27:
                    colum_SRNA = list(SRNA_d[k][0][i] for i in range_column)
                    colum_NRNA = list(NRNA_d[k][0][i] for i in range_column)
                    colum_NDNA = list(NDNA_d[k][0][i] for i in range_column)
                    colum_SDNA = list(SDNA_d[k][0][i] for i in range_column)
                    writefile(colum_NDNA, colum_SDNA,
                              colum_NRNA, colum_SRNA, outpt)


def LOH_exome(NDNA_d, SDNA_d):
    with open(base + 'Events_LOH.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SDNA_d.keys()
        for k in sorted(keys):
            chr_sample = SDNA_d[k][0][0]
            pos_sample = SDNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            if float(NDNA_d[k][0][14]) >= 23 and float(SDNA_d[k][0][13]) >= 23:
                colum_NDNA = list(NDNA_d[k][0][i] for i in range_column)
                colum_SDNA = list(SDNA_d[k][0][i] for i in range_column)
                outpt.write('\t'.join(colum_NDNA) + '\n')
                outpt.write('\t'.join(colum_SDNA) + '\n')


def LOH(SRNA_d, NRNA_d, NDNA_d, SDNA_d):
    with open(base + 'Events_LOH.tsv', 'w') as outpt:
        outpt.write('\t'.join(header_up) + '\n')
        keys = SRNA_d.keys()
        for k in sorted(keys):
            chr_sample = SRNA_d[k][0][0]
            pos_sample = SRNA_d[k][0][1]
            key_sample = chr_sample + ":" + pos_sample
            NDNA_SNP_count = int(NDNA_d[k][0][9])
            NDNA_Ref_count = int(NDNA_d[k][0][10])
            SDNA_SNP_count = int(SDNA_d[k][0][9])
            SDNA_Ref_count = int(SDNA_d[k][0][10])
            NRNA_SNP_count = int(NRNA_d[k][0][9])
            NRNA_Ref_count = int(NRNA_d[k][0][10])
            SRNA_SNP_count = int(SRNA_d[k][0][9])
            SRNA_Ref_count = int(SRNA_d[k][0][10])
            # if (NDNA_Ref_count >= 1 and NDNA_SNP_count >= 1 and SDNA_Ref_count == 0 and SDNA_SNP_count>=1):
            # if (NRNA_Ref_count >= 3 or NRNA_SNP_count >= 3 or NRNA_Ref_count == 0 or NRNA_SNP_count == 0 ):
            # if (SRNA_SNP_count>= 1 and SRNA_Ref_count== 0):
            if float(NDNA_d[k][0][14]) >= 23 and float(SDNA_d[k][0][13]) >= 23:
                if (float(NRNA_d[k][0][14]) >= 5 or float(NRNA_d[k][0][13]) >= 5 or float(NRNA_d[k][0][15]) >= 5):
                    if float(SRNA_d[k][0][13]) >= 5:
                        colum_SRNA = list(SRNA_d[k][0][i]
                                          for i in range_column)
                        colum_NRNA = list(NRNA_d[k][0][i]
                                          for i in range_column)
                        colum_NDNA = list(NDNA_d[k][0][i]
                                          for i in range_column)
                        colum_SDNA = list(SDNA_d[k][0][i]
                                          for i in range_column)
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
