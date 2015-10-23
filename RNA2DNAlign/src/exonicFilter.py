#!/bin/env python27
import os
import csv
import sys
from operator import itemgetter
import operator
d = {}
header_dict = {}
header1 = []
header2 = []
header3 = []
header4 = []

from version import VERSION
VERSION = '1.0.0 (%s)' % (VERSION,)

from optparse_gui import OptionParser
parser = OptionParser(version=VERSION)

parser.add_option("--exons", type="file", dest="exons", default=None,
                  help="Exonic coordinates (sored). Required.", notNone=True,
                  filetypes=[("Exonic Coordinates", "*.txt")])
parser.add_option("--input", type="file", dest="input", default=None,
                  help="Input SNPs in VCF format. Required",
                  filetypes=[("Input SNP File", "*.vcf")])
parser.add_option("--output", type="savefile", dest="output", default=None,
                  help="Output file. Required",
                  filetypes=[("Output SNP File", "*.vcf")])

opt, args = parser.parse_args()


def ReadVCF(file_name):
    d[file_name] = []
    header_dict = []
    with open((file_name), 'r') as f:
        reader = csv.reader((f), delimiter='\t')
        for row in reader:
            if row[0].startswith("##"):
                if "NDNA" in file_name:
                    header1.append(row[0])
                    continue
                if "SDNA" in file_name:
                    header2.append(row[0])
                    continue
                if "NRNA" in file_name:
                    header3.append(row[0])
                    continue
                if "SRNA" in file_name:
                    header4.append(row[0])
                    continue
            if '#' not in row[0]:
                if row[0] != 'MT':
                    d[file_name].append(row)

exoncoords = opt.exons
folder_input = opt.input
outfile = opt.output
l_pathname = []


file = folder_input
print "the file", file
if file.endswith('.vcf') or file.endswith('.txt'):
    print
    print "Now parsing file", file,
    print "Please wait, thanks...."
    ReadVCF(file)
    print
    print "done parsing for "
    print file


def chrorder(chr):
    if chr != "X":
        if chr != "Y":
            return int(chr[0:])
    if chr == "X":
        return 23
    if chr == "Y":
        return 24

#=========================================================================
# A function called "all_filteration" which does take one argument; that is, a list of all the variants present
# within the vcf file. This function does the filteration process on the vcf file based the following aspect:
 # The variants chromosomal positions from vcf file have to be within an exonic regions
#=========================================================================
headers = ['#CHROM', 'POS',  'ID', 'REF',   'ALT',
           'QUAL',  'FILTER', 'INFO', 'FORMAT', 'MCF7']


def all_filteration(d):

    # Opening, Reading the Exonic coordinates from ensemble database, and
    # opening the output file to write to
    with open(exoncoords, 'r') as csvfile:
        #  f_list = f.split('/')
     # file = f_list[1]
        #dir = f_list[0]
        # print "DIR", dir
        # print file
        file = f.split('/')[-1]
        with open(outfile, 'w') as outpt:
            if 'NDNA' in file:
                outpt.write("\n".join(header1) + '\n')
            if 'SDNA' in file:
                outpt.write("\n".join(header2) + '\n')
            if 'NRNA' in file:
                outpt.write("\n".join(header3) + '\n')
            if 'SRNA' in file:
                outpt.write("\n".join(header4) + '\n')
            outpt.write("\t".join(headers) + '\n')
            # Initiating a counter to looping around the list of variants (vcf
            # file)
            variant_count = 0
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:   # This loop belongs to the exonic coordinates file
                if row[0] != 'X'and row[0] != 'Y':
                    chrom_exonic_coord = int(row[0])
                else:
                    chrom_exonic_coord = (row[0])
                exonic_start_pos = int(row[1]) + 1
                exonic_end_pos = int(row[2])
                # This is to check that counter doesn't outrun the end of the
                # vcf list.
                if (variant_count < len(d)):
                    chrom_variant = d[variant_count][0]
                    # The following is to parse the chromosome which is of type string with two steps process:
                    # If it is numerical then convert the string to int , otherwise leave it as string. This will enable to do
                    # numerical comparision instead of string one.
                    if 'MT' not in d[variant_count][0]:
                        if chrom_variant != 'X' and chrom_variant != 'Y':
                            chrom_variant = int(d[variant_count][0])
                        else:
                            chrom_variant = d[variant_count][0]
                    variant_pos = int(d[variant_count][1])

                    reads_var_unsplit = d[variant_count][7]
                    # We want to loop and insure that the chromosome of both exonic coords and variants are equal or the chromosome of exonic coord is greater than the variants'
                    # So we can get the next chromosome from the exonic coordinates
                    # The only addition I added is to say that do the while loop when exonic_end_pos is greater than variant_pos to insure it is within the exonic
                    # so that it will stop the loop once the vairant pos is greater than exonic_end_pos and take the next exonic coords. Also, we want to start
                    # what we left off and currntly the program is not doing this. Because we Initiating the variant_count to 0 so every time it takes another exonic region
                    # the counter is set to 0 which means that it will start again from the first position in the vcf file which we don't want, we want the next variant not the very
                    # first one.
                    while ((chrom_exonic_coord == chrom_variant or chrom_exonic_coord > chrom_variant) and (exonic_end_pos > variant_pos)):
                        # Below is just to extract the two digits (forward and reverse variant reads) from "DP4" of the "Info" column in VCF file
                        # print"ya kareem", chrom_variant, " ", variant_pos
                        variant_count += 1

                        split_rd = reads_var_unsplit.split(';')

                        for dp4 in split_rd:
                            # print "DP4", dp4

                            if 'DP4=' in dp4:

                                split_DP4_Colmn = dp4[4:].split(',')
                        if (chrom_exonic_coord == chrom_variant):

                            # if (int(split_DP4_Colmn[2]) >=1 and
                            # int(split_DP4_Colmn[3]) >=1):

                            if (exonic_start_pos <= variant_pos <= exonic_end_pos):

                                # writing to the file the passed (filtered) variants
                                # print "FILTER"
                                # chrom_exonic_coord == chrom_variant

                                # print chrom_exonic_coord  ,":", chrom_variant
                                # print "Check in One:", exonic_start_pos, " ",
                                # variant_pos, " ", exonic_end_pos
                                outpt.writelines(
                                    "\t".join(d[variant_count - 1]) + '\n')

                            # print "Not in ", chrom_exonic_coord  ,":", chrom_variant
                            # print "NOT IN:", exonic_start_pos, " ", variant_pos, " ", exonic_end_pos
                               # print " Done FILTER"
                        # When we are done reading the variants list , then
                        # exit
                        if (variant_count >= len(d)):
                            break
                        if 'MT' not in d[variant_count][0]:
                            if d[variant_count][0] != 'X' and d[variant_count][0] != 'Y':
                                chrom_variant = int(d[variant_count][0])
                            else:
                                chrom_variant = d[variant_count][0]
                        variant_pos = int(d[variant_count][1])

                        reads_var_unsplit = d[variant_count][7]
                    while (variant_pos > exonic_end_pos and chrom_variant < chrom_exonic_coord):
                            # print variant_count, len(d)
                        if (variant_count == len(d) - 1):
                            break
                        variant_count += 1
                        # print  d[variant_count][0], d[variant_count][1],
                        # variant_count
                        if d[variant_count][0] != 'X':
                            if d[variant_count][0] != 'Y':
                                chrom_variant = int(d[variant_count][0])
                        # if d[variant_count][0] != 'X' or d[variant_count][0]
                        # != 'X':
                        else:
                            chrom_variant = d[variant_count][0]
                        variant_pos = int(d[variant_count][1])
                        # print"hi ya rab", chrom_variant, " ", variant_pos,
                        # "and", chrom_exonic_coord," ", exonic_end_pos

                        if (chrom_exonic_coord == chrom_variant):

                          # if (int(split_DP4_Colmn[2]) >=1 and
                          # int(split_DP4_Colmn[3]) >=1):

                            if (exonic_start_pos <= variant_pos <= exonic_end_pos):

                               #     print "hi"
                                # writing to the file the passed (filtered) variants
                                # print "Check in:", exonic_start_pos, " ",
                                # variant_pos, " ", exonic_end_pos
                                outpt.writelines(
                                    "\t".join(d[variant_count - 1]) + '\n')

#sorted_d = sorted(d.items(), key=operator.itemgetter(1))
# for i in range(len(sorted_d[1][0])):
# print "i",sorted_d
print
print "Now Filtering Data"
for f in d:
            # print
    print "the F", f
    # for path in l_pathname:
    #  print "path", path
    # print "check", d[f]

    # for i in d[f]:
    #	print i

    # print "CHECK", f
    #  print "now we sort"
   # print
    sort_l = sorted(d[f], key=lambda l: chrorder(l[0]))
    d[f] = sort_l
    # for i in d[f]:
    #	print i
    # print
    # print "f", f
    all_filteration(d[f])
    # print all_filteration(d[f])
    print "done filtering for"
    print f
print
print "Done All with Parsing and Filtering"
print
