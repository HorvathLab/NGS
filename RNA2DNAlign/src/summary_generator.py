#!/bin/env python27
from operator import itemgetter
from collections import defaultdict, Counter
import os,sys
import csv
def read_events(file):
        cosmic_dic = defaultdict(list)
        count = 0
        count_sam = 0
        l = []
        total = Counter()
        l_sam = []
        count_A_T = 0; count_A_G = 0;count_A_C = 0
        count_C_T = 0; count_C_G = 0;count_C_A= 0
        count_G_T = 0; count_G_C = 0;count_G_A= 0
        count_T_A = 0; count_T_G = 0;count_T_C = 0

        with open(file ,'Ur') as f:
          reader=csv.reader((f),delimiter='\t')
          for row in reader:
           if row:
             if 'AlignedReads' not in row[0]:
                if (row[1], row[2]) not in l:
                        if row[3] == 'A' and row[4] == 'T':
                                count_A_T = count_A_T+1
                        elif  row[3] == 'A' and row[4] == 'G':
                                count_A_G = count_A_G+1
                        elif row[3] == 'A' and row[4] == 'C':
                                count_A_C = count_A_C +1
                        elif row[3] == 'C' and row[4] == 'T':
                                count_C_T = count_C_T +1
                        elif  row[3] == 'C' and row[4] == 'G':
                                count_C_G = count_C_G+1
                        elif  row[3] == 'C' and row[4] == 'A':
                                count_C_A = count_C_A+1
                        elif row[3] == 'G' and row[4] == 'A':
                                count_G_A = count_G_A + 1
                        elif row[3] == 'G' and row[4] == 'C':
                                count_G_C = count_G_C + 1
                        elif row[3] == 'G' and row[4] == 'T':
                                count_G_T = count_G_T + 1
                        elif row[3] == 'T' and row[4] == 'A':
                                         count_T_A = count_T_A + 1
                        elif row[3] == 'T' and row[4] == 'G':
                                count_T_G = count_T_G + 1
                        elif row[3] == 'T' and row[4] == 'C':
                                count_T_C = count_T_C + 1
      	                l.append((row[1], row[2]))
                        count = count + 1
                        total[row[1]] += 1
                if row[0] not in l_sam:
                        count_sam = count_sam + 1
                        l_sam.append(row[0])
        with open("summary_result.txt", 'a') as out:
                                out.write("##Summary Result of"+" "+ file)
                                out.write("\n")
                                out.write("#Number of Sample Analyzed:"+ " "+str(count_sam))
                                out.write("\n")
                                out.write("#Number of Total SNPs:"+ " "+ str(count))
                                out.write("\n")
                                out.write("#Most Frequent Chromosome:"+ " chr"+ " " +str(max(total.iteritems(), key=itemgetter(1))[0]))
                                out.write("\n")
                                out.write("#Least Frequent Chromosome:"+ " chr"+ " " +str(min(total.iteritems(), key=itemgetter(1))[0]))
                                out.write("\n")
                                out.write("#Change A > T: " + str(count_A_T) + "\t" + "#Change A > G: " + str(count_A_G) + "\t" + "#Change A > C: " + str(count_A_C))
					
				out.write("\n")
                                out.write("#Change C > T: " + str(count_C_T) + "\t" + "#Change C > G: " + str(count_C_G) + "\t" + "#Change C > A: " + str(count_C_A))
                                out.write("\n")
                                out.write("#Change G > T: "+ str(count_C_T)+ "\t" + "#Change G > C: " + str(count_G_C) + "\t" + "#Change G > A: " + str(count_G_A))
                                out.write("\n")
                                out.write("#Change T > G: "+ str(count_T_G) + "\t" + "#Change T > C: " + str(count_T_C) + "\t" + "#Change T > A: " + str(count_T_A))

				out.write("\n")
				out.write("\n")


							
for subdir, dirs, files in os.walk('./'):
    for file in files:
        if file.startswith('Events'):
            read_events(file)
