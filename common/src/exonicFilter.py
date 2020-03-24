#!/bin/env python2.7
import os
import csv
import sys
from operator import itemgetter
import operator

from os.path import join, dirname, realpath, split
try:
    scriptdir = dirname(realpath(__file__))
except NameError:
    scriptdir = dirname(realpath(sys.argv[0]))
sys.path.append(join(scriptdir, '..', '..', 'common', 'src'))

from version import VERSION
VERSION = '1.0.4 (%s)' % (VERSION,)

from optparse_gui import OptionParser
parser = OptionParser(version=VERSION)

parser.add_option("--exons", type="file", dest="exons", default=None,
                  help="Exonic coordinates (sored). Required.", notNone=True,
                  filetypes=[("Exonic Coordinates", "*.txt")])
parser.add_option("--input", type="file", dest="input", default=None,
                  help="Input SNVs. Required",
                  filetypes=[("Input SNV File", "*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("--output", type="savefile", dest="output", default=None,
                  help="Output file. Required",
                  filetypes=[("Output SNV File", "*.vcf;*.tsv")])

opt, args = parser.parse_args()

from dataset import XLSFileTable, CSVFileTable, TSVFileTable, XLSXFileTable, TXTFileTable

def ReadVCF(file_name):
    fileheader = []
    rows = []
    chrom = set()
    
    with open((file_name), 'r') as f:
        reader = csv.reader((f), delimiter='\t')
        inheader = True
        for row in reader:
            if row[0].startswith("#") and inheader:
                fileheader.append('\t'.join(row))
            else:
                inheader = False
                chrom.add(row[0])
                rows.append(row)
    return fileheader, chrom, rows

def ReadTSV(filename):
    snvheaders = [_f for _f in """CHROM POS REF ALT""".split() if _f]  
    base, extn = filename.rsplit('.', 1)
    extn = extn.lower()
    if extn == 'csv':
        snvs = CSVFileTable(filename=filename)
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

    assert(snvs.headers()[:4] == snvheaders)

    chrom = set()
    snvdata = []
    for r in snvs:
        ri = list(map(r.get,snvs.headers()))
        chrom.add(ri[0])
        snvdata.append(ri)

    return ["\t".join(snvs.headers())], chrom, snvdata
                                                                                                
exoncoords = opt.exons
filename = opt.input
outfile = opt.output

if outfile.endswith('.vcf'):
    assert(filename.endswith('.vcf'))
if not filename.endswith('.vcf'):
    assert(outfile.endswith('.tsv'))

print("Reading", filename, "...", end=' ')
sys.stdout.flush()

if filename.rsplit('.',1)[-1].lower() == 'vcf':
    fileheader,chrlab,snvdata = ReadVCF(filename)
else:
    fileheader,chrlab,snvdata = ReadTSV(filename)

from chromreg import ChromLabelRegistry
chrreg = ChromLabelRegistry()
chrreg.add_labels(filename,chrlab)

for i in range(len(snvdata)-1,-1,-1):
    chrlab = snvdata[i][0]
    chrom = chrreg.label2chrom(filename,chrlab)
    if not chrreg.isnumberedchrom(chrom) and \
       not chrreg.issexchrom(chrom):
        del snvdata[i]
        continue
    snvdata[i][0] = chrom
    snvdata[i][1] = int(snvdata[i][1])

# Expected chromosome labels from UCSC file...
# Extra numeric chromosomes might make this more robust for organisms
# with more chromosomes? Doesn't hurt normal (human) case either way.
# Should we add mitochondria here too? What notation does UCSC use?
exonlabels = list(map(str,list(range(1,100)))) + ["X","Y","MT"]
chrreg.add_labels(exoncoords,exonlabels)

chrreg.default_chrom_order()
chrorder = chrreg.chrom_order

snvdata.sort(key=lambda sd: (chrorder(sd[0]),sd[1]))

print("done")

#=========================================================================
# A function called "all_filteration" which does take one argument;
# that is, a list of all the variants present within the vcf
# file. This function does the filteration process on the vcf file
# based the following aspect: The variants chromosomal positions from
# vcf file have to be within an exonic regions
#=========================================================================

def all_filteration(d):

    #
    # This code assumes that the exon file is sorted by chromosome and
    # start and end position. Chromosome sort order (for human) is
    #     1, 2, ..., 9, 10, 11, ...., 19, 20, 21, 22, X, Y.
    # Position sort order is as integers.
    #

    # Opening, Reading the Exonic coordinates from EMBL database, and
    # opening the output file to write to
    with open(exoncoords, 'r') as csvfile:
        with open(outfile, 'w') as outpt:
            outpt.write("\n".join(fileheader) + '\n')

            variant_count = 0
            last_exonic_coords = (-1e+20,-1e+20,-1e+20)
            reader = csv.reader(csvfile, delimiter='\t')
            
            for row in reader:   # This loop belongs to the exonic coordinates file

                chrom_exonic_coord = chrreg.label2chrom(exoncoords,row[0])
                assert chrom_exonic_coord != None, "Unexpected chromosome label in exon coordinates file"
                exonic_start_pos = int(row[1]) + 1
                exonic_end_pos = int(row[2])

                assert (chrorder(chrom_exonic_coord),exonic_start_pos,exonic_end_pos) >= last_exonic_coords, \
                        "Exon coordinates file is not correctly ordered by chromosome and start/end positions"

                if (variant_count < len(d)):
                    
                    chrom_variant = d[variant_count][0]
                    variant_pos = d[variant_count][1]
                    
                    while ((chrom_exonic_coord == chrom_variant) and (exonic_end_pos > variant_pos)) or (chrorder(chrom_exonic_coord) > chrorder(chrom_variant)):
        
                        if (chrom_exonic_coord == chrom_variant):

                            if (exonic_start_pos <= variant_pos <= exonic_end_pos):

                                chrlab = chrreg.chrom2label(filename,d[variant_count][0])
                                chrpos = str(d[variant_count][1])
                                outpt.writelines("\t".join([chrlab,chrpos] + d[variant_count][2:]) + '\n')

                        variant_count += 1

                        if (variant_count >= len(d)):
                            break

                        chrom_variant = d[variant_count][0]
                        variant_pos = d[variant_count][1]

                last_exonic_coords = (chrorder(chrom_exonic_coord),exonic_start_pos,exonic_end_pos)

print("Filtering", filename, "...", end=' ')
sys.stdout.flush()

all_filteration(snvdata)

print("done")
