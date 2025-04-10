#!/bin/env python3.12

import sys
import pysam

bam_file = pysam.AlignmentFile(sys.argv[1], "rb")
header = bam_file.header.to_dict()
current_names = [r["SN"] for r in header["SQ"]]

for i, ref_name in enumerate(current_names):
    print("%s\t%02d"%(ref_name,i+1,))
