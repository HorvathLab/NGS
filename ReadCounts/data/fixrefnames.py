#!/bin/env python3.12

import sys
import pysam

mapping = dict(map(str.split,map(str.strip,open(sys.argv[1]))))

for f in sys.argv[2:]:
  base,extn = f.rsplit('.')
  nf = base+".out."+extn
  if extn == 'bam':

    bam_file = pysam.AlignmentFile(f, "rb")
    header = bam_file.header.to_dict()

    current_names = [r["SN"] for r in header["SQ"]]

    for i, ref_name in enumerate(current_names):
        if ref_name in mapping:
            header["SQ"][i]["SN"] = mapping[ref_name]

    output_bam = pysam.AlignmentFile(nf, "wb", header=header)

    for read in bam_file.fetch():
        output_bam.write(read)

    bam_file.close()
    output_bam.close()

  else:

    out = open(nf,'w')
    for l in open(f):
        if l[0] in ("@","#"):
            out.write(l)
            continue
        sl = l.split(None,1)
        if len(sl) < 1:
            out.write(l)
            continue
        slzerolen = len(sl[0])
        sl[0] = mapping[sl[0]]
        if l[slzerolen] == "\t":
            nl = "\t".join(sl)
        else:
            nl = " ".join(sl)
        out.write(nl)
    out.close()
