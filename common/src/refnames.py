#!/bin/env python2.7
import sys
from os.path import join, dirname, realpath
try:
    sys.path.append(join(dirname(realpath(__file__)),
                         '..', '..', 'common', 'src'))
except NameError:
    pass

from pysamimport import pysam

for filename in sys.argv[1:]:
    samfile = pysam.Samfile(filename, "rb")
    assert samfile._hasIndex(), "Cannot open BAM index for file %s"%filename
    print filename,"\n  "+", ".join(samfile.references)+"."
