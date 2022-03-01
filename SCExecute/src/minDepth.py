#!/bin/env python3
import sys
import os
import os.path
import traceback
from os.path import join, dirname, realpath
try:
    sys.path.append(join(dirname(realpath(__file__)),
                         '..', '..', 'common', 'src'))
except NameError:
    pass

from intervalscan import MinDepthBAM, Intervals

outputbam = False
while sys.argv[1].startswith('-'):
    if sys.argv[1] == "-b":
        outputbam  = True
        sys.argv.pop(1)
    else:
        raise RuntimeError("Bad argument")
    
bamfilename = sys.argv[1]
depth = 1
if len(sys.argv) >= 3:
    depth = int(sys.argv[2])

if not outputbam:
  scan = Intervals(bamfilename,mindepth=depth)
else:
  scan = MinDepthBAM(bamfilename,mindepth=depth)
scan.process()
