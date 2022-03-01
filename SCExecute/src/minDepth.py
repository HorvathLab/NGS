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

depth = int(sys.argv[2])
bamfilename = sys.argv[1]

# scan = Intervals(bamfilename,mindepth=depth)
scan = MinDepthBAM(bamfilename,mindepth=depth)
scan.process()
