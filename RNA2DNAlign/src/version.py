#!/bin/env python27
VERSION = '1.0.9'
PROGRAMS = 'readCounts.py LoH.py RNA2DNAlign.py exonicFilter.py snv_computation.py'
if __name__ == '__main__':
    import sys
    print eval(sys.argv[1])
