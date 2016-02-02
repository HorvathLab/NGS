#!/bin/env python27
VERSION = '1.0.10'
PROGRAMS = 'readCounts.py LoH.py RNA2DNA.py exonicFilter.py snv_computation.py'
if __name__ == '__main__':
    import sys
    print eval(sys.argv[1])
