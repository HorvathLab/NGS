#!/bin/env python2.7
VERSION = '2.1.2'
PROGRAMS = 'readCounts.py LoH.py RNA2DNAlign.py exonicFilter.py snv_computation.py'
INCLUDES = 'common ReadCounts'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
