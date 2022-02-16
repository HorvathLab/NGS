#!/bin/env python
RELEASE = "SCReadCounts"
VERSION = '1.1.10'
PROGRAMS = 'readCountsMatrix.py scReadCounts.py readCounts.py'
INCLUDES = 'common ReadCounts'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
