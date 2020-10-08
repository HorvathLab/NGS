#!/bin/env python
RELEASE = "SCReadCounts"
VERSION = '1.0.0'
PROGRAMS = 'readCountsMatrix.py scReadCounts.py readCounts.py'
INCLUDES = 'common ReadCounts'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
