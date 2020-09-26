#!/bin/env python
VERSION = '2.2.0'
PROGRAMS = 'readCounts.py phasedReadCounts.py readCountsMatrix.py scReadCounts.py'
INCLUDES = 'common'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
