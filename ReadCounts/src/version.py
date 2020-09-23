#!/bin/env python2.7
VERSION = '2.1.3'
PROGRAMS = 'readCounts.py phasedReadCounts.py readCountsMatrix.py'
INCLUDES = 'common'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
