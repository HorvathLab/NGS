#!/bin/env python2.7
VERSION = '2.0.0'
PROGRAMS = 'readCounts.py phasedReadCounts.py'
INCLUDES = 'common'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
