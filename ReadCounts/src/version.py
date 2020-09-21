#!/bin/env python2.7
VERSION = '2.1.1'
PROGRAMS = 'readCounts.py phasedReadCounts.py'
INCLUDES = 'common'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
