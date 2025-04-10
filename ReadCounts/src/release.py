#!/bin/env python
RELEASE = "ReadCounts"
VERSION = '2.5.3'
PROGRAMS = 'readCounts.py phasedReadCounts.py varLoci.py'
INCLUDES = 'common'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
