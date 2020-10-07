#!/bin/env python
VERSION = '1.0.0'
PROGRAMS = 'readCountsMatrix.py scReadCounts.py'
INCLUDES = 'common ReadCounts'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
