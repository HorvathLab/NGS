#!/bin/env python
RELEASE = "ReadCounts"
VERSION = '2.3.1'
PROGRAMS = 'readCounts.py phasedReadCounts.py'
INCLUDES = 'common'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
