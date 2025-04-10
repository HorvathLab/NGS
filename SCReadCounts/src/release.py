#!/bin/env python
RELEASE = "SCReadCounts"
VERSION = '1.4.1'
PROGRAMS = 'readCountsMatrix.py scReadCounts.py readCounts.py varLoci.py scVarLoci.py rcio.py scCountLoci.py exonicFilter.py'
INCLUDES = 'common ReadCounts'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
