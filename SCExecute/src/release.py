#!/bin/env python
RELEASE = "SCExecute"
VERSION = '1.1.0'
PROGRAMS = 'scExecute.py'
INCLUDES = 'common'
if __name__ == '__main__':
    import sys
    print(eval(sys.argv[1]))
