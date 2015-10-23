#!/bin/env python27
import tempfile
import os
import summary_analysis
from subprocess import Popen, PIPE, STDOUT, check_output, CalledProcessError
import sys
import os
import os.path
import glob
import copy
import traceback
import time
import re
import csv
import tempfile
import urllib
import shutil
import atexit
import subprocess
import time
import math
from collections import defaultdict, Counter
from os.path import join, dirname, realpath, split
try:
    scriptdir = dirname(realpath(__file__))
except NameError:
    scriptdir = dirname(realpath(sys.argv[0]))
sys.path.append(join(scriptdir, '..', '..', 'common', 'src'))
try:
    scriptextn = "." + os.path.split(sys.argv[0])[1].rsplit('.', 1)[1]
except:
    scriptextn = ""
from optparse_gui import OptionParser, OptionGroup, GUI, UserCancelledError, ProgressText
from util import *
from fisher import *
from operator import itemgetter

from version import VERSION
VERSION = '1.0.0 (%s)' % (VERSION,)


def excepthook(etype, value, tb):
    traceback.print_exception(etype, value, tb)
    print >>sys.stderr, "Type <Enter> to Exit...",
    sys.stderr.flush()
    raw_input()

toremove = []


def cleanup():
    for d in toremove:
        shutil.rmtree(d, ignore_errors=True)
atexit.register(cleanup)

if not GUI() and len(sys.argv) == 2 and sys.argv[1] == '--GUI':
    from optparse_gui.needswx import *
    sys.exit(1)

if GUI() and len(sys.argv) == 1:
    from optparse_gui import OptionParserGUI
    parser = OptionParserGUI(version=VERSION)
    error_kwargs = {'exit': False}
    sys.excepthook = excepthook
else:
    parser = OptionParser(version=VERSION)
    error_kwargs = {}

exfilt = OptionGroup(parser, "SNP Filtering")
readcounts = OptionGroup(parser, "Read Counting")
# advanced = OptionGroup(parser, "Advanced")
snpannot = OptionGroup(parser, "SNP Annotation")
parser.add_option("-s", "--snps", type="files", dest="snps", default=None,
                  help="Single-Nucleotide-Polymophisms. Required.", name="SNPs",
                  notNone=True, remember=True,
                  filetypes=[("SNPs", "*.vcf")])
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignments in BAM/SAM format. Required.", name="Read Alignments",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignments (BAM/SAM Format)", "*.bam;*.sam")])
exfilt.add_option("-e", "--exoncoords", type="file", dest="exoncoords", default=None,
                  help="Exon coordinates for SNP filtering. Optional.", name="Exon Coords.",
                  remember=True,
                  filetypes=[("Exonic Coordinates", "*.txt")])
snpannot.add_option("-d", "--darned", type="file", dest="darned", default="",
                    help="DARNED Annotations. Optional.", remember=True,
                    filetypes=[("DARNED Annotations", "*.txt")])
snpannot.add_option("-c", "--cosmic", type="file", dest="cosmic", default="",
                    help="COSMIC Annotations. Optional.", remember=True,
                    filetypes=[("COSMIC Annotations", "*.tsv;*.tsv.gz")])
readcounts.add_option("-m", "--minreads", type="int", dest="minreads", default=10, remember=True,
                      help="Minimum number of good reads at SNP locus per alignment file. Default=10.", name="Min. Reads")
readcounts.add_option("-f", "--alignmentfilter", action="store_false", dest="filter", default=True, remember=True,
                      help="(Turn off) alignment filtering by length, edits, etc.", name="Filter Alignments")
readcounts.add_option("-U", "--uniquereads", action="store_true", dest="unique", default=False, remember=True,
                      help="Consider only distinct reads.", name="Unique Reads")
readcounts.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, remember=True,
                      help="Quiet.", name="Quiet")

parser.add_option("-o", "--output", type="savedir", dest="output", remember=True,
                  help="Output folder. Required.", default=None, notNone=True,
                  name="Output Folder")

parser.add_option_group(exfilt)
parser.add_option_group(readcounts)
parser.add_option_group(snpannot)

opt = None
while True:
    if 'exit' in error_kwargs:
        try:
            opt, args = parser.parse_args(opts=opt)
        except UserCancelledError:
            sys.exit(0)
    else:
        opt, args = parser.parse_args()

    break


def makedirs(d):
    if os.path.isdir(d):
        return
    os.makedirs(d)


def execprog(prog, *args, **kw):
    progpath = os.path.join(scriptdir, prog + scriptextn)
    assert os.path.exists(progpath), "%s does not exist" % (progpath,)
    if kw.get('verbose', True):
        argstr = " ".join(
            map(lambda a: a if " " not in a else '"%s"' % a, args))
        print >>sys.stderr, "Executing:\n  %s %s" % (prog + scriptextn, argstr)
    if progpath.endswith('.py'):
        sys.argv = [progpath] + list(args)
        execfile(progpath, globals())
    else:
        status = subprocess.call([progpath] + list(args))
        assert(status == 0)
    return True

mainopt = copy.deepcopy(opt)

# Apply exonic filter on SNPs if desired...
snpfiles = []
for snpfile in mainopt.snps:
    if mainopt.exoncoords:
        base, extn = snpfile.rsplit('.', 1)
        basedir, basename = split(base)
        outfile = join(mainopt.output, basename + '.filtered.' + extn)
        if not os.path.exists(outfile):
            makedirs(mainopt.output)
            execprog("exonicFilter", "--exons", mainopt.exoncoords,
                     "--input", snpfile, "--output", outfile)
        snpfiles.append(outfile)
    else:
        snpfiles.append(snpfile)

# Apply readCounts to SNPs and aligned reads. Pass on mainoptions as needed...
outfile = join(mainopt.output, "readCounts.tsv")
if not os.path.exists(outfile):

    args = ["-F",
            "-r", " ".join(mainopt.alignments),
            "-s", " ".join(snpfiles),
            "-o", outfile]
    args.extend(["-m", str(mainopt.minreads)])
    if not mainopt.filter:
        args.append("-f")
    if mainopt.unique:
        args.append("-U")
    if mainopt.quiet:
        args.append("-q")

    makedirs(mainopt.output)
    execprog("readCounts", *args)

# Set up and apply snp_computation.py
args = ["--counts", outfile]
if mainopt.darned:
    args.extend(["--darned", mainopt.darned])
if mainopt.cosmic:
    args.extend(["--cosmic", mainopt.cosmic])
execprog("snp_computation", *args)

# Summarize events
if os.path.exists(join(mainopt.output, "summary_result.txt")):
    os.unlink(join(mainopt.output, "summary_result.txt"))
for f in glob.glob(join(mainopt.output, "Events*.tsv")):
    summary_analysis.read_events(f)
