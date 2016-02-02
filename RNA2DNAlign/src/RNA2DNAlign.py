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
VERSION = '1.0.4 (%s)' % (VERSION,)


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
						
exfilt = OptionGroup(parser, "Filtering")
readcounts = OptionGroup(parser, "Read Counting")
# advanced = OptionGroup(parser, "Advanced")
regexs = OptionGroup(parser, "Filename Matching")
snvannot = OptionGroup(parser, "SNV Annotation")
parser.add_option("-s", "--snvs", type="files", dest="snvs", default=None,
                  help="Single-Nucleotide-Variant files. Required.", name="SNV Files",
                  notNone=True, remember=True,
                  filetypes=[("SNV Files", "*.vcf")])
parser.add_option("-r", "--readalignments", type="files", dest="alignments", default=None,
                  help="Read alignment files in indexed BAM format. Required.", name="Read Alignment Files",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignment Files (Indexed BAM)", "*.bam")])
exfilt.add_option("-e", "--exoncoords", type="file", dest="exoncoords", default=None,
                  help="Exon coordinates for SNV filtering. Optional.", name="Exon Coords.",
                  remember=True,
                  filetypes=[("Exonic Coordinates", "*.txt")])
regexs.add_option("--normalexomere", type="str", dest="normalexomere", default='GDNA',
                  help="Germline/Normal exome filename regular expression. Default: GDNA.",
                  remember=True, name="Germline Exome RE")
regexs.add_option("--normaltransre", type="str", dest="normaltransre", default='NRNA',
                  help="Normal transcriptome filename regular expression. Default: NRNA.",
                  remember=True, name="Normal Transcr. RE")
regexs.add_option("--tumorexomere", type="str", dest="tumorexomere", default='SDNA',
                  help="Somatic/Tumor exome filename regular expression. Default: SDNA.",
                  remember=True, name="Somatic Exome RE")
regexs.add_option("--tumortransre", type="str", dest="tumortransre", default='TRNA',
                  help="Tumor transcriptome filename regular expression. Default: TRNA.",
                  remember=True, name="Tumor Transcr. RE")
snvannot.add_option("-d", "--darned", type="file", dest="darned", default="",
                    help="DARNED Annotations. Optional.", remember=True,
                    filetypes=[("DARNED Annotations", "*.txt")])
snvannot.add_option("-c", "--cosmic", type="file", dest="cosmic", default="",
                    help="COSMIC Annotations. Optional.", remember=True,
                    filetypes=[("COSMIC Annotations", "*.tsv;*.tsv.gz")])
readcounts.add_option("-m", "--minreads", type="int", dest="minreads", default=10, remember=True,
                      help="Minimum number of good reads at SNV locus per alignment file. Default=10.", name="Min. Reads")
readcounts.add_option("-f", "--alignmentfilter", action="store_false", dest="filter", default=True, remember=True,
                      help="(Turn off) alignment filtering by length, edits, etc.", name="Filter Alignments")
readcounts.add_option("-U", "--uniquereads", action="store_true", dest="unique", default=False, remember=True,
                      help="Consider only distinct reads.", name="Unique Reads")
readcounts.add_option("-t", "--threadsperbam", type="int", dest="tpb", default=1, remember=True,
                    help="Worker threads per alignment file. Indicate no threading with 0. Default=1.",
                      name="Threads/BAM")
readcounts.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, remember=True,
                      help="Quiet.", name="Quiet")

parser.add_option("-o", "--output", type="savedir", dest="output", remember=True,
                  help="Output folder. Required.", default=None, notNone=True,
                  name="Output Folder")

parser.add_option_group(exfilt)
parser.add_option_group(readcounts)
parser.add_option_group(regexs)                 
parser.add_option_group(snvannot)

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
        execfile(progpath,{})
    else:
        status = subprocess.call([progpath] + list(args))
        assert(status == 0)
    return True

# Apply exonic filter on SNVs if desired...
snvfiles = []
for snvfile in opt.snvs:
    if opt.exoncoords:
        base, extn = snvfile.rsplit('.', 1)
        if extn != 'vcf':
            extn = 'tsv'
        basedir, basename = split(base)
        outfile = join(opt.output, "." + basename + '.filtered.' + extn)
        if not os.path.exists(outfile):
            makedirs(opt.output)
            execprog("exonicFilter", "--exons", opt.exoncoords,
                     "--input", snvfile, "--output", outfile)
        snvfiles.append(outfile)
    else:
        snvfiles.append(snvfile)

# Apply readCounts to SNVs and aligned reads. Pass on options as needed...
outfile = join(opt.output, "readCounts.tsv")
if not os.path.exists(outfile):

    args = ["-F",
            "-r", " ".join(opt.alignments),
            "-s", " ".join(snvfiles),
            "-o", outfile]
    args.extend(["-m", str(opt.minreads)])
    args.extend(["-t", str(opt.tpb)])
    if not opt.filter:
        args.append("-f")
    if opt.unique:
        args.append("-U")
    if opt.quiet:
        args.append("-q")

    makedirs(opt.output)
    execprog("readCounts", *args)

# Set up and apply snv_computation.py
args = ["--counts", outfile]
if opt.darned:
    args.extend(["--darned", opt.darned])
if opt.cosmic:
    args.extend(["--cosmic", opt.cosmic])
args.extend(["--normalexomere",opt.normalexomere])
args.extend(["--normaltransre",opt.normaltransre])
args.extend(["--tumorexomere",opt.tumorexomere])
args.extend(["--tumortransre",opt.tumortransre])

execprog("snv_computation", *args)

# Summarize events
if os.path.exists(join(opt.output, "summary_result.txt")):
    os.unlink(join(opt.output, "summary_result.txt"))
for event in "RNAed T-RNAed VSE T-VSE VSL T-VSL SOM LOH".split():
    f = join(opt.output, "Events_%s.tsv"%(event,))
    if os.path.exists(f):
        summary_analysis.read_events(f)
