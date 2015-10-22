#!/bin/env python27
import tempfile,os, summary_analysis
from subprocess import Popen, PIPE, STDOUT, check_output, CalledProcessError
from tempfile import TemporaryFile
import sys, os, os.path, glob, copy, traceback, time, re, csv, tempfile, urllib, shutil, atexit, subprocess, time, math
from collections import defaultdict, Counter
from os.path import join, dirname, realpath
try:
    sys.path.append(join(dirname(realpath(__file__)),'..','..','common','src'))
except NameError:
    pass
from optparse_gui import OptionParser, OptionGroup, GUI, UserCancelledError, ProgressText
from util import *
from fisher import *
from operator import itemgetter

VERSION='1.0.0'

def excepthook(etype,value,tb):
    traceback.print_exception(etype,value,tb)
    print >>sys.stderr, "Type <Enter> to Exit...",
    sys.stderr.flush()
    raw_input()
sys.excepthook = excepthook

toremove = []
def cleanup():
    for d in toremove:
        shutil.rmtree(d,ignore_errors=True)
atexit.register(cleanup)

if not GUI() and len(sys.argv) == 2 and sys.argv[1] == '--GUI':
    from optparse_gui.needswx import *
    sys.exit(1)

if GUI() and len(sys.argv) == 1:
    from optparse_gui import OptionParserGUI
    parser = OptionParserGUI(version=VERSION)
    error_kwargs = {'exit': False}
else:
    parser = OptionParser(version=VERSION)
    error_kwargs = {}
advanced = OptionGroup(parser, "Advanced")
parser.add_option("-s","--snps",type="files",dest="snps",default=None,
                  help="Single-Nucleotide-Polymophisms. Required.", name="SNPs",
                  notNone=True, remember=True,
                  filetypes=[("SNPs","*.vcf;*.csv;*.tsv;*.xls;*.xlsx;*.txt")])
parser.add_option("-r","--readalignments",type="files",dest="alignments",default=None,
                  help="Read alignments in BAM/SAM format. Required.", name="Read Alignments",
                  notNone=True, remember=True,
                  filetypes=[("Read Alignments (BAM/SAM Format)","*.bam;*.sam")])
					
advanced.add_option("-f","--filtrsnps",action="store_true",dest="filtr",default=False,remember=True,
                  help="Exonic Position snps. Default=False.",name="Exonic snps")
advanced.add_option("-m","--minreads",action="store_true",dest="minreads",default=10,remember=True,
                  help="Minimum number of good reads at SNP locus per alignment file. Default=10.",name="Min. Reads")
advanced.add_option("-F","--full",action="store_true",dest="full",default=False,remember=True,
                  help="Output extra diagnostic read count fields. Default=False.",name="All Fields")
advanced.add_option("-U","--uniquereads",action="store_true",dest="unique",default=False,remember=True,
                  help="Consider only distinct reads.",name="Unique Reads")
advanced.add_option("-q","--quiet",action="store_true",dest="quiet",default=False,remember=True,
                  help="Quiet.",name="Quiet")
parser.add_option("-o","--output",type="savefile",dest="output",remember=True,
                  help="Output file. Leave empty for console ouptut.",default="",
                  name="Output File", filetypes=[("All output formats","*.xlsx;*.xls;*.csv;*.tsv;*.txt"),
                                                 ("Excel","*.xlsx"),("Excel2003","*.xls"),
                                                 ("CSV","*.csv"),("TSV","*.tsv"),("Text","*.txt")])
parser.add_option_group(advanced)
opt = None
while True:
    if 'exit' in error_kwargs:
        try:
            opt,args = parser.parse_args(opts=opt)
        except UserCancelledError:
            sys.exit(0);
    else:
        opt,args = parser.parse_args()

    break

if not opt.output:
    opt.quiet = True
				
if opt.filtr == True:					
  for i in enumerate(opt.snps):		
    P = Popen('./exonicFilter.py' + ' ' + i[1], stdout=PIPE , stderr=STDOUT, shell=True)	
    P.communicate()[0]				
  l_file = []
  for subdir, dirs, files in os.walk('./'):
    for file in files:		
      if 'Filtr' in file and file.endswith('.vcf'):
        	l_file.append(file)
  if opt.output:			
   p_rds = Popen('./readCounts.py' + ' --'   + ' --'.join(l)+ ' -r "' + ' '.join(opt.alignments) + '" -s'  + ' "' + ' '.join(l_file) + '" -o' + ' ' + opt.output, shell=True)
   p_rds.communicate()[0]
   p3 = Popen('./snp_computation.py' + ' ' + opt.output ,stdout=PIPE,shell=True)
   p3.communicate()[0]
   for subdir, dirs, files in os.walk('./'):
    for file in files:
      if file.startswith('Event'):
         summary_analysis.read_events(file)														
   print "Done processing"
  else:											
      f = open("result_readscount.tsv", 'w')												
      p2= Popen('./readCounts.py' + ' --'   + ' --'.join(l)+ ' -r "' + ' '.join(opt.alignments) + '" -s'  + ' "' + ' '.join(l_file)  + '" ',shell=True,stdout=f)
      p2.communicate()[0]			
      f.close()											
      outfile = open('result_readscounts_final.tsv', 'w')
      infile = open('result_readscount.tsv')
      for line in infile:
        if line.startswith('CHROM'):
          outfile.write(line)
        else:									
          fields = line.strip().split()
          print >>outfile, '\t'.join(fields)
      outfile.close()
      infile.close()

      p3 = subprocess.Popen('./snp_computation.py' + ' ' +  "result_readscounts_final.tsv",shell=True)
      p3.communicate()[0]
      for subdir, dirs, files in os.walk('./'):
       for file in files:
	if file.startswith('Event'):
         summary_analysis.read_events(file)
      print "Done processing" 
else:
 if opt.output:																					               
	p_rds = Popen('./readCounts.py' +' --'   + ' --'.join(l)+ ' -r "' + ' '.join(opt.alignments) + '" -s'  + ' "'.join(opt.snps) + '" -o' + ' ' + opt.output,shell=True)
	p_rds.communicate()[0]
   	p3 = Popen('./snp_computation.py' + ' ' + opt.output ,stdout=PIPE,shell=True)
   	output = p3.communicate()[0]
   	for subdir, dirs, files in os.walk('./'):
    	 for file in files:
	  if file.startswith('Event'):
	    summary_analysis.read_events(file)		
   	print "Done processing"
		
 else:
																	
      f = open("result_readscount.tsv", 'w')												
      p2= Popen('./readCounts.py' + ' --'   + ' --'.join(l)+ ' -r "' + ' '.join(opt.alignments) + '" -s'  + ' "' + ' '.join(opt.snps)  + '" ',shell=True,stdout=f)
      p2.communicate()[0]
      f.close()
      with open('result_readscounts_final.tsv', 'w') as outfile:
       with open('result_readscount.tsv') as infile:
         for line in infile:
          if line.startswith('CHROM'):
             outfile.write(line)
          else:
             fields = line.strip().split()
             print >>outfile, '\t'.join(fields)

      p3 = subprocess.Popen('./snp_computation.py' + ' ' +  "result_readscounts_final.tsv",shell=True)
      p3.communicate()[0]
      for subdir, dirs, files in os.walk('./'):
       for file in files:
        if file.startswith('Event'):
         summary_analysis.read_events(file)
      print "Done processing"
