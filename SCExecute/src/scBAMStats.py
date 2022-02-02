#!/bin/env python3
import sys
import os
import os.path
import traceback
from os.path import join, dirname, realpath
try:
    sys.path.append(join(dirname(realpath(__file__)),
                         '..', '..', 'common', 'src'))
except NameError:
    pass

from collections import defaultdict
from pysamimport import pysam

nreads = 0
cellumi2readcount = defaultdict(lambda: defaultdict(int))
samfile = pysam.AlignmentFile(sys.argv[1], "rb", require_index=False)
for al in samfile.fetch(until_eof=True):
    cellbc = al.get_tag('CB')
    umibc = al.get_tag('UB')
    cellumi2readcount[cellbc][umibc] += 1

umidepthfreq = defaultdict(lambda: defaultdict(int))
readdepthfreq = defaultdict(lambda: defaultdict(int))
for pileupcolumn in samfile.pileup(stepper='all'):
    umis = defaultdict(set)
    reads = defaultdict(int)
    for pileupread in pileupcolumn.pileups:
      al = pileupread.alignment
      cellbc = al.get_tag('CB')
      umibc = al.get_tag('UB')
      umis[cellbc].add(umibc)
      reads[cellbc] += 1
    for cell in umis:
      depth = len(umis[cell])
      umidepthfreq[cell][depth] += 1
      readdepthfreq[cell][reads[cell]] += 1

def summarize(label,values):
    metrics = {}
    metrics['distinct'] = len(values)
    metrics['count'] = sum(values.values())
    metrics['min'] = min(values)
    metrics['max'] = max(values)
    prevvalue = metrics['min']
    cumindex = 0
    p5 = p10 = p90 = p95 = None
    q1 = None
    med = None
    q3 = None
    for v,f in sorted(values.items()):
        cumindex += f
        if p5 == None and cumindex/metrics['count'] > 0.05:
            p5 = prevvalue
        if p10 == None and cumindex/metrics['count'] > 0.1:
            p10 = prevvalue
        if q1 == None and cumindex/metrics['count'] > 0.25:
            q1 = prevvalue
        if med == None and cumindex/metrics['count'] > 0.50:
            med = prevvalue
        if q3 == None and cumindex/metrics['count'] > 0.75:
            q3 = prevvalue
        if p90 == None and cumindex/metrics['count'] > 0.90:
            p90 = prevvalue
        if p95 == None and cumindex/metrics['count'] > 0.95:
            p95 = prevvalue
        prevvalue = v
    metrics['median'] = med
    metrics['q1'] = q1
    metrics['q3'] = q3
    metrics['p5'] = p5
    metrics['p10'] = p10
    metrics['p90'] = p90
    metrics['p95'] = p95
    metrics['mean'] = round(sum(v*f for v,f in values.items())/metrics['count'],2)
    print("%s:"%(label,),end="")    
    for m in ('count','min','q1','mean','median','q3','p90','p95','max'):
        print(" %s %s"%(m,metrics[m]),end="")
    print()

# print("Cells:",len(cells))
for cell in cellumi2readcount:
    readperumifreq = defaultdict(int)
    for umi in cellumi2readcount[cell]:
       readperumifreq[cellumi2readcount[cell][umi]] += 1
    summarize(cell+" reads per UMI",readperumifreq)
    summarize(cell+" UMI depth per position",umidepthfreq[cell])
    summarize(cell+" read depth per position",readdepthfreq[cell])
