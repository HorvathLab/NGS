#!/bin/env python3

from pysamimport import pysam
from heapq import heappop, heappush
from collections import defaultdict
import os, os.path

from util import BasicFilter, BadRead

class VariantLoci(object):
    def __init__(self, filename, **kw):
        self.filename = filename
        self.region = kw.get('region',None)
        self.mincnt = kw.get('mincnt',None)

class VariantLociByFetch(VariantLoci):
    def loci(self):
        samfile = pysam.AlignmentFile(self.filename, "rb", require_index=True)
        self.loci = defaultdict(lambda: defaultdict(int))
        self.minlocus = None
        for al in samfile.fetch(region=self.region):
            while self.minlocus != None and self.minlocus[:2] < (al.reference_id,al.reference_start):
                for rname, rpos, ref, alt, freq in self.process(self.minlocus):
                    yield rname, rpos, ref, alt, freq
            if al.is_duplicate or al.is_qcfail or al.is_secondary or al.is_unmapped or al.get_tag('NM') == 0:
                continue
            for qpos,rpos,rbase in filter(lambda t: t[2] in 'acgt', al.get_aligned_pairs(matches_only=True,with_seq=True)):
                rpos = (al.reference_id,rpos,al.reference_name)
                if rpos not in self.loci:
                    self.loci[rpos]['ref'] = rbase.upper()
                    if self.minlocus == None or self.minlocus > rpos:
                        self.minlocus = rpos
                qbase = al.query_sequence[qpos]
                self.loci[rpos][qbase] += 1
        while self.minlocus != None:
            for rname, rpos, ref, alt, freq in self.process(self.minlocus):
                    yield rname, rpos, ref, alt, freq
    
    def process(self,locus):
        d = self.loci[locus]
        del self.loci[locus]
        if locus == self.minlocus:
            if len(self.loci) == 0:
                self.minlocus = None
            else:
                self.minlocus = min(self.loci)
        for alt in d:
            if alt == 'ref':
                ref = d['ref']
                continue
            if d[alt] < self.mincnt:
                continue
            yield locus[2],locus[1],ref,alt,{alt: d[alt]}

class VariantLociByPileup(VariantLoci):

    def loci(self):
        extract_base = BasicFilter().extract_base
        samfile = pysam.AlignmentFile(self.filename, "rb", require_index=True)
        pileup_kwargs = dict(stepper = 'all', # stepper='nofilter',
                             ignore_overlaps=True,
                             min_base_quality=0)
        for pileupcolumn in samfile.pileup(region=self.region, truncate=True, **pileup_kwargs):
             if pileupcolumn.n < self.mincnt:
                 continue
             freq = defaultdict(int) 
             refnuc = None
             for pileupread in pileupcolumn.pileups:
                 try:
                     al, qpos, nuc = extract_base(pileupread)
                 except BadRead:
                     continue
                 if refnuc == None:
                     for t in filter(lambda t: t[0] == qpos,al.get_aligned_pairs(matches_only=True,with_seq=True)):
                         refnuc = t[2].upper()
                         break
                 freq[nuc] += 1
             if refnuc == None:
                 continue
             if len(freq) == 0 or (len(freq) == 1 and refnuc in freq):
                 continue
             sortednucs = sorted(freq.items(),key=lambda t: t[1],reverse=True)
             ref,rcnt = refnuc,freq[refnuc]
             for alt,acnt in sortednucs:
                 if alt == ref:
                     continue
                 if acnt < self.mincnt:
                     break
                 yield pileupcolumn.reference_name, pileupcolumn.reference_pos, ref, alt, freq
