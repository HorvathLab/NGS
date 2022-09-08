#!/bin/env python3

from pysamimport import pysam
from heapq import heappop, heappush
from collections import defaultdict
import os, os.path, copy, traceback

from util import BasicFilter, BadRead

class VariantLoci(object):
    def __init__(self, filename, **kw):
        self.filename = filename
        self.region = kw.get('region',None)
        self.mincnt = kw.get('minreads',0)
        self.minqual = kw.get('minqual',0)
        self.mindist = kw.get('mindist',0)
        self._loci = defaultdict(lambda: defaultdict(int))
        self._minlocus = None
        self._outqueue = []

    def loci(self):
        samfile = pysam.AlignmentFile(self.filename, "rb", require_index=True)
        has_NM = None; has_nM = None
        for al in samfile.fetch(region=self.region):
            while self._minlocus != None and self._minlocus[:2] < (al.reference_id,al.reference_start):
                for rname, rpos, ref, alt, freq in self.process(self._minlocus):
                    yield rname, rpos+1, ref, alt, freq
            if al.is_duplicate or al.is_qcfail or al.is_secondary or al.is_unmapped:
                continue
            if has_NM == None:
                has_NM = al.has_tag('NM')
                has_nM = al.has_tag('nM')
            if has_NM and al.get_tag('NM') == 0:
                continue
            if has_nM and al.get_tag('nM') == 0:
                continue
            first = True
            almd = None
            for qpos,rpos,rbase in filter(lambda t: t[2] in 'acgt', al.get_aligned_pairs(matches_only=True,with_seq=True)):
                # print(qpos,rpos,rbase,al.query_sequence[qpos])
                if al.query_qualities[qpos] < self.minqual:
                    continue
                if first:
                    almd = self.getalmd(al)
                rpos = (al.reference_id,rpos,al.reference_name)
                if rpos not in self._loci:
                    self._loci[rpos]['ref'] = rbase.upper()
                    if self._minlocus == None or self._minlocus > rpos:
                        self._minlocus = rpos
                qbase = al.query_sequence[qpos]
                assert(qbase != rbase.upper())
                self.countread(rpos,qbase,almd)
        while self._minlocus != None:
            for rname, rpos, ref, alt, freq in self.process(self._minlocus):
                yield rname, rpos+1, ref, alt, freq
        for rname, rpos, ref, alt, freq in self.process((None,None,None)):
            yield rname, rpos+1, ref, alt, freq

    def getalmd(self,al):
        return None

    def process(self,locus):
        assert len(self._outqueue) in (0,1,2,3)
        if len(self._outqueue) == 0 or locus[2] != self._outqueue[0][0]:
            for i in range(1,len(self._outqueue)-1):
                if (self._outqueue[i][1] - self._outqueue[i-1][1]) >= self.mindist and \
                   (self._outqueue[i+1][1] - self._outqueue[i][1]) >= self.mindist:
                    yield self._outqueue[i]
            i = len(self._outqueue)-1
            if i >= 1 and abs(self._outqueue[i][1] - self._outqueue[i-1][1]) >= self.mindist:
                yield self._outqueue[i]
            if locus[2] == None:
                return
            self._outqueue = [ (locus[2],0) ]

        d = self._loci[locus]
        del self._loci[locus]
        if locus == self._minlocus:
            if len(self._loci) == 0:
                self._minlocus = None
            else:
                self._minlocus = min(self._loci)
        for alt in d:
            if alt == 'ref':
                ref = d['ref']
                continue
            if not self.accept(d[alt]):
                continue
            self._outqueue.append((locus[2],locus[1],ref,alt,{alt: copy.deepcopy(d[alt])}))
            if len(self._outqueue) == 3:
                if (self._outqueue[1][1] - self._outqueue[0][1]) >= self.mindist and \
                   (self._outqueue[2][1] - self._outqueue[1][1]) >= self.mindist:
                    yield self._outqueue[1]
                self._outqueue.pop(0)

class SCVariantLoci(VariantLoci):
    def __init__(self, filename, **kw):
        super().__init__(filename,**kw)
        self._loci = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
        self.cellbarcodes = kw.get('cellbarcodes',None)
        self.umis = kw.get('umis',None)
        self.mincells = kw.get('mincells',None)
        self.minumipercell = kw.get('minumipercell',None)

class VariantLociByFetch(VariantLoci):

    def countread(self,rpos,qbase,almd):
        self._loci[rpos][qbase] += 1

    def accept(self,cnt):
        return cnt >= self.mincnt

class SCVariantLociByFetch(SCVariantLoci):

    def getalmd(self,al):
        return self.cellbarcodes.group(al),self.umis.group(al)

    def countread(self,rpos,qbase,almd):
        cb,umi = almd
        if not cb or not umi:
            return
        self._loci[rpos][qbase][cb].add(umi)

    def accept(self,cells):
        goodcb = 0
        for cb in cells:
            if len(cells[cb]) >= self.minumipercell:
                goodcb += 1
                if goodcb >= self.mincells:
                    return True
        return False
    
class VariantLociByPileup(VariantLoci):

    def loci(self):
        extract_base = BasicFilter().extract_base
        samfile = pysam.AlignmentFile(self.filename, "rb", require_index=True)
        pileup_kwargs = dict(stepper = 'nofilter',
                             max_depth = 100000,
                             ignore_overlaps=True,
                             ignore_orphans=False,
                             min_base_quality=0,
                             flag_filter=0)
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
                 yield pileupcolumn.reference_name, pileupcolumn.reference_pos+1, ref, alt, copy.deepcopy(freq)
