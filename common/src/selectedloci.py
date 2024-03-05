#!/bin/env python3

from pysamimport import pysam
from heapq import heappop, heappush
from collections import defaultdict
from operator import itemgetter
import os, os.path, copy, traceback, re, time

from util import BasicFilter, BadRead, ReadGroupFactory

class SelectedLoci(object):
    def __init__(self, filename, selectedloci, **kw):
        self.filename = filename
        self.minbasequal = kw.get('minbasequal',25)
        self.minmappingqual = kw.get('minmappingqual',25)
        self.gapsize = kw.get('gapsize',1000)
        self.maxedits = kw.get('maxedits',3)
        self.quickcnt = kw.get('quickcnt',True)
        self.selloci = []
        self.regions = []
        self.refalt = {}
        last = (None,None)
        for ch,pos,ref,alt in sorted(map(lambda t: (t[0],int(t[1])-1,t[2],t[3]),selectedloci)):
            if ch != last[0]:
                if last[0] != None:
                    self.regions[-1][2] = last[1]
                self.regions.append([ch,pos,None])
                self.selloci.append(set([pos]))
            elif (pos - last[1]) > self.gapsize:
                self.regions[-1][2] = last[1]
                self.regions.append([ch,pos,None])
                self.selloci.append(set([pos]))
            else:
                self.selloci[-1].add(pos)
            if (ch,pos) not in self.refalt:
                self.refalt[ch,pos] = [ref,alt]
            else: 
                self.refalt[ch,pos].append(alt)
            last = (ch,pos)
        self.regions[-1][2] = last[1]

    def loci(self):
        samfile = pysam.AlignmentFile(self.filename, "rb", require_index=True)
        starttime = time.time()
        totalloci = 0
        allmatch = re.compile(r'^\d+(S\d+)?M(\d+S)?$')
        badcigar = re.compile(r'[NDI]')
        for (ch,start,end),selected in zip(self.regions,self.selloci):
          #print(">"+ch,start+1,end+1,len(selected))
          self.clear()
          self._minlocus = 0
          nselected = len(selected)
          maxgap = 0
          if nselected > 1:
              sellist = sorted(selected)
              maxgap = max(sellist[i+1]-sellist[i] for i in range(nselected-1))
          totalalcnt = 0
          matchalcnt = 0
          emptyalcnt = 0
          quickalcnt = 0
          quickemptycnt = 0
          for al in samfile.fetch(contig=ch,start=start,end=end+1):
            totalalcnt += 1
            while self._minlocus < al.reference_start and self._minlocus:
                for rpos, freq in self.process_minlocus():
                    yield ch, rpos+1, self.refalt[ch,rpos], freq
            if al.is_duplicate or al.is_qcfail or al.is_secondary or al.is_unmapped:
                continue
            if al.mapping_quality < self.minmappingqual:
                continue
            if al.get_tag('NM') > self.maxedits:
                continue
            if maxgap > al.reference_length and not any(map(lambda p: al.reference_start <= p < al.reference_end, selected)):
                continue
            if self.quickcnt and not badcigar.search(al.cigarstring) and allmatch.search(al.cigarstring):
                almd = None
                anymatch = False
                for rpos in selected:
                    qpos = al.query_alignment_start+(rpos-al.reference_start)
                    if al.query_alignment_start <= qpos < al.query_alignment_end:
                        anymatch = True
                        if al.query_qualities[qpos] < self.minbasequal:
                            continue
                        if not almd:
                            almd = self.getalmd(al)
                        if rpos not in self._loci:
                            if rpos < self._minlocus or self._minlocus == 0:
                                self._minlocus = rpos
                        qbase = al.query_sequence[qpos]
                        self.countread(rpos,qbase,almd)
                if anymatch:
                    quickalcnt += 1
                else:
                    quickemptycnt += 1
                continue

            if self.quickcnt and not badcigar.search(al.cigarstring):
                print(al)

            almd = None
            anymatch = False
            for qpos,rpos in filter(lambda t: t[1] in selected, al.get_aligned_pairs(matches_only=True)):
                anymatch = True
                if al.query_qualities[qpos] < self.minbasequal:
                    continue
                if not almd:
                    almd = self.getalmd(al)
                if rpos not in self._loci:
                    if rpos < self._minlocus or self._minlocus == 0:
                        self._minlocus = rpos
                qbase = al.query_sequence[qpos]
                self.countread(rpos,qbase,almd)
            if anymatch:
                matchalcnt += 1
            else:
                emptyalcnt += 1
          while self._minlocus:
            for rpos, freq in self.process_minlocus():
                yield ch, rpos+1, self.refalt[ch,rpos], freq
          # elapsed = (time.time()-starttime)
          # totalloci += len(selected)
          # print(">> %.2f%% %.2f%% %.2f%% %.2f%% %.2f%% %d reads %d loci %d sec %.2f loci/sec"%(100*quickalcnt/totalalcnt,100*quickemptycnt/totalalcnt,100*matchalcnt/totalalcnt,100*emptyalcnt/totalalcnt,100*(quickalcnt+quickemptycnt+matchalcnt+emptyalcnt)/totalalcnt,totalalcnt,totalloci,elapsed,totalloci/elapsed))

    def getalmd(self,al):
        return None

    def process_minlocus(self):
        rpos = self._minlocus
        d = self._loci[self._minlocus]
        del self._loci[self._minlocus]
        if len(self._loci) == 0:
            self._minlocus = 0
        else:
            self._minlocus = min(self._loci)
        if self.accept(d):
            yield rpos,d

class SelectedLociByFetch(SelectedLoci):

    def __init__(self,*args,**kw):
        self.mincnt = kw.get('minreads',5)
        super().__init__(*args,**kw)

    def clear(self):
        self._loci = defaultdict(lambda: defaultdict(int))

    def getalmd(self,al):
        # return ("R" if al.is_reverse else "F")
        return None

    def countread(self,rpos,qbase,almd):
        self._loci[rpos][qbase] += 1

    def accept(self,freq):
        return sum(freq[nuc] for nuc in "ACGT") >= self.mincnt

class SCSelectedLociByFetch(SelectedLoci):

    def __init__(self, *args, **kw):
        self.cellbarcodes = kw.get('cellbarcodes',None)
        self.umis = kw.get('umibarcodes',None)
        self.mincells = kw.get('mincells',None)
        self.minumipercell = kw.get('minumipercell',None)
        self.minreadspercell = kw.get('minreadspercell',None)
        super().__init__(*args,**kw)

    def clear(self):
        self._loci = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))

    def getalmd(self,al):
        return self.cellbarcodes.group(al),self.umis.group(al)

    def countread(self,rpos,qbase,almd):
        cb,umi = almd
        if not cb or not umi:
            return
        self._loci[rpos][cb][qbase][umi] += 1

    def accept(self,freq):
        goodcb = 0
        for cb in freq:
            if (self.minumipercell > 0) and (sum([len(freq[cb][nuc]) for nuc in 'ACGT']) < self.minumipercell):
               continue
            if (self.minreadspercell > self.minumipercell) and (sum([sum(freq[cb][nuc].values()) for nuc in 'ACGT']) < self.minreadspercell):
               continue
            goodcb += 1
            if goodcb >= self.mincells:
                return True
        return False

from rcio import scrcbinwriter

if __name__ == "__main__":
    import sys, time

    bf = sys.argv[1]
    sf = sys.argv[2]
    kwargs = {'minreads': 5, 'minbasequal': 25, 'minmappingqual': 60, 'gapsize': 1000, 'maxedits': 3, 
              'cellbarcodes': 'STARsolo', 'mincells': 3, 'minvarumipercell': 3, 'minumipercell': 3}
    for kv in sys.argv[3:]:
        k,v = kv.split('=')
        try:
           v = eval(v)
        except:
           pass
        kwargs[k] = v
    snvs = [ l.split()[:4] for l in open(sf) ]

    groupFactory = ReadGroupFactory()
    groupOptions = [t[1] for t in groupFactory.list(type="CellBarcode")] + [""]
    groupMap = {}
    for s,n,d in sorted(groupFactory.list(type="CellBarcode")):
        groupMap[n] = s
    
    umiOptions = [t[1] for t in groupFactory.list(type="UMI")]
    umiMap = {}
    for s,n,d in sorted(groupFactory.list(type="UMI")):
        umiMap[n] = s

    sc = False
    if kwargs.get('cellbarcodes',None) != None:
        sc = True
        if not kwargs.get('umibarcodes',None):
            kwargs['umibarcodes'] = kwargs['cellbarcodes']
        if kwargs.get('acceptlist') != None:
            if kwargs.get('acceptlist') in ("","None","-"):
                readgroupparam = "*:acceptlist=None"
            else:
                readgroupparam = "*:acceptlist='%s'"%(kwargs.get('acceptlist',))
        else:
            readgroupparam = ""
        kwargs['cellbarcodes'] = groupFactory.get(groupMap[kwargs['cellbarcodes']],readgroupparam)
        kwargs['umibarcodes'] = groupFactory.get(umiMap[kwargs['umibarcodes']])

    start = time.time()
    if not sc:
        locifreq = SelectedLociByFetch(bf,snvs,**kwargs)
    else:
        locifreq = SCSelectedLociByFetch(bf,snvs,**kwargs)
    writer = scrcbinwriter(bf.rsplit('.',1)[0]+'.rcbin')
    for chr,pos,refalt,freq in locifreq.loci():
        ref = refalt[0]
        if not sc:
            for alt in refalt[1:]:
                print(chr,pos,refalt[0],alt,freq[refalt[0]],freq[alt])
        else:
            for alt in refalt[1:]:
                for cb in freq:
                    allumi = set()
                    vals = dict()
                    for nuc in 'ACGT':
                        allumi.update(freq[cb][nuc])
                        vals[nuc] = len(freq[cb][nuc])
                    if len(allumi) < kwargs.get('minumipercell'):
                        continue
                    if vals[alt] < kwargs.get('minvarumipercell'):
                        continue
                    writer.writecounts(chr,pos,ref,alt,None,cb,vals[ref],vals[alt])
    writer.close()
    elapsed = time.time() - start
    print("Elapsed: %d:%02d.%02d"%(int(elapsed//60),int(elapsed%60),int(round(100*(elapsed-int(elapsed))))))
