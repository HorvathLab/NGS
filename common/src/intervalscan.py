#!/bin/env python3

from pysamimport import pysam
from heapq import heappop, heappush
from collections import defaultdict

class IntervalScan(object):
    def __init__(self, filename, **kw):
        self.filename = filename

    def process(self):
        current_name = None
        current_ends = []
        alindex = 0
        self.samfile = pysam.AlignmentFile(self.filename, "rb", require_index=False)
        self.process_init()
        for al in self.samfile.fetch(until_eof=True):
            if al.reference_name != current_name:
                while len(current_ends) > 0:
                    pos,ind,al1 = heappop(current_ends)
                    self.process_end(al1)
                # reset
                current_name = al.reference_name
                current_ends = []
                last_start = -1
            assert al.reference_start >= last_start
            current_start = al.reference_start 
            while (len(current_ends) > 0 and current_ends[0][0] < current_start):
                pos,ind,al1 = heappop(current_ends)
                assert pos < current_start
                self.process_end(al1)
            heappush(current_ends,(al.reference_end,alindex,al))
            alindex += 1
            self.process_start(al)
            last_start = current_start
        while len(current_ends) > 0:
            pos,ind,al1 = heappop(current_ends)
            self.process_end(al1)
        self.process_done()

    def process_init(self):
        pass
 
    def process_start(self,alignment):
        pass

    def process_end(self,alignment):
        pass

    def process_done(self):
        pass

class Intervals(IntervalScan):
    def __init__(self,filename,**kw):
        self.mindepth = kw.get('mindepth',1)
        super(Intervals,self).__init__(filename,**kw)

    def process_init(self):
        self.depth = 0
        self.current = set()

    def process_start(self,al):
        self.depth += 1
        self.current.add(al)
        if self.depth < self.mindepth:
            return
        if self.depth == self.mindepth:
            self.interval_start = al.reference_start
            self.interval_current = set(self.current)
            self.maxdepth = self.depth
        else:
            self.maxdepth = max(self.maxdepth,self.depth)
            self.interval_current.add(al)

    def process_end(self,al):
        if self.depth == self.mindepth:
            interval_end = al.reference_end
            print(al.reference_name,
                  self.interval_start,
                  interval_end,
                  len(self.interval_current),
                  interval_end-self.interval_start,
                  self.maxdepth)
        self.depth -= 1
        self.current.remove(al)

class MinDepthBAM(Intervals):

    def process_init(self):
        super(MinDepthBAM,self).process_init()
        self.outreads = set()

    def process_end(self,al):
        if self.depth == self.mindepth:
            self.outreads.update(self.interval_current)
            # for al1 in sorted(self.interval_current,key=lambda al1: (al1.reference_start,al1.reference_end)):
            #     self.outbam.write(al1)
        super(MinDepthBAM,self).process_end(al)

    def process_done(self):
        outfilename = os.path.split(self.filename.rsplit('.',1)[0] + ".d%d.bam"%(self.mindepth,))[1]
        outbam = pysam.AlignmentFile(outfilename,'wb',template=self.samfile)
        for al1 in sorted(self.outreads,key=lambda al1: (al1.reference_id,al1.reference_start,al1.reference_end)):
            outbam.write(al1)
        outbam.close()
        pysam.index(outfilename)
        super(MinDepthBAM,self).process_done()
