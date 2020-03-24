
import threading
from collections import Counter
from pysamimport import pysam
from util import BadRead
from Queue import Queue
import time, math, sys

class Pileups(object):
    def __init__(self,loci,samfiles,filter,chrreg):
        self.loci = loci
        self.samfiles = samfiles
        self.filter = filter
	self.chrreg = chrreg

class SerialPileups(Pileups):

    def iterator(self):
        samfiles = []
        chrommap = []
        for al in self.samfiles:
            samfile = pysam.Samfile(al, "rb")
	    assert samfile.has_index(), "Cannot open BAM index for file %s"%sf
            samfiles.append(samfile)
            chrommap.append(self.chrreg.chrommap(al))
            
        for snvchr, snvpos, ref, alt, snvextra in self.loci:
            cnts = Counter()
            total = Counter()
            reads = []
            snvpos1 = snvpos - 1
            for i, samfile in enumerate(samfiles):
		try:
		    snvlabel = chrommap[i](snvchr)
		    if snvlabel != None:
                      for pileupcolumn in samfile.pileup(snvlabel, snvpos1, snvpos1 + 1, truncate=True):
                        total[i] += pileupcolumn.n
                        for pileupread in pileupcolumn.pileups:
                            try:
                                al, pos, base, nseg = self.filter.test(pileupread)
                            except BadRead, e:
                                cnts[(i, e.message)] += 1
                                continue
                            reads.append((al, pos, base, i))
                            cnts[(i, 'Good')] += 1
		except ValueError:
		    pass
	    total[i] -= cnts[(i,"GapInQueryAtSNVLocus")]
	    # del cnts[(i,"GapInQueryAtSNVLocus")]
            yield (snvchr, snvpos, ref, alt, total, reads, cnts)

class ThreadedPileups(Pileups):
    def __init__(self,*args,**kw):
        super(ThreadedPileups,self).__init__(*args)
        self.tpb = kw.get('threadsperbam',1)
        self.nb = len(self.samfiles)
        self.nt = self.tpb*self.nb
        self._queue = []
        k = 0;
        for j in range(self.tpb):
            for i in range(self.nb):
                self._queue.append(Queue(20))
                t = threading.Thread(target=self.worker,args=(i,j,k))
		t.daemon = True
                t.start()
                k += 1
	    time.sleep(1)

    def worker(self,i,j,k):
        samfile = pysam.Samfile(self.samfiles[i], "rb")
	assert samfile.has_index(), "Cannot open BAM index for file %s"%sf
	chrommap = self.chrreg.chrommap(self.samfiles[i])
	# blocksize = int(math.ceil(len(self.loci)/self.tpb))
	# for l in range(j*blocksize,min((j+1)*blocksize,len(self.loci))):
	for l in range(j,len(self.loci),self.tpb):
	    snvchr, snvpos, ref, alt, snvextra = self.loci[l]
            cnts = Counter()
            total = Counter()
            reads = []
            snvpos1 = snvpos - 1
	    try:
		snvlabel = chrommap(snvchr)
		if snvlabel != None:
                  for pileupcolumn in samfile.pileup(snvlabel, snvpos1, snvpos1 + 1, truncate=True):
                    total[i] += pileupcolumn.n
                    for pileupread in pileupcolumn.pileups:
                        try:
                            al, pos, base, nseg = self.filter.test(pileupread)
                        except BadRead, e:
                            cnts[(i, e.message)] += 1
                            continue
                        reads.append((al, pos, base, i))
                        cnts[(i, 'Good')] += 1
	    except ValueError, e:
	        pass # raise e
	    total[i] -= cnts[(i,"GapInQueryAtSNVLocus")]
	    # del cnts[(i,"GapInQueryAtSNVLocus")]
	    # print >>sys.stderr, (snvchr, snvpos, ref, alt, total, cnts)
            self._queue[k].put((snvchr, snvpos, ref, alt, total, reads, cnts))
        return
        
    def iterator(self):
        k = 0
	# for i in range(len(self.loci)):
        for snvchr, snvpos, ref, alt, snvextra in self.loci:
            cnts = Counter()
            total = Counter()
            reads = []
            for i in range(len(self.samfiles)):
                snvchri, snvposi, refi, alti, totali, readsi, cntsi = self._queue[k].get()
                assert(snvchri == snvchr and snvposi == snvpos)
                assert(i in totali or len(totali.keys()) == 0)
                reads.extend(readsi)
                cnts.update(cntsi)
                total.update(totali)
                self._queue[k].task_done()
                k = (k+1)%self.nt
            yield (snvchr, snvpos, ref, alt, total, reads, cnts)

