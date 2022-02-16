
import threading
import multiprocessing
from collections import Counter, namedtuple
from pysamimport import pysam
from util import BadRead
import queue
import time, math, sys

class Pileups(object):
    def __init__(self,loci,samfiles,filter,chrreg,readgroups):
        self.loci = loci
        self.samfiles = samfiles
        self.filter = filter
        self.chrreg = chrreg
        self.readgroups = readgroups

class SerialPileups(Pileups):

    def iterator(self):
        samfiles = []
        chrommap = []
        for al in self.samfiles:
            samfile = pysam.Samfile(al, "rb")
            assert samfile.has_index(), "Cannot open BAM index for file %s"%sf
            samfiles.append(samfile)
            chrommap.append(self.chrreg.chrommap(al))

        pileup_kwargs = self.filter.pileup_kwargs()
        for snvchr, snvpos, ref, alt, snvextra in self.loci:
            cnts = Counter()
            total = Counter()
            reads = []
            snvpos1 = snvpos - 1
            for i, samfile in enumerate(samfiles):
                try:
                    snvlabel = chrommap[i](snvchr)
                    if snvlabel != None:
                        for pileupcolumn in samfile.pileup(snvlabel, snvpos1, snvpos1 + 1, truncate=True, **pileup_kwargs):
                            self.filter.pileup_start(pileupcolumn)
                            total[i] = pileupcolumn.n
                            for pileupread in pileupcolumn.pileups:
                                if self.readgroups != None:
                                    grp = (i,self.readgroups.group(pileupread.alignment))
                                    if grp[1] == None:
                                        continue
                                    total[grp] += 1
                                else:
                                    grp = i
                                try:
                                    al, pos, base = self.filter.extract_base(pileupread)
                                except BadRead as e:
                                    cnts[(grp, e.args[0])] += 1
                                    continue
                                reads.append((al, pos, base, grp))
                                cnts[(grp, 'Good')] += 1
                            self.filter.pileup_end(pileupcolumn)
                except ValueError:
                    pass
            yield (snvchr, snvpos, ref, alt, total, reads, cnts)

class ThreadedPileups(Pileups):
    def __init__(self,*args,**kw):
        super(ThreadedPileups,self).__init__(*args)
        self.nt = kw.get('threads',1)
        self.nb = len(self.samfiles)
        if self.nt < self.nb or self.nt % self.nb != 0:
            raise RuntimeError("Number of threads should be a (non-negative integer) multiple of the number BAM files")
        self.tpb = self.nt//self.nb
        self._queue = []
        k = 0;
        for j in range(self.tpb):
            for i in range(self.nb):
                self._queue.append(queue.Queue(20))
                t = threading.Thread(target=self.worker,args=(i,j,k))
                t.daemon = True
                t.start()
                k += 1
            time.sleep(1)

    def worker(self,i,j,k):
        samfile = pysam.Samfile(self.samfiles[i], "rb")
        assert samfile.has_index(), "Cannot open BAM index for file %s"%sf
        chrommap = self.chrreg.chrommap(self.samfiles[i])
        pileup_kwargs = self.filter.pileup_kwargs()
        for l in range(j,len(self.loci),self.tpb):
            snvchr, snvpos, ref, alt, snvextra = self.loci[l]
            cnts = Counter()
            total = Counter()
            reads = []
            snvpos1 = snvpos - 1
            try:
                snvlabel = chrommap(snvchr)
                if snvlabel != None:
                  for pileupcolumn in samfile.pileup(snvlabel, snvpos1, snvpos1 + 1, truncate=True, **pileup_kwargs):
                    self.filter.pileup_start(pileupcolumn)
                    total[i] += pileupcolumn.n
                    for pileupread in pileupcolumn.pileups:
                        if self.readgroups != None:
                            grp = (i,self.readgroups.group(pileupread.alignment))
                            if grp[1] == None:
                                continue
                            total[grp] += 1
                        else:
                            grp = i
                        try:
                            al, pos, base = self.filter.extract_base(pileupread)
                        except BadRead as e:
                            cnts[(grp, e.args[0])] += 1
                            continue
                        reads.append((al, pos, base, grp))
                        cnts[(grp, 'Good')] += 1
                    self.filter.pileup_end(pileupcolumn)
            except ValueError as e:
                pass
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
                assert(i in totali or len(list(totali.keys())) == 0)
                reads.extend(readsi)
                cnts.update(cntsi)
                total.update(totali)
                self._queue[k].task_done()
                k = (k+1)%self.nt
            yield (snvchr, snvpos, ref, alt, total, reads, cnts)

PileupAlignment = namedtuple('PileupAlignment',['seq','is_reverse'])

class MultiprocPileups(Pileups):
    # A python class for alignments since the C++ wrapped data-structure
    # doesn't survive the multiprocess communication process...
    def __init__(self,*args,**kw):
        super(MultiprocPileups,self).__init__(*args)
        self.nt = kw.get('procs',1)
        self.nb = len(self.samfiles)
        if self.nt < self.nb or self.nt % self.nb != 0:
            raise RuntimeError("Number of threads should be a (non-negative integer) multiple of the number BAM files")
        self.tpb = self.nt//self.nb
        self._queue = []
        k = 0;
        for j in range(self.tpb):
            for i in range(self.nb):
                self._queue.append(multiprocessing.Queue(20))
                t = multiprocessing.Process(target=self.worker,args=(i,j,k))
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
        pileup_kwargs = self.filter.pileup_kwargs()
        for l in range(j,len(self.loci),self.tpb):
            snvchr, snvpos, ref, alt, snvextra = self.loci[l]
            cnts = Counter()
            total = Counter()
            reads = []
            snvpos1 = snvpos - 1
            try:
                snvlabel = chrommap(snvchr)
                if snvlabel != None:
                  for pileupcolumn in samfile.pileup(snvlabel, snvpos1, snvpos1 + 1, truncate=True, **pileup_kwargs):
                    self.filter.pileup_start(pileupcolumn)
                    total[i] += pileupcolumn.n
                    for pileupread in pileupcolumn.pileups:
                        if self.readgroups != None:
                            grp = (i,self.readgroups.group(pileupread.alignment))
                            if grp[1] == None:
                                continue
                            total[grp] += 1
                        else:
                            grp = i
                        try:
                            al, pos, base = self.filter.extract_base(pileupread)
                        except BadRead as e:
                            cnts[(grp, e.args[0])] += 1
                            continue
                        reads.append((PileupAlignment(al.seq,al.is_reverse), pos, base, grp))
                        cnts[(grp, 'Good')] += 1
                    self.filter.pileup_end(pileupcolumn)
            except ValueError as e:
                pass 
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
                assert(i in totali or len(list(totali.keys())) == 0)
                reads.extend(readsi)
                cnts.update(cntsi)
                total.update(totali)
                # self._queue[k].task_done()
                k = (k+1)%self.nt
            yield (snvchr, snvpos, ref, alt, total, reads, cnts)

