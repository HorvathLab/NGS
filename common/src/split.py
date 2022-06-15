
from pysamimport import pysam
import re, os, hashlib
from collections import defaultdict

class SplitBAM(object):
    def __init__(self,bamfile,readgroups,batchsize=10,directory='.',index=False,regions=None,limit=None):
        self.bamfile = os.path.realpath(bamfile)
        self.bambase,self.bamextn = self.bamfile.rsplit('.',1)
        self.bambase = os.path.split(self.bambase)[1]
        self.readgroups = readgroups
        self.batchsize = batchsize
        self.directory = directory
        self.index = index
        if regions == None:
            regions = [(None,None,None)]
        self.regions = regions
        if limit == None or limit <= 0:
            limit = 1e+20
        self.limit = limit

    def normalize_readgroup(self,rg):
        return re.sub(r'[^A-Z0-9.]','_',rg)

    def readgroup_filename(self,rg):
        uniqstr = hashlib.md5(rg.encode('ascii')).hexdigest().lower()[:5]
        return os.path.join(self.directory,self.bambase + '.' + self.normalize_readgroup(rg) + "." + uniqstr + "." + self.bamextn)

    def iterator(self):
        allrg = {}
        rgindex = 0
        limit = self.limit
        seenrg = set()
        while True:
            outal = defaultdict(set)
            more=False
            samfile = pysam.AlignmentFile(self.bamfile, "rb", require_index=False)
            for reg in sorted(self.regions):
              try:
                aliter =  samfile.fetch(contig=reg[0],start=reg[1],stop=reg[2],until_eof=True)
              except ValueError:
                # Bad contig names...
                continue
              for al in aliter:
                rg = self.readgroups.group(al)
                if rg and rg not in allrg:
                    allrg[rg] = rgindex
                    rgindex += 1
                if rg and rg not in seenrg:
                    if rg not in outal:
                        if len(outal) >= min(self.batchsize,limit):
                            more = True
                            continue
                    outal[rg].add(al)
            for rg in outal:
                rgfilename = self.readgroup_filename(rg)
                outsam = pysam.AlignmentFile(rgfilename, "wb", template=samfile)
                for al in sorted(outal[rg],key=lambda al: (al.reference_id,al.reference_start,al.reference_end)):
                    outsam.write(al)
                outsam.close()
                if self.index:
                    pysam.index(rgfilename)
                yield rg,rgfilename,allrg[rg],min(self.limit,len(allrg))
            seenrg.update(outal)
            limit = limit - len(outal)
            outal = defaultdict(set)
            if not more or limit <= 0:
                break
