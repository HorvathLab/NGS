
import sys
from collections import defaultdict

class IntervalLookup(object):
    def __init__(self,regions=[],binsize=100):
        self._occupied = defaultdict(lambda: defaultdict(list))
        self._regions = list(regions)
        self._binsize = binsize
        self.preprocess()

    def regions(self):
        return iter(self._regions)

    def nregions(self):
        return len(self._regions)

    def preprocess(self):
        self._regions.sort()
        for i,(contig,start,stop) in enumerate(self._regions):
            if start == None or stop == None:
                continue
            startbin = start//self._binsize
            stopbin = stop//self._binsize
            for bin in range(startbin,stopbin+1):
                self._occupied[contig][bin].append(i)

    def overlap(self,contig,start,stop):
        startbin = start//self._binsize
        stopbin = stop//self._binsize
        for bin in range(startbin,stopbin+1):
            for i in self._occupied[contig][bin]:
                if start <= self._regions[i][2] and stop >= self._regions[i][1]:
                    return self._regions[i]
        return None
