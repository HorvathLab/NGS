
from pysamimport import pysam
import re, sys
from collections import defaultdict

class ChromLabels(object):

    def __init__(self,labels):

        self.labels = list(labels)
        self.label2chrom = {}
        self.chrom2label = {}
        self.labelorder = {}
        self.scheme = self.extract_labels(self.labels)

    def addlabel(self,index,label,chrom):
        assert(label not in self.label2chrom), "Repeated chromosome label %s for chromosome %s"%(label,chrom)
        self.label2chrom[label] = chrom
        self.labelorder[label] = index
        assert(chrom not in self.chrom2label), "Repeated chromosome %s with label %s"%(chrom,label)
        self.chrom2label[chrom] = label

    def extract_labels(self,labels):
        scheme1 = 0; scheme2 = 0; noscheme = 0;
        for i,label in enumerate(labels):
            m = re.search(r"^chr(\d\d?|[XY]|MT?)$",label,re.IGNORECASE)
            if m:
                scheme1 += 1
                try:
                    chr = int(m.group(1))
                except ValueError:
                    chr = m.group(1).upper()
                if chr == "MT":
                    chr = "M"
                self.addlabel(i,label,chr)
                continue
            m = re.search(r"^(\d|\d\d|[XY]|MT?)$",label,re.IGNORECASE)
            if m:
                scheme2 += 1
                try:
                    chr = int(m.group(1))
                except ValueError:
                    chr = m.group(1).upper()
                if chr == "MT":
                    chr = "M"
                self.addlabel(i,label,chr)
                continue
            noscheme += 1
            self.addlabel(i,label,label)

        assert max(scheme1,scheme2) > 0, "No chromosome names in %s match expected naming schemes"%(filename,)
        assert min(scheme1,scheme2) == 0, "Multiple chromosome naming schemes in %s"%(filename,)
        if scheme1 > 0:
            return 1
        elif scheme2 > 0:
            return 2
        return 0

class ChromLabelRegistry(object):
    def __init__(self):
	self._reg = {}
        self._bam = []

    def add_labels(self,filename,labels):
        self._reg[filename] = ChromLabels(labels)

    def add_bamlabels(self,filename):
        samfile = pysam.Samfile(filename, "rb")
	assert samfile._hasIndex(), "Cannot open BAM index for file %s"%filename
	self._reg[filename] = ChromLabels(samfile.references)
        self._bam.append(filename)

    def label2chrom(self, filename, label):
        return self._reg[filename].label2chrom.get(label)

    def chrom2label(self, filename, chrom):
        return self._reg[filename].chrom2label.get(chrom)

    def chrommap(self, filename):
        return self._reg[filename].chrom2label.get

    def label2order(self, filename, label):
        return self._reg[filename].labelorder.get(label)

    def labels(self, filename):
        return self._reg[filename].labels

    def chroms(self, filename):
	for lab in self.labels(filename):
	    yield self.label2chrom(filename,lab)

    def label2label(self, filename1, filename2, label1):
        chr1 = self.label2chrom(filename1,label1)
        if not chr1:
            return None
        return self.chrom2label(filename2, chr1)

    def chrom_order(self, chrom):
        assert hasattr(self,'_chrom_order')
        return self._chrom_order.get(chrom,1+20)

    def determine_chrom_order(self):
        if self.consistent_bamfile_order():
            self.bamfile_chrom_order()
        else:
	    # print >>sys.stderr, "Warning: inconsistent chromosome name order"
            self.default_chrom_order()

    def consistent_bamfile_order(self):
        for i1 in range(len(self._bam)):
            bf1 = self._bam[i1]
            for i2 in range(len(self._bam)):
                bf2 = self._bam[i2]
                if i1 >= i2:
                    continue
                labels = filter(lambda l: self.label2label(bf2,bf1,l) != None, self.labels(bf2))
                for j2 in range(1,len(labels)):
                    ord0 = self.label2order(bf1,self.label2label(bf2,bf1,labels[j2-1]))
                    ord1 = self.label2order(bf1,self.label2label(bf2,bf1,labels[j2]))
                    if ord0 >= ord1:
                        return False
        return True

    def bamfile_chrom_order(self):
	allchrom = set()
	inedges = defaultdict(set)
	outedges = defaultdict(set)
	s = set()
	for bf in self._bam:
	    bfchrs = list(self.chroms(bf))
	    allchrom.update(bfchrs)
	    for i in range(1,len(bfchrs)):
		inedges[bfchrs[i]].add(bfchrs[i-1])
		outedges[bfchrs[i-1]].add(bfchrs[i])
	s = set()
	for chr in allchrom:
	    if len(inedges[chr]) == 0:
	        s.add(chr)
	if len(s) == 0:
	    s.add(min(allchrom))
	nextordinal = 1
	labeled = set()
	self._chrom_order = {}
	while len(s) > 0:
            u = s.pop()
            self._chrom_order[u] = nextordinal
            labeled.add(u)
            nextordinal += 1
            for v in outedges[u]:
                 inedges[v].remove(u)
                 if len(inedges[v]) == 0 and v not in labeled:
                     s.add(v)     
        return

    def default_chrom_order(self):
        allchrom = set()
        for bf in self._bam:
            allchrom.update(self._reg[bf].chrom2label.keys())
        self._chrom_order = dict((chr,i) for i,chr in enumerate(sorted(allchrom)))
        for chr in allchrom:
            if isinstance(chr,int):
                self._chrom_order[chr] = -1000+chr
            if chr == "X":
                self._chrom_order[chr] = -3
            if chr == "Y":
                self._chrom_order[chr] = -2
            if chr == "M":
                self._chrom_order[chr] = -1
        return 

    
        
