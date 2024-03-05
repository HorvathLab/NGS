
import zipfile, shutil, sys, re
from operator import itemgetter
from collections import defaultdict
from csv import DictReader

# 4 ints, initial byte gives number of bytes for each integer
# two bits encode 0,1,2,3 bytes (0 bytes indicates a zero).
# 1,2,3 bytes encode the number depending on magnitude.
# least significant bit first. Max number stored is 2^24 ~= 16M

def varint403(x):
    result = []
    while x > 0:
        result.append(x&0xff)
        x = (x>>8)
    return result

def compactenc403(*args):
    bl = 0
    bs = []
    for i,a in enumerate(args):
        vi = varint403(a)
        assert(len(vi) <= 3)
        bs.extend(vi)
        bl = bl | (len(vi)<<((3-i)*2))
    return bytes([ bl ] + bs)

def compactdec403(buf,p):
    bl = buf[p];
    lens = []
    for i in range(4):
        lens.append((bl>>((3-i)*2))&0x3)
    p += 1
    vals = []
    for l in lens:
        val = 0
        for b in range(l):
            val += (int(buf[p])<<(8*b))
            p += 1
        vals.append(val)
    return p,vals
        
class scrcbinwriter(object):
    def __init__(self,filename=None):
        if filename:
            self.open(filename)

    def open(self,filename):
        self.zf = zipfile.ZipFile(filename, mode='w', compression=zipfile.ZIP_DEFLATED)
        self.zfcnts = self.zf.open('counts.403.bin', mode='w')
        self.locusindex = {}
        self.locuscount = 0
        self.cbindex = {}
        self.cbcount = 0
        self.buffer = b""

    def addloci(self,loci):
        for l in map(tuple,loci):
            if l not in self.locusindex:
                self.locusindex[l] = self.locuscount
                self.locuscount += 1

    def addbarcodes(self,barcodes):
        for cb in map(tuple,barcodes):
            if cb not in self.cbindex:
                self.cbindex[cb] = self.cbcount
                self.cbcount += 1

    def writecounts(self,chr,pos,ref,alt,fn,cb,refcnt,altcnt):
        if altcnt == 0 and refcnt == 0:
            return
        if (chr,pos,ref,alt) not in self.locusindex:
            self.locusindex[(chr,pos,ref,alt)] = self.locuscount
            self.locuscount += 1
        if (fn,cb) not in self.cbindex:
            self.cbindex[(fn,cb)] = self.cbcount
            self.cbcount += 1
        self.buffer += compactenc403(self.locusindex[(chr,pos,ref,alt)],
                                     self.cbindex[(fn,cb)],
                                     refcnt,altcnt)
        if len(self.buffer) > 4096:
            self.zfcnts.write(self.buffer)
            self.buffer = b""

    def close(self):
        if len(self.buffer) > 0:
            self.zfcnts.write(self.buffer)
            self.buffer = b""
            self.zfcnts.close()
        loci = self.zf.open('loci.tsv', mode='w')
        for (chr,pos,ref,alt),v in sorted(self.locusindex.items(),key=itemgetter(1)):
            loci.write(bytes("%s:%s_%s>%s\n"%(chr,pos,ref,alt),encoding='ascii'))
        loci.close()
        cbs = self.zf.open('barcodes.tsv', mode='w')
        for (fn,cb),v in sorted(self.cbindex.items(),key=itemgetter(1)):
            if fn:
                cbs.write(bytes("%s:%s\n"%(fn,cb),encoding='ascii'))
            else:
                cbs.write(bytes("%s\n"%(cb),encoding='ascii'))
        cbs.close()
        self.zf.close()

class scrctxtwriter(object):
    headers = """
       CHROM POS REF ALT ReadGroup SNVCount RefCount VAF
    """.split()
    def __init__(self, file=None, filename=None):
        self.needclose = False
        self.wroteheaders = False
        if filename:
            self.open(filename)
        elif file:
            self.fh = file
        else:
            self.fh = sys.stdout

    def open(self,filename,mode='wt'):
        self.fh = open(filename,mode)
        self.needclose = True

    def writecounts(self,chr,pos,ref,alt,fn,cb,refcnt,altcnt):
        if not self.wroteheaders:
            self.fh.write("\t".join(self.headers)+"\n")
            self.wroteheaders = True
        if altcnt == 0 and refcnt == 0:
            return
        vaf = "%.3f"%(altcnt/(refcnt+altcnt),)
        self.fh.write("\t".join(map(str,[chr,pos,ref,alt,cb if not fn else fn+":"+cb,altcnt,refcnt,vaf]))+"\n")

    def close(self):
        if self.needclose:
            self.fh.close()

class scrcbinreader(object):
    def __init__(self,filename=None):
        if filename:
            self.open(filename)

    def open(self,filename):
        self.zf = zipfile.ZipFile(filename, mode='r')
        self.cblookup = []
        cbtsv = self.zf.open('barcodes.tsv')
        for l in map(lambda b: bytes.decode(b).strip(),cbtsv):
            if ':' in l:
                fn,cb = l.rsplit(':',1)
                self.cblookup.append((fn,cb))
            else: 
                self.cblookup.append((None,l))
        cbtsv.close()
        locitsv = self.zf.open('loci.tsv')
        self.locilookup = []
        for l in map(lambda b: bytes.decode(b).strip(),locitsv):
            locus,alleles=l.split('_')
            chr,pos = locus.split(':')
            pos = int(pos)
            ref,alt = alleles.split('>')
            self.locilookup.append((chr,pos,ref,alt))
        locitsv.close()

    def itercounts(self):
        self.zfcnts = self.zf.open('counts.403.bin', mode='r')
        buffer = self.zfcnts.read(4096)
        bufpos = 0
        while True:
            if (len(buffer)-bufpos) < 20:
                theread = self.zfcnts.read(4096)
                if len(theread) == 0:
                    break
                buffer = buffer[bufpos:]+theread
                bufpos = 0
            bufpos,vals = compactdec403(buffer,bufpos)
            yield tuple(list(self.locilookup[vals[0]]) + list(self.cblookup[vals[1]]) + vals[2:4])
        while bufpos != len(buffer):
            bufpos,vals = compactdec403(buffer,bufpos)
            yield tuple(list(self.locilookup[vals[0]]) + list(self.cblookup[vals[1]]) + vals[2:4])

    def loci(self):
        return self.locilookup

    def barcodes(self):
        return self.cblookup

    def close(self):
        self.zf.close()

class scrctxtreader(object):
    headers = """
       CHROM POS REF ALT ReadGroup SNVCount RefCount
    """.split()
    def __init__(self, file=None, filename=None):
        self.needclose = False
        if filename:
            self.open(filename)
        elif file:
            self.fh = file
        else:
            self.fh = sys.stdin

    def open(self,filename):
        self.fh = open(filename)
        self.needclose = True

    def itercounts(self):
        for row in DictReader(self.fh,dialect='excel-tab'):
            values = list(map(row.get,self.headers))
            values.insert(4,None)
            for i in (1,6,7):
                values[i] = int(values[i])
            yield tuple(values)

    def close(self):
        if self.needclose:
            self.fh.close()

class scmtxwriter(object):
    def __init__(self,filename=None,value=None):
        if filename:
            self.writemtx(filename,value)

    def writemtx(self,filename,value):
        base = check(filename,True,'brc')
        check(base+'-%s-matrix.mtx'%(value,),False,'mtx')
        reader = scrcbinreader(filename)

        if not os.path.exists(base + '-barcodes.tsv'):
            wh = open(base + '-barcodes.tsv','wt')
        else:
            wh = None
        barcodemap = {}
        ncb = 0
        for cb in reader.barcodes():
            if wh:
                wh.write(":".join(filter(None,cb))+"\n")
            ncb += 1
            barcodemap[cb] = ncb
        if wh:
            wh.close()

        if not os.path.exists(base + '-features.tsv'):
            wh = open(base + '-features.tsv','wt')
        else:
            wh = None
        locimap = {}
        nlo = 0
        for lo in reader.loci():
            if wh:
                wh.write("%s:%s_%s>%s"%lo+"\n")
            nlo += 1
            locimap[lo] = nlo
        if wh:
            wh.close()

        values = []
        valtype = "integer"
        if value.lower() == "vaf":
            valtype = "real"
        for chr,pos,ref,alt,fn,cb,refcnt,altcnt in reader.itercounts():
            loind = locimap[(chr,pos,ref,alt)]
            cbind = barcodemap[(fn,cb)]
            if value.lower() == "snvcount" and altcnt > 0:
                values.append((loind,cbind,altcnt))
            elif value.lower() == "refcount" and refcnt > 0:
                values.append((loind,cbind,refcnt))
            elif value.lower() == "totalcount" and (altcnt+refcnt) > 0:
                values.append((loind,cbind,altcnt+refcnt))
            elif value.lower() == "vaf" and altcnt > 0:
                values.append((loind,cbind,"%.6f"%(altcnt/(altcnt+refcnt),)))
        wh = open(base + '-%s-matrix.mtx'%(value,),'wt')
        wh.write("\n".join(["%%MatrixMarket matrix coordinate %s general"%(valtype,),
                            "%s %s %s"%(nlo,ncb,len(values))])+"\n")
        for t in values:
            wh.write(" ".join(map(str,t))+"\n")
        wh.close()

class scmatwriter(object):
    def __init__(self,filename=None,value=None):
        if filename:
            self.writemat(filename,value)

    def writemat(self,filename,value):
        base = check(filename,True,'brc')
        check(base+'-%s.tsv'%(value,),False,'tsv')
        reader = scrcbinreader(filename)

        if value.lower() == "vaf":
            values = defaultdict(lambda: "NA")
        else:
            values = defaultdict(lambda: 0)

        barcodemap = {}
        ncb = 0
        for cb in reader.barcodes():
            barcodemap[cb] = ncb
            ncb += 1

        locimap = {}
        nlo = 0
        for lo in reader.loci():
            locimap[lo] = nlo
            nlo += 1

        for chr,pos,ref,alt,fn,cb,refcnt,altcnt in reader.itercounts():
            if value.lower() == "snvcount" and altcnt > 0:
                values[locimap[chr,pos,ref,alt],barcodemap[fn,cb]] = altcnt
            elif value.lower() == "refcount" and refcnt > 0:
                values[locimap[chr,pos,ref,alt],barcodemap[fn,cb]] = refcnt
            elif value.lower() == "totalcount" and (altcnt+refcnt) > 0:
                values[locimap[chr,pos,ref,alt],barcodemap[fn,cb]] = (altcnt+refcnt)
            elif value.lower() == "vaf" and (altcnt+refcnt) > 0:
                values[locimap[chr,pos,ref,alt],barcodemap[fn,cb]] = "%.6f"%(altcnt/(altcnt+refcnt),)

        wh = open(base + '-%s.tsv'%(value,),'wt')
        for cb in reader.barcodes():
            wh.write("\t"+":".join(filter(None,cb)))
        wh.write("\n")
        for i,lo in enumerate(reader.loci()):
            wh.write("%s:%s_%s>%s"%lo)
            for j in range(ncb):
                wh.write("\t"+str(values[i,j]))
            wh.write("\n")
        wh.close()

class scfilter(object):
    def __init__(self,filename=None,constraint=None,outfile=None,inplace=False):
        if filename and constraint:
            self.filter(filename,constraint,outfile,inplace)

    def filter(self,filename,constraint,outfile=None,inplace=False):
        reader = scrcbinreader(filename)
        if not outfile:
            outfile1 = filename.rsplit('.',1)[0]+'.flt.brc'
        else:
            outfile1 = outfile
        writer = scrcbinwriter(outfile1)
        reader = scrcbinreader(filename)
        for counts in reader.itercounts():
            if eval(constraint,dict(refcount=counts[6],
                                    snvcount=counts[7],
                                    totalcount=counts[6]+counts[7])):
                writer.writecounts(*counts)
        writer.close()
        reader.close()
        if inplace and not outfile:
            shutil.move(outfile1,filename)

class scrccompress(object):
    def __init__(self,filename=None,outfile=None,inplace=False):
        if filename:
            self.compress(filename,outfile,inplace)

    def compress(self,filename,outfile=None,inplace=False):
        reader = scrcbinreader(filename)
        locusfreq = defaultdict(int)
        barcodefreq = defaultdict(int)
        for counts in reader.itercounts():
            locus = counts[:4]
            barcode = counts[4:6]
            locusfreq[locus] += 1
            barcodefreq[barcode] += 1
        reader.close()
        if not outfile:
            outfile1 = filename.rsplit('.',1)[0]+'.cmp.brc'
        else:
            outfile1 = outfile
        writer = scrcbinwriter(outfile1)
        loci = sorted(locusfreq.keys(),key=locusfreq.get,reverse=True)
        bcs = sorted(barcodefreq.keys(),key=barcodefreq.get,reverse=True)
        writer.addloci(loci)
        writer.addbarcodes(bcs)
        reader = scrcbinreader(filename)
        for counts in reader.itercounts():
            writer.writecounts(*counts)
        writer.close()
        reader.close()
        if inplace and not outfile:
            shutil.move(outfile1,filename)

class scrcsorted(object):
    def __init__(self,filename=None,outfile=None,inplace=False):
        if filename:
            self.sort(filename,outfile,inplace)

    def sort(self,filename,outfile=None,inplace=False):
        reader = scrcbinreader(filename)
        if not outfile:
            outfile1 = filename.rsplit('.',1)[0]+'.srt.brc'
        else:
            outfile1 = outfile
        writer = scrcbinwriter(outfile1)
        for counts in sorted(reader.itercounts()):
            writer.writecounts(*counts)
        writer.close()
        reader.close()
        if inplace and not outfile:
            shutil.move(outfile1,filename)

from random import randint

class scrcgenerator(object):
    chrom = list(map(str,range(1,23))) + ['X','Y']
    dna = 'ACGT'
    
    def __init__(self,nloci,ncb):
        self.allcbs = []
        for i in range(ncb):
            self.allcbs.append("".join([self.dna[randint(0,3)] for l in range(20)]))
        self.nloci = nloci
        self.ncb = ncb

    def itercounts(self):
        for i in range(self.nloci):
            chr = self.chrom[randint(0,len(self.chrom)-1)]
            pos = randint(0,10000000)
            ref = self.dna[randint(0,3)]
            alt = ref
            while alt == ref:
                alt = self.dna[randint(0,3)]
            ncb = randint(self.ncb//4,2*self.ncb//3)
            for j in range(ncb):
                cb = self.allcbs[randint(0,self.ncb-1)]
                refcnt = self.makecnt()
                altcnt = self.makecnt()
                yield chr,pos,ref,alt,None,cb,refcnt,altcnt

    def makecnt(self):
        bs = randint(0,3)
        val = 0
        for i in range(bs):
            val += (randint(0,255)<<(8*i))
        return val

def fnsplit(filename):
    return filename.rsplit('.',1)

def check(filename,exists,extn):
    base,ex = fnsplit(filename)
    assert(ex == extn), "File %s does not have extention .%s"%(filename,extn)
    if exists:
        assert(os.path.exists(filename)), "File %s does not exist"%(filename,)
    else:
        assert(not os.path.exists(filename)), "File %s should not exist"%(filename,)
    return base

if __name__ == "__main__":
    import sys, os, math

    if len(sys.argv) > 1:
        cmd = sys.argv[1]
        sys.argv.pop(1)
        args = sys.argv[1:]
    else:
        cmd = 'help'

    if cmd == "compress" and len(args) == 1:
        check(args[0],exists=True,extn='brc')
        compress = scrccompress(filename=args[0],inplace=True)

    elif cmd == "sort" and len(args) == 1:
        check(args[0],exists=True,extn='brc')
        sort = scrcsorted(filename=args[0],inplace=True)

    elif cmd == "tsv2bin" and len(args) in (1,2):
        base = check(args[0],exists=True,extn='tsv')
        if len(args) == 1:
            args.append(base+'.brc')
        else:
            check(args[1],exists=False,extn='brc')
        reader = scrctxtreader(filename=args[0])
        writer = scrcbinwriter(filename=args[1])
        for counts in reader.itercounts():
            writer.writecounts(*counts)
        writer.close()
        reader.close()

    elif cmd == "bin2tsv" and len(args) in (1,2):
        base = check(args[0],exists=True,extn='brc')
        if len(args) == 1:
            args.append(base+'.tsv')
        else:
            check(args[1],exists=False,extn='tsv')
        reader = scrcbinreader(filename=args[0])
        writer = scrctxtwriter(filename=args[1])
        for counts in reader.itercounts():
            writer.writecounts(*counts)
        writer.close()
        reader.close()

    elif cmd == "cat" and len(args) == 1:
        base = check(args[0],exists=True,extn='brc')
        reader = scrcbinreader(filename=args[0])
        writer = scrctxtwriter()
        for counts in reader.itercounts():
            writer.writecounts(*counts)
        writer.close()
        reader.close()

    elif cmd == "bin2mtx" and len(args) == 2:
        base = check(args[0],exists=True,extn='brc')
        assert args[1].lower() in ("vaf","snvcount","refcount","totalcount")
        mtx = scmtxwriter(filename=args[0],value=args[1])

    elif cmd == "bin2mat" and len(args) == 2:
        base = check(args[0],exists=True,extn='brc')
        assert args[1].lower() in ("vaf","snvcount","refcount","totalcount")
        mat = scmatwriter(filename=args[0],value=args[1])

    elif cmd == "filter" and len(args) in (3,4):
        base = check(args[0],exists=True,extn='brc')
        assert args[-2].lower() in ("snvcount","refcount","totalcount")
        sa1 = list(filter(None,map(str.strip,re.split(r'(-)',args[-1]))))
        if sa1[0] == '-':
            const = '%s <= %s'%(args[-2],sa1[1])
        elif sa1[1] == '-':
            const = '%s >= %s'%(args[-2],sa1[0])
        else:
            raise RuntimeError("Bad filter criteria "+args[-1])
        if len(args) == 4:
            check(args[1],exists=False,extn='brc')
            filt = scfilter(args[0],const,outfile=args[1])
        else:
            filt = scfilter(args[0],const,inplace=True)

    elif cmd == "merge" and len(args) >= 2:
        for a in args[:-1]:
            check(a,exists=True,extn='brc')
        check(args[-1],exists=False,extn='brc')
        writer = scrcbinwriter(filename=args[-1])
        for fn in args[:-1]:
            base = check(fn,exists=True,extn='brc')
            reader = scrcbinreader(filename=fn)
            for counts in reader.itercounts():
                values = list(counts)
                values[4] = base
                writer.writecounts(*values)
            reader.close()
        writer.close()

    else:
        print("""
Usage: 
    rcio.py tsv2bin <data>.tsv [ <data>.brc ]
    rcio.py bin2tsv <data>.brc [ <data>.tsv ]
    rcio.py bin2mtx <data>.brc ( snvcount | refcount | totalcount | vaf )
    rcio.py bin2mat <data>.brc ( snvcount | refcount | totalcount | vaf )
    rcio.py filter <data>.brc [ <out>.brc ] ( snvcount | refcount | totalcount ) ( n- | -n )
    rcio.py cat <data>.brc
    rcio.py sort <data>.brc
    rcio.py compress <data>.brc
    rcio.py merge <data1>.brc [ <data2>.brc ... ] <merge>.brc
        """.strip());

    # reader = scrcgenerator(nloci=100,ncb=100)
    # writer = scrcbinwriter('test.rcbin')
    # for counts in reader.itercounts():
    #     writer.writecounts(*counts)
    # writer.close()


    # reader = scrcbinreader('test.rcbin.cmp')
    # writer = scrctxtwriter()
    # for counts in reader.itercounts():
    #     writer.writecounts(*counts)
    # writer.close()
    # reader.close()
