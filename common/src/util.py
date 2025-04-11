
import sys, os, os.path, textwrap, hashlib
from pysamimport import pysam
import re
import inspect
from configparser import ConfigParser as SafeConfigParser
from collections import defaultdict


class BadRead(RuntimeError):
    def __init__(self):
        RuntimeError.__init__(self, self.header)


class IsBadRead(BadRead):
    header = "BadRead"


class IsDuplicate(BadRead):
    header = "Alignment:IsDuplicate"


class IsQCFail(BadRead):
    header = "Alignment:IsQCFail"


class IsSecondary(BadRead):
    header = "Alignment:IsSecondary"


class IsUnmapped(BadRead):
    header = "Alignment:IsUnmapped"


class TooShort(BadRead):
    header = "TooShort"


class TooManyHits(BadRead):
    header = "MultipleAlignments"


class BadCigar(BadRead):
    header = "BadCIGAROperation"


class IndelAtSNV(BadRead):
    header = "QueryIndelAtSNVLocus"


class GapAtSNV(BadRead):
    header = "GapInQueryAtSNVLocus"


class SNVPadding(BadRead):
    header = "SNVLocusAtEndOfRead"


class SNVEditPadding(BadRead):
    header = "SubstitutionNearSNVLocus"


class TooManyEdits(BadRead):
    header = "TooManyEdits"


class TooManyEditsOtherThanSNV(BadRead):
    header = "TooManyEditsOtherThanSNV"


class TooManyQueryGaps(BadRead):
    header = "TooManyQueryGaps"


class MappingQualityTooLow(BadRead):
    header = "MappingQualityTooLow"

class BaseQualityTooLow(BadRead):
    header = "BaseQualityTooLow"

class OrphanRead(BadRead):
    header = "OrphanRead"

class OverlapRead(BadRead):
    header = "OverlapRead"

class DuplicateRead(BadRead):
    header = "DuplicateRead"

BadRead.allheaders = [cls[1].header for cls in inspect.getmembers(sys.modules[
                         __name__], lambda member: inspect.isclass(member) and issubclass(member, BadRead) and member != BadRead)]

class OtherError(RuntimeError):
    def __init__(self):
        RuntimeError.__init__(self, self.msg) 

class NoNMTag(OtherError):
    msg = "No NM tag provided for alignments, cannot filter based on edit distance."

class NoNHTag(OtherError):
    msg = "No NH tag provided for alignments, cannot filter based on number of hits."

class NoMDTag(OtherError):
    msg = "No MD tag provided for alignments, cannot filter based on position of edits."

BAM_CMATCH = 0
BAM_CREF_SKIP = 3

class ReadFilter(object):

    def pileup_kwargs(self):
        return dict(stepper='nofilter',
                    ignore_overlaps=True,
                    ignore_orphans=False,
                    flag_filter=0,
                    min_base_quality=0,
                    max_depth=100000)

    def pileup_start(self,pileupcolumn):
        pass

    def pileup_end(self,pileupcolumn):
        pass

    @staticmethod
    def extract_base_(pileupread):
        al = pileupread.alignment
        readbase = al.query_sequence[pileupread.query_position]
        return al, pileupread.query_position, readbase

class BasicReadFilter(ReadFilter):
    NONH = "Warning: Tag NH missing from alignments"
    NONM = "Warning: Tag NM missing from alignments"
    NOMD = "Warning: Tag MD missing from alignments"

    def __init__(self, maxsegments=1, minlength=45,
                 maxhits=1, maxedits=1, mapq=4,
                 warnings=set([NONM, NOMD])):
        self.minlength = minlength
        self.maxsegments = maxsegments
        self.maxhits = maxhits
        self.maxedits = maxedits
        self.warnings = warnings
        self.mapq = mapq
        if self.warnings == None:
            self.warnings = set()

    def segments(self, al):
        if al.is_duplicate:
            raise IsDuplicate()
        if al.is_qcfail:
            raise IsQCFail()
        if al.is_secondary:
            raise IsSecondary()
        if al.is_unmapped:
            raise IsUnmapped()
        if al.qlen < self.minlength:
            raise TooShort()
        if al.mapq < self.mapq:
            raise MappingQualityTooLow()
        try:
            if int(al.opt('NH')) > self.maxhits:
                raise TooManyHits()
        except KeyError:
            if self.NONH in self.warnings:
                print(self.NONH + \
                    ".\n         Cannot filter out reads that align to mutiple loci.", file=sys.stderr)
                self.warnings.remove(self.NONH)
        if any([t[0] not in (BAM_CMATCH, BAM_CREF_SKIP) for t in al.cigartuples]):
            raise BadCigar()
        segments = [t[1] for t in al.cigar if t[0] == BAM_CMATCH]
        if len(segments) > self.maxsegments:
            raise TooManyQueryGaps()
        try:
            if int(al.get_tag('NM')) > self.maxedits:
                raise TooManyEdits()
        except KeyError:
            if self.NONM in self.warnings:
                print(self.NONM + \
                    ".\n         Cannot filter out reads with too many substitutions.", file=sys.stderr)
                self.warnings.remove(self.NONM)
        return segments


class SNVPileupReadFilter(BasicReadFilter):

    def __init__(self, minpad=3, minsubstdist=3, maxedits=1, **kw):
        kw['maxedits'] = (maxedits + 1)
        self.maxedits = maxedits
        self.minpad = minpad
        self.minsubstdist = minsubstdist
        super(SNVPileupReadFilter,self).__init__(**kw)

    def findseg(self, pos, segments):
        i = 0
        while True:
            if (pos <= segments[i]):
                return i, pos
            pos -= segments[i]
            i += 1
        return None

    def extract_base(self, pileupread):
        if pileupread.indel != 0:
            raise IndelAtSNV()
        if pileupread.is_del:
            raise GapAtSNV()
        al = pileupread.alignment
        segments = self.segments(al)
        qpos = pileupread.query_position
        seg, qpos = self.findseg(qpos, segments)
        if qpos < self.minpad or (segments[seg] - qpos) < self.minpad:
            raise SNVPadding()
        try:
            edits = re.split(r'(\d+)', al.get_tag('MD'))[1:-1]
            substs = dict()
            reference = None
            for i in range(0, len(edits) - 1, 2):
                pos = int(edits[i])
                substs[pos] = (edits[i + 1], al.query_sequence[pileupread.query_position])
                if pos == pileupread.query_position:
                    reference = edits[i + 1]
                elif abs(pos - pileupread.query_position) < self.minsubstdist:
                    raise SNVEditPadding()
            try:
                if int(al.get_tag('NM')) > (self.maxedits + (0 if (reference) else 1)):
                    raise TooManyEditsOtherThanSNV()
            except KeyError:
                if self.NONM in self.warnings:
                    print(self.NONM + \
                        ".\n         Cannot filter out reference reads with one too many\n         substitutions.", file=sys.stderr)
                    self.warnings.remove(self.NONM)

        except KeyError:
            if self.NOMD in self.warnings:
                print(self.NOMD + \
                    ".\n         Cannot filter out reads with edits too close to the SNV locus\n         or reference reads with one too many substitutions.", file=sys.stderr)
                self.warnings.remove(self.NOMD)

        readbase = al.query_sequence[pileupread.query_position]
        return al, pileupread.query_position, readbase

class BasicFilter(ReadFilter):

    def __init__(self,
                 skip_duplicate=True,
                 skip_qcfail=True,
                 skip_secondary=True,
                 skip_unmapped=True):
        self._skip_duplicate = skip_duplicate
        self._skip_qcfail = skip_qcfail
        self._skip_secondary = skip_secondary
        self._skip_unmapped = skip_unmapped

    def extract_base(self, pileupread):
        if pileupread.is_del or pileupread.is_refskip:
            raise GapAtSNV()
        al = pileupread.alignment
        if self._skip_duplicate and al.is_duplicate:
            raise IsDuplicate()
        if self._skip_qcfail and al.is_qcfail:
            raise IsQCFail()
        if self._skip_secondary and al.is_secondary:
            raise IsSecondary()
        if self._skip_unmapped and al.is_unmapped:
            raise IsUnmapped()
        readbase = al.query_sequence[pileupread.query_position]
        return al, pileupread.query_position, readbase

class BaseQualityFilter(ReadFilter):

    def __init__(self, min_base_quality=None):
        self._min_base_quality = min_base_quality

    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        if self._min_base_quality != None and alignment.query_qualities[query_pos] < self._min_base_quality:
            raise BaseQualityTooLow()
        return alignment, query_pos, readbase

class MappingQualityFilter(ReadFilter):

    def __init__(self, min_mapping_quality=None):
        self._min_mapping_quality = min_mapping_quality

    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        if self._min_mapping_quality != None and alignment.mapping_quality < self._min_mapping_quality:
            raise MappingQualityTooLow()
        return alignment, query_pos, readbase

class ReadLengthFilter(ReadFilter):
    def __init__(self, min_length=None):
        self._min_length = min_length

    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        if self._min_length != None and alignment.query_length < self._min_length:
            raise TooShort()
        return alignment, query_pos, readbase

class EditsFilter(ReadFilter):
    def __init__(self, max_edits=None):
        self._max_edits = max_edits

    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        try:
            if self._max_edits != None and int(alignment.get_tag('NM')) > self._max_edits:
                raise TooManyEdits()
        except KeyError:
            raise NoNMTag()
        return alignment, query_pos, readbase

class HitsFilter(ReadFilter):
    def __init__(self, max_hits=None):
        self._max_hits = max_hits

    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        try:
            if self._max_hits != None and int(alignment.get_tag('NH')) > self._max_hits:
                raise TooManyHits()
        except KeyError:
            raise NoNHTag()
        return alignment, query_pos, readbase

class SegmentsFilter(ReadFilter):
    def __init__(self, max_segments=None):
        self._max_segments = max_segments

    @staticmethod
    def segments(alignment):
        if any([t[0] not in (BAM_CMATCH, BAM_CREF_SKIP) for t in alignment.cigar]):
            raise BadCigar()
        segments = [t[1] for t in alignment.cigar if t[0] == BAM_CMATCH]
        return segments

    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        segments = self.segments(alignment)
        if self._max_segments != None and len(segments) > self.maxsegments:
            raise TooManyQueryGaps()
        return alignment, query_pos, readbase

class EditPositionFilter(ReadFilter):

    def __init__(self, min_edge_dist=None, min_subst_dist=None, max_other_edits=None):
        self._min_edge_dist = min_edge_dist
        self._min_subst_dist = min_subst_dist
        self._max_other_edits = max_other_edits

    def findseg(self, pos, segments):
        i = 0
        while True:
            if (pos <= segments[i]):
                return i, pos
            pos -= segments[i]
            i += 1
        return None
    
    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        segments = SegmentsFilter.segments(alignment)
        seg, qpos = self.findseg(query_pos, segments)

        if self._min_edge_dist != None and \
               (qpos < self._min_edge_dist or (segments[seg] - qpos) < self._min_edge_dist):
            raise SNVEdgeDist()
        
        try:
            edits = re.split(r'(\d+)', alignment.get_tag('MD'))[1:-1]
        except KeyError:
            raise NoMDTag()
        reference = False
        for i in range(0, len(edits) - 1, 2):
            pos = int(edits[i])
            if pos == query_pos:
                reference = True
            elif self._min_subst_dist != None and abs(pos - query_pos) < self._min_subst_dist:
                raise SNVOtherEditDist()
        try:
            if self._max_other_edits != None and \
                   int(al.get_tag('NM')) > self._max_other_edits + (0 if reference else 1):
                raise TooManyEditsOtherThanSNV()
        except KeyError:
            raise NoNMTag()
        return alignment, query_pos, readbase

class OrphanFilter(ReadFilter):

    def __init__(self, remove=False):
        self._remove = remove

    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        if self._remove and alignment.is_paired and (not alignment.is_proper_pair):
            raise OrphanRead()
        return alignment, query_pos, readbase

class OverlapFilter(ReadFilter):

    def __init__(self, remove=False):
        self._remove = remove

    def pileup_start(self,pileupcolumn):
        self._seen = dict()
                
    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        if self._remove and \
               alignment.is_paired and \
               alignment.is_proper_pair and \
               alignment.query_name in self._seen:
            raise OverlapRead()
        self._seen[alignment.query_name] = True
        return alignment, query_pos, readbase
    
class UniqueReads(ReadFilter):

    def __init__(self, remove_dups=False):
        self._remove = remove_dups

    def pileup_start(self,pileupcolumn):
        self._seen = set()
                
    def extract_base(self, pileupread):
        alignment, query_pos, readbase = self.extract_base_(pileupread)
        seqhash = hashlib.md5(alignment.query_sequence.encode('utf8')).hexdigest().lower()
        if self._remove and seqhash in self._seen:
            raise DuplicateRead()
        self._seen.add(seqhash)
        return alignment, query_pos, readbase
    
class CompoundMethod(object):

    def __init__(self):
        self._elements = []
        self._description = ""
        self._name = ""
        self._type = ""
        self._specification = []

    def set_special_params(self,method,**params):
        raise NotImplemented()

    def add_element(self,element):
        self._elements.append(element)

    def add_desc(self,desc):
        self._description = desc

    def add_name(self,name):
        self._name = name

    def add_type(self,type):
        self._type = type

    def add_spec(self,spec):
        self._specification.append(spec)

    def tostr(self):
        lines = []
        if self._description:
            lines.append("Description:")
            for line in textwrap.wrap(self._description,50):
                lines.append("    "+line)
        lines.append("Specification:")
        for line in self._specification:
            lines.append("    "+line)
        return "\n".join(lines)

class CompoundFilter(CompoundMethod,ReadFilter):

    def __init__(self):
        CompoundMethod.__init__(self)
        # Defaults, unlikely to need to be changed...
        self._pileup_params = dict(stepper='nofilter',
                                   ignore_overlaps=True,
                                   ignore_orphans=False,
                                   flag_filter=0,
                                   min_base_quality=0,
                                   max_depth=100000)

    def set_special_params(self,method,**params):
        assert method == "Pileup"
        self._pileup_params = dict(params.items())

    def extract_base(self, pileupread):
        for f in self._elements:
            alignment, query_pos, readbase = f.extract_base(pileupread)
        return alignment, query_pos, readbase

    def pileup_kwargs(self):
        return self._pileup_params

    def pileup_start(self,pileupcolumn):
        for f in self._elements:
            f.pileup_start(pileupcolumn)

    def pileup_end(self,pileupcolumn):
        for f in self._elements:
            f.pileup_end(pileupcolumn)

class ReadGroup(object):

    def __init__(self, acceptlist=None, missing=None):
        self._default = missing 
        self._type = type
        self.set_acceptlist(acceptlist)

    def default(self):
        return self._default

    def type(self):
        return self._type

    validbcregex = re.compile(r'[ACGT]{6}')
    def accept(self,value):
        if not self.validbcregex.search(value):
            return False
        return (self._acceptlist == None or value in self._acceptlist)

    def set_acceptlist(self,acceptlist=None):
        self._acceptlist = None
        if acceptlist:
            try:
                self._acceptlist = set(open(acceptlist).read().split())
            except IOError:
                raise IOError("Can't open read group acceptlist file: %s"%(acceptlist))

    def group(self, alignment):
        return None

class CompoundGroup(CompoundMethod,ReadGroup):

    def group(self, alignment):
        if len(self._elements) == 1:
            return self._elements[0].group(alignment)

        grp = None
        for rg in self._elements:
            grp = rg.group(alignment)
            if grp != None:
                break
        return grp

class MethodFactory(object):
    specialMethods = []
    def __init__(self):
        iniPath = []

        progdir = os.path.split(__file__)[0]
        iniPath.append(os.path.join(progdir,self.iniFile))

        progdir = os.path.split(progdir)[0]
        iniPath.append(os.path.join(progdir,self.iniFile))

        progdir = os.path.split(progdir)[0]
        progdir = os.path.split(progdir)[0]
        iniPath.append(os.path.join(progdir,self.iniFile))

        iniPath.append(os.path.join(os.path.expanduser("~"),self.iniFile))
        iniPath.append(os.path.join(os.path.expanduser("~"),"."+self.iniFile))
        iniPath.append(os.path.join(os.getcwd(),self.iniFile))
        self.config = SafeConfigParser()
        self.config.optionxform = str
        self.config.read(iniPath)
        if len(self.config.sections()) == 0:
            print("Configuration file path:",iniPath,file=sys.stderr)
            raise RuntimeError("Can't find configuration file %s for %s"%(self.iniFile,self.__class__.__name__))

    def tovalue(self,vstr):
        vstr = vstr.strip()
        if vstr.startswith('"') and vstr.endswith('"'):
            return str(vstr[1:-1])
        if vstr.startswith("'") and vstr.endswith("'"):
            return str(vstr[1:-1])
        if vstr in ('True','False','None'):
            v = eval(vstr)
        else:
            try:
                v = str(vstr)
                v = float(vstr)
                v = int(vstr)
            except ValueError:
                pass
        return v
    
    NAME='Name'
    TYPE='Type'
    DESC='Description'
    def list(self,type=None):
        methods = []
        for sec in self.config.sections():
            if type != None:
                if not self.config.has_option(sec,self.TYPE):
                    continue
                if self.config.get(sec,self.TYPE) != type:
                    continue
            if self.config.has_option(sec,self.NAME):
                name = self.config.get(sec,self.NAME)
            else:
                name = sec
            if self.config.has_option(sec,self.DESC):
                desc = self.config.get(sec,self.DESC)
                methods.append((sec,name,desc))
            else:
                methods.append((sec,name,self.defaultDesc%(sec,)))
        methods.sort()
        return methods

    def get(self,name,params=""):

        if not self.config.has_section(name):
            raise LookupError(self.nomethodError%(name,))

        paramsbyopt = defaultdict(str)
        for line in map(str.strip,params.split(';')):
            if line == "":
                continue
            opt,rest = line.split(':',1)
            opt = opt.strip()
            value = rest.strip()
            paramsbyopt[opt] = (paramsbyopt[opt] + " " + value).strip()

        method = self.compoundMethodClass()
        for opt,value in self.config.items(name):
            opt = opt.strip()
            value = value.strip()
            if opt == self.NAME:
                method.add_name(value)
                continue
            if opt == self.TYPE:
                method.add_type(value)
                continue
            if opt == self.DESC:
                method.add_desc(value)
                continue
            kwargs = dict()
            kvpairs = []
            if paramsbyopt[opt]:
                kvpairs += re.split(r'\s+(\w+)=',' '+paramsbyopt[opt])[1:]
            elif paramsbyopt["*"]:
                kvpairs += re.split(r'\s+(\w+)=',' '+paramsbyopt["*"])[1:]                
            kvpairs += re.split(r'\s+(\w+)=',' '+value)[1:]
            seenk = set()
            for i in range(0,len(kvpairs),2):
                k = kvpairs[i]
                if k in seenk:
                    continue
                seenk.add(k)
                v = self.tovalue(kvpairs[i+1])
                vstr = str(v)
                if self.tovalue(kvpairs[i+1]) != self.tovalue(vstr):
                    vstr = '"%s"'%(v,)
                if i == 1:
                    method.add_spec("%s: %s=%s"%(opt,k,vstr))
                else:
                    method.add_spec("%s  %s=%s"%(" "*len(opt),k,vstr))
                kwargs[k] = v
            if len(kvpairs) == 1:
                method.add_spec("%s:"%(opt,))
            if opt in self.specialMethods:
                method.set_special_params(opt,**kwargs)
            else:
                try:
                    methodcls = getattr(sys.modules[__name__], opt)
                except AttributeError:
                    raise LookupError(self.noelementError%(opt,name))
                if not issubclass(methodcls,self.baseMethodClass):
                    raise LookupError(self.noelementError%(opt,name))
                try:
                    method.add_element(methodcls(**kwargs))
                except TypeError as e:
                    msg = e.args[0]
                    msg = msg.replace("__init__()",self.paramError%(opt,name))
                    e.args = tuple([msg] + list(e.args[1:]))
                    raise
                
        return method


class ReadFilterFactory(MethodFactory):

    baseMethodClass = ReadFilter
    compoundMethodClass = CompoundFilter
    iniFile = 'filter.ini'
    defaultDesc = 'Aligned read filter: %s.'
    nomethodError = "Can\'t find named read filter: %s."
    noelementError = "Can\'t find element %s of read filter %s."
    paramError = "Element %s of read filter %s"
    specialMethods = ['Pileup']

class ReadGroupFactory(MethodFactory):

    baseMethodClass = ReadGroup
    compoundMethodClass = CompoundGroup
    iniFile = 'group.ini'
    defaultDesc = 'Read group method: %s.'
    nomethodError = "Can\'t find named read group method: %s."
    noelementError = "Can\'t find element %s of read group method %s."
    paramError = "Element %s of read group method %s"

class ReadNameRegex(ReadGroup):

    def __init__(self, regex, regexgrp=1, **kw):
        super(ReadNameRegex,self).__init__(**kw)
        self._regex = re.compile(regex)
        self._regexgrp = int(regexgrp)
        

    def group(self, alignment):
        name = alignment.query_name
        m = self._regex.search(name)
        if m:
            try:
                value = m.group(self._regexgrp)
                if self.accept(value):
                    return value
            except IndexError:
                pass
        return self.default()

class ReadNameWord(ReadGroup):

    def __init__(self, field_index, field_sep='_', **kw):
        super(ReadNameWord,self).__init__(**kw)
        self._index = field_index
        self._sep = field_sep

    def group(self, alignment):
        name = alignment.query_name
        words = name.split(self._sep)
        try:
            value = words[self._index]
            if self.accept(value):
                return value
        except IndexError:
            pass
        return self.default()

class ReadTagValue(ReadGroup):

    def __init__(self, tag, **kw):
        super(ReadTagValue,self).__init__(**kw)
        self._tag = tag

    def group(self, alignment):
        try:
            value = str(alignment.get_tag(self._tag))
            if self.accept(value):
                return value
        except KeyError:
            pass
        return self.default()

class RGTag(ReadTagValue):

    def __init__(self,**kw):
        kw['tag'] = "RG"
        super(RGTag,self).__init__(**kw)

class NoValue(ReadGroup):

    def __init__(self, **kw):
        super(NoValue,self).__init__(**kw)

    def group(self, alignment):
        return self.default()



