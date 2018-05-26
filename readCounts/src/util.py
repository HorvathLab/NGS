
import sys
from pysamimport import pysam
import re
import inspect


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

BadRead.allheaders = map(lambda cls: cls[1].header, inspect.getmembers(sys.modules[
                         __name__], lambda member: inspect.isclass(member) and issubclass(member, BadRead) and member != BadRead))

BAM_CMATCH = 0
BAM_CREF_SKIP = 3


class ReadFilter(object):
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

    def test(self, al):
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
                print >>sys.stderr, self.NONH + \
                    ".\n         Cannot filter out reads that align to mutiple loci."
                self.warnings.remove(self.NONH)
        if any(map(lambda t: t[0] not in (BAM_CMATCH, BAM_CREF_SKIP), al.cigar)):
            raise BadCigar()
        segments = [t[1] for t in al.cigar if t[0] == BAM_CMATCH]
        if len(segments) > self.maxsegments:
            raise TooManyQueryGaps()
        try:
            if int(al.opt('NM')) > self.maxedits:
                raise TooManyEdits()
        except KeyError:
            if self.NONM in self.warnings:
                print >>sys.stderr, self.NONM + \
                    ".\n         Cannot filter out reads with too many substitutions."
                self.warnings.remove(self.NONM)
        return segments


class SNVPileupReadFilter(ReadFilter):

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

    def test(self, pileupread):
        if pileupread.indel != 0:
            raise IndelAtSNV()
        if pileupread.is_del:
            raise GapAtSNV()
        al = pileupread.alignment
        segments = super(SNVPileupReadFilter,self).test(al)
        qpos = pileupread.query_position
        seg, qpos = self.findseg(qpos, segments)
        if qpos < self.minpad or (segments[seg] - qpos) < self.minpad:
            raise SNVPadding()
        try:
            edits = re.split(r'(\d+)', al.opt('MD'))[1:-1]
            substs = dict()
            reference = None
            for i in range(0, len(edits) - 1, 2):
                pos = int(edits[i])
                substs[pos] = (edits[i + 1], al.seq[pileupread.query_position])
                if pos == pileupread.query_position:
                    reference = edits[i + 1]
                elif abs(pos - pileupread.query_position) < self.minsubstdist:
                    raise SNVEditPadding()
            try:
                if int(al.opt('NM')) > (self.maxedits + (0 if (reference) else 1)):
                    raise TooManyEditsOtherThanSNV()
            except KeyError:
                if self.NONM in self.warnings:
                    print >>sys.stderr, self.NONM + \
                        ".\n         Cannot filter out reference reads with one too many\n         substitutions."
                    self.warnings.remove(self.NONM)

        except KeyError:
            if self.NOMD in self.warnings:
                print >>sys.stderr, self.NOMD + \
                    ".\n         Cannot filter out reads with edits too close to the SNV locus\n         or reference reads with one too many substitutions."
                self.warnings.remove(self.NOMD)

        readbase = al.seq[pileupread.query_position]
        return al, pileupread.query_position, readbase, segments


class NoFilter:

    def __init__(self):
        pass

    def test(self, pileupread):
        al = pileupread.alignment
        readbase = al.seq[pileupread.query_position]
        return al, pileupread.query_position, readbase, 1


class BasicFilter:

    def __init__(self):
        pass

    def test(self, pileupread):
        if pileupread.indel != 0:
            raise IndelAtSNV()
        if pileupread.is_del:
            raise GapAtSNV()
        al = pileupread.alignment
        if al.is_duplicate:
            raise IsDuplicate()
        if al.is_qcfail:
            raise IsQCFail()
        if al.is_secondary:
            raise IsSecondary()
        if al.is_unmapped:
            raise IsUnmapped()
        readbase = al.seq[pileupread.query_position]
        return al, pileupread.query_position, readbase, 1


class AllFilter:

    def __init__(self):
        pass

    def test(self, pileupread):
        raise IsBadRead()
