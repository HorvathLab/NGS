
from pkg_resources import resource_stream
from ConfigParser import SafeConfigParser
import os.path, sys, inspect
from operator import itemgetter

class Event(object):
    """
    Base class for all events, and for events to
    be detected with all four sample types.
    """
    config = None
    GDNA = None
    SDNA = None
    NRNA = None
    TRNA = None
    keys = None
    cosmic = False
    darned = False
    def __init__(self,**kw):
        if not Event.config:
            iniPath = os.path.split(__file__)[0]
            for i in range(2):
                iniFile = os.path.join(iniPath,'event.ini')
                if os.path.exists(iniFile):
                    break
                iniPath = os.path.split(iniPath)[0]
            iniFile = open(iniFile,'r')
            Event.config = SafeConfigParser()
            Event.config.optionxform = str
            Event.config.readfp(iniFile)
        self.goodkeys = set()
        self.params = dict(self.config.items(self.__class__.__name__))
        for p,vstr in self.params.items():
            try:
                self.params[p] = vstr
                self.params[p] = float(vstr)
                self.params[p] = int(vstr)
            except ValueError, TypeError:
                pass
        for k,v in kw.items():
            assert k in self.params, "Bad keyword parameter %s"%(k,)
            self.params[k] = v

    @staticmethod
    def setCounts(GDNA,SDNA,NRNA,TRNA):
        Event.GDNA = GDNA
        Event.SDNA = SDNA
        Event.NRNA = NRNA
        Event.TRNA = TRNA
        Event.keys = set(Event.GDNA)
        Event.keys.update(Event.SDNA)
        Event.keys.update(Event.NRNA)
        Event.keys.update(Event.TRNA)

    def param(self,p):
        return self.params[p]

    def apply(self,k):
        if "GDNA" in self.requires and k not in self.GDNA:
            return False
        if "SDNA" in self.requires and k not in self.SDNA:
            return False
        if "NRNA" in self.requires and k not in self.NRNA:
            return False
        if "TRNA" in self.requires and k not in self.TRNA:
            return False
        return True

    @classmethod
    def listall(cls):
        abbrevs = []
        for cls1 in cls.allEvents():
            abbrevs.append(cls1.abbrev)
        abbrevs.sort()
        return abbrevs
            
    @classmethod
    def testall(cls):
        cls.events = []
        for cls1 in cls.allEvents():
            e = cls1()
            for k in cls.keys:
                if e.apply(k):
                    e.test(k)
            cls.events.append(e)

    @classmethod
    def allEvents(cls):
        for cls1 in map(itemgetter(1),inspect.getmembers(sys.modules[__name__])):
            if inspect.isclass(cls1) and issubclass(cls1,cls) and hasattr(cls1,'keep'):
                yield cls1

    @classmethod
    def getEvent(cls):
        for ev in cls.events:
            if isinstance(ev,cls):
                return ev
        return None

class AllSamplesEvent(Event):
    requires = set(["GDNA","SDNA","NRNA","TRNA"])

    def test(self,k):
        kw = dict(gdna = self.GDNA[k],
                  sdna = self.SDNA[k],
                  nrna = self.NRNA[k],
                  trna = self.TRNA[k])
        if self.keep(**kw):
            self.goodkeys.add(k)

class NormalOnlyEvent(Event):
    requires = set(["GDNA","NRNA"])

    def test(self,k):
        kw = dict(gdna = self.GDNA[k],
                  nrna = self.NRNA[k])
        if self.keep(**kw):
            self.goodkeys.add(k)

class DNAOnlyEvent(Event):
    requires = set(["GDNA","SDNA"])

    def test(self,k):
        kw = dict(gdna = self.GDNA[k],
                  sdna = self.SDNA[k])
        if self.keep(**kw):
            self.goodkeys.add(k)

class NoGDNAEvent(Event):
    requires = set(["SDNA","NRNA","TRNA"])

    def test(self,k):
        kw = dict(sdna = self.SDNA[k],
                  nrna = self.NRNA[k],
                  trna = self.TRNA[k])
        if self.keep(**kw):
            self.goodkeys.add(k)

class SomaticVariantEvent(Event):
    abbrev = "SOM"
    cosmic = True

class RNAEditingEvent(Event):
    abbrev = "RNAed"
    darned = True
    
class TumorRNAEditingEvent(Event):
    abbrev = "T-RNAed"
    darned = True

class VSEEvent(Event):
    abbrev = "VSE"
    
class TumorVSEEvent(Event):
    abbrev = "T-VSE"

class LoHEvent(Event):
    abbrev = "LOH"

class VSLEvent(Event):
    abbrev = "VSL"

class TumorVSLEvent(Event):
    abbrev = "T-VSL"

class Somatic_DNAOnly(DNAOnlyEvent,SomaticVariantEvent):

    def keep(self,gdna,sdna):

        cond  = (gdna['HomoRefSc'] >= self.param('GDNA_HomoRefSc'))

        cond &= (sdna['HomoVarSc'] >= self.param('SDNA_HomoVarSc') or \
                 sdna['HetSc']     >= self.param('SDNA_HetSc'))

        return cond
        
class TSS_AllSamples(AllSamplesEvent,SomaticVariantEvent):

    def keep(self,gdna,sdna,nrna,trna):

        cond  = (gdna['HomoRefSc'] >= self.param('GDNA_HomoRefSc'))

        cond &= (sdna['HomoVarSc'] >= self.param('SDNA_HomoVarSc') or \
                 sdna['HetSc']     >= self.param('SDNA_HetSc'))

        cond &= (nrna['HomoRefSc'] >= self.param('NRNA_HomoRefSc'))

        cond &= (trna['HomoVarSc'] >= self.param('TRNA_HomoVarSc') or \
                 trna['HetSc']     >= self.param('TRNA_HetSc'))

        return cond
        
class RNA_Editing_AllSamples(AllSamplesEvent,RNAEditingEvent):

    def keep(self,gdna,sdna,nrna,trna):

        cond  = (gdna['HomoRefSc'] >= self.param('GDNA_HomoRefSc') and \
                 sdna['HomoRefSc'] >= self.param('SDNA_HomoRefSc'))

        cond &= (trna['HetSc']     >= self.param('TRNA_HetSc') or \
                 trna['HomoVarSc'] >= self.param('TRNA_HomoVarSc'))
        
        cond &= (nrna['HetSc']     >= self.param('NRNA_HetSc') or \
                 nrna['HomoVarSc'] >= self.param('NRNA_HomoVarSc'))
        
        return cond

class RNA_Editing_NormalOnly(NormalOnlyEvent,RNAEditingEvent):

    def keep(self,gdna,nrna):

        cond  = (gdna['RefCount']  >= self.param('GDNA_RefCount') and \
                 gdna['SNVCount']  == self.param('GDNA_SNVCount') and \
                 nrna['SNVCount']  >= self.param('NRNA_SNVCount'))

        cond &= (gdna['HomoRefSc'] >= self.param('GDNA_HomoRefSc') and \
                 nrna['HetSc']     >= self.param('NRNA_HetSc'))
        
        return cond

class Tumor_RNA_Edit_AllSamples(AllSamplesEvent,TumorRNAEditingEvent):

    def keep(self,gdna,sdna,nrna,trna):

        cond  = (trna['HomoVarSc'] >= self.param('TRNA_HomoVarSc') or \
                 trna['HetSc']     >= self.param('TRNA_HetSc'))

        cond &= (nrna['HomoRefSc'] >= self.param('NRNA_HomoRefSc'))

        cond &= (gdna['HomoRefSc'] >= self.param('GDNA_HomoRefSc'))

        cond &= (sdna['HomoRefSc'] >= self.param('SDNA_HomoRefSc'))
        
        return cond

class Tumor_RNA_Edit_NoGDNA(NoGDNAEvent,TumorRNAEditingEvent):

    def keep(self,sdna,nrna,trna):

        cond  = (trna['HomoVarSc'] >= self.param('TRNA_HomoVarSc') or \
                 trna['HetSc']     >= self.param('TRNA_HetSc'))

        cond &= (nrna['HomoRefSc'] >= self.param('NRNA_HomoRefSc'))

        cond &= (sdna['HomoRefSc'] >= self.param('SDNA_HomoRefSc'))
        
        return cond

class VSE_AllSamples(AllSamplesEvent,VSEEvent):

    def keep(self,gdna,sdna,nrna,trna):

        cond  = (trna['HomoVarSc'] >= self.param('TRNA_HomoVarSc'))

        cond &= (nrna['HomoVarSc'] >= self.param('NRNA_HomoVarSc'))

        cond &= (gdna['HetSc']     >= self.param('GDNA_HetSc'))

        cond &= (sdna['HetSc']     >= self.param('SDNA_HetSc'))
        
        return cond    

class VSE_NormalOnly(NormalOnlyEvent,VSEEvent):

    def keep(self,gdna,nrna):

        cond  = (nrna['HomoVarSc'] >= self.param('NRNA_HomoVarSc'))

        cond &= (gdna['HetSc']     >= self.param('GDNA_HetSc'))
        
        return cond    

class Tumor_VSE_AllSamples(AllSamplesEvent,TumorVSEEvent):

    def keep(self,gdna,sdna,nrna,trna):

        cond  = (trna['HomoVarSc'] >= self.param('TRNA_HomoVarSc'))

        cond &= (nrna['HetSc']     >= self.param('NRNA_HetSc'))

        cond &= (gdna['HetSc']     >= self.param('GDNA_HetSc'))

        cond &= (sdna['HetSc']     >= self.param('SDNA_HetSc'))
        
        return cond    

class Tumor_VSE_NoGDNA(NoGDNAEvent,TumorVSEEvent):

    def keep(self,sdna,nrna,trna):

        cond  = (trna['HomoVarSc'] >= self.param('TRNA_HomoVarSc'))

        cond &= (nrna['HetSc']     >= self.param('NRNA_HetSc'))

        cond &= (sdna['HetSc']     >= self.param('SDNA_HetSc'))
        
        return cond    

class LOH_AllSamples(AllSamplesEvent,LoHEvent):

    def keep(self,gdna,sdna,nrna,trna):

        cond  = ((gdna['HetSc']     >= self.param('GDNA_HetSc')) and \
                 (sdna['HomoVarSc'] >= self.param('SDNA_HomoVarSc')))

        cond &= ((nrna['HetSc']     >= self.param('NRNA_HetSc')) or \
                 (nrna['HomoVarSc'] >= self.param('NRNA_HomoVarSc')) or \
                 (nrna['HomoRefSc'] >= self.param('NRNA_HomoRefSc')))

        cond &= (trna['HomoVarSc']  >= self.param('TRNA_HomoVarSc'))

        return cond    
                        
class LOH_DNAOnly(DNAOnlyEvent,LoHEvent):

    def keep(self,gdna,sdna):

        cond  = ((gdna['HetSc']     >= self.param('GDNA_HetSc')) and \
                 (sdna['HomoVarSc'] >= self.param('SDNA_HomoVarSc')))

        return cond    

class VSL_AllSamples(AllSamplesEvent,VSLEvent):

    def keep(self,gdna,sdna,nrna,trna):

        cond  = ((trna['HomoRefSc'] >= self.param('TRNA_HomoRefSc')) and \
                 (nrna['HomoRefSc'] >= self.param('NRNA_HomoRefSc')))

        cond &= ((gdna['HetSc']     >= self.param('GDNA_HetSc')) and \
                 (sdna['HetSc']     >= self.param('SDNA_HetSc')))

        return cond

class VSL_NormalOnly(NormalOnlyEvent,VSLEvent):

    def keep(self,gdna,nrna):

        cond  = (nrna['HomoRefSc'] >= self.param('NRNA_HomoRefSc'))

        cond &= (gdna['HetSc']     >= self.param('GDNA_HetSc'))

        return cond

class Tumor_VSL_AllSamples(AllSamplesEvent,TumorVSLEvent):

    def keep(self,gdna,sdna,nrna,trna):

        cond  = ((trna['HomoRefSc'] >= self.param('TRNA_HomoRefSc')) and \
                 (nrna['HetSc']     >= self.param('NRNA_HetSc')))

        cond &= ((gdna['HetSc']     >= self.param('GDNA_HetSc')) and \
                 (sdna['HetSc']     >= self.param('SDNA_HetSc')))

        return cond

class Tumor_VSL_NoGDNA(NoGDNAEvent,TumorVSLEvent):

    def keep(self,sdna,nrna,trna):

        cond  = ((trna['HomoRefSc'] >= self.param('TRNA_HomoRefSc')) and \
                 (nrna['HetSc']     >= self.param('NRNA_HetSc')))

        cond &=  (sdna['HetSc']     >= self.param('SDNA_HetSc'))

        return cond

