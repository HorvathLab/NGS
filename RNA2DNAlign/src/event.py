
from pkg_resources import resource_stream
from configparser import SafeConfigParser
import os.path, sys, inspect
from operator import itemgetter

class Event(object):
    """
    Base class for all events
    """
    config = None
    GDNA = None
    SDNA = None
    NRNA = None
    TRNA = None
    keys = None
    cosmic = False
    darned = False
    def __init__(self):
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
        self.conditions = dict(self.config.items(self.__class__.__name__))

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
            if inspect.isclass(cls1) and issubclass(cls1,cls) and cls != cls1:
                yield cls1

    @classmethod
    def getEvent(cls):
        for ev in cls.events:
            if isinstance(ev,cls):
                return ev
        return None

    def test(self,k):
        keep = True
        for r in self.requires:
            keep &= eval(self.conditions.get(r,'True'),getattr(self,r)[k])
        if keep:
            self.goodkeys.add(k)

class AllSamplesEvent(Event):
    requires = set(["GDNA","SDNA","NRNA","TRNA"])

class NormalOnlyEvent(Event):
    requires = set(["GDNA","NRNA"])

class DNAOnlyEvent(Event):
    requires = set(["GDNA","SDNA"])

class NoGDNAEvent(Event):
    requires = set(["SDNA","NRNA","TRNA"])

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
    pass
        
class TSS_AllSamples(AllSamplesEvent,SomaticVariantEvent):
    pass
        
class RNA_Editing_AllSamples(AllSamplesEvent,RNAEditingEvent):
    pass

class RNA_Editing_NormalOnly(NormalOnlyEvent,RNAEditingEvent):
    pass

class Tumor_RNA_Edit_AllSamples(AllSamplesEvent,TumorRNAEditingEvent):
    pass

class Tumor_RNA_Edit_NoGDNA(NoGDNAEvent,TumorRNAEditingEvent):
    pass

class VSE_AllSamples(AllSamplesEvent,VSEEvent):
    pass

class VSE_NormalOnly(NormalOnlyEvent,VSEEvent):
    pass

class Tumor_VSE_AllSamples(AllSamplesEvent,TumorVSEEvent):
    pass

class Tumor_VSE_NoGDNA(NoGDNAEvent,TumorVSEEvent):
    pass

class LOH_AllSamples(AllSamplesEvent,LoHEvent):
    pass
                        
class LOH_DNAOnly(DNAOnlyEvent,LoHEvent):
    pass

class VSL_AllSamples(AllSamplesEvent,VSLEvent):
    pass

class VSL_NormalOnly(NormalOnlyEvent,VSLEvent):
    pass

class Tumor_VSL_AllSamples(AllSamplesEvent,TumorVSLEvent):
    pass

class Tumor_VSL_NoGDNA(NoGDNAEvent,TumorVSLEvent):
    pass
