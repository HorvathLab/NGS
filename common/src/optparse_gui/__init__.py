'''
A drop-in replacement for optparse ( "import optparse_gui as optparse" )
Provides an identical interface to optparse(.OptionParser),
But displays an automatically generated wx dialog in order to enter the 
options/args, instead of parsing command line arguments
'''

import sys, os, os.path, fnmatch, types, time
import re, copy, StringIO, csv, glob
import math, optparse
from optparse import OptionGroup
from datetime import timedelta

__version__ = 0.1
__revision__ = '$Id: $'

def check_multichoice(option, opt, value):
    if not value:
	return value
    for v in value.split(','):
        if v not in option.multichoices:
	    choices = ", ".join(map(repr, option.multichoices))
            raise optparse.OptionValueError(
                "option %s: invalid choice: %r (choose one or more from %s)"
                % (opt, value, choices))
    return value

def check_file(option, opt, value):
    value = value.strip('"')
    if not value:
	return value
    value = os.path.expanduser(value)
    value = os.path.expandvars(value)
    value1 = glob.glob(value)
    # value1 += glob.glob(value+'.gz')
    # value1 += glob.glob(value+'.bz2')
    if len(value1) > 1:
       raise optparse.OptionValueError(
          "option %s: Too many files selected: %s" % (opt, value))
    if len(value1) == 0:
        raise optparse.OptionValueError(
            "option %s: File does not exist: %s" % (opt, value))
    value = value1[0]
    if option.filetypes:
        match = False
        for name,globlst in option.filetypes:
            for gl in globlst.split(';'):
	      for cmp in ('','.gz','.bz2'):
                if fnmatch.fnmatch(os.path.split(value)[1],gl+cmp):
                    match = True
                    break
              if match:
		  break
            if match:
                break
        if not match:
            raise optparse.OptionValueError(
                "option %s: File %s does not match required filetypes: %s" % (opt, value, ', '.join([ "%s (%s)"%(nm,ft) for nm,ft in option.filetypes])))
    return value

def check_files(option, opt, ssv):
    s = StringIO.StringIO(ssv)
    rd = csv.reader(s,delimiter=' ',quotechar='"')
    try:
        files = iter(rd).next()
    except StopIteration:
        files = []
    s.close()
    files1 = []
    for value in files:
	value = os.path.expanduser(value)
        value = os.path.expandvars(value)
	gv = glob.glob(value)
	# gv += glob.glob(value+'.gz')
	# gv += glob.glob(value+'.bz2')
	if len(gv) == 0 and '*' not in value and '?' not in value:
	    raise optparse.OptionValueError(
                "option %s: File does not exist: %s" % (opt, value))
	files1.extend(gv)
    if len(files1) == 0 and ssv.strip():
	raise optparse.OptionValueError(
            "option %s: No files match pattern(s): %s" % (opt, ssv))
    for value in files1:
        if not os.path.isfile(value):
            raise optparse.OptionValueError(
                "option %s: File does not exist: %s" % (opt, value))
        if option.filetypes:
            match = False
            for name,glb in option.filetypes:
		for glbi in glb.split(';'):
	          for cmp in ('','.gz','.bz2'):
                    if fnmatch.fnmatch(os.path.split(value)[1],glbi+cmp):
                        match = True
                        break
		  if match:
		      break
		if match:
		    break
            if not match:
                raise optparse.OptionValueError(
                    "option %s: File %s does not match required filetypes: %s" % (opt, value, ', '.join([ "%s (%s)"%(nm,ft) for nm,ft in option.filetypes])))
    return files1

def check_savefile(option, opt, value):
    value = value.strip('"')
    if not option.notNone and not value:
        return value
    if os.path.exists(value) and not os.path.isfile(value):
        raise optparse.OptionValueError(
            "option %s: Can't overwrite path: %s" % (opt, value))
    if option.filetypes:
        match = False
        for name,glb in option.filetypes:
	  for glbi in glb.split(';'):
            if fnmatch.fnmatch(os.path.split(value)[1],glbi):
                match = True
                break
	  if match:
	    break
        if not match:
            raise optparse.OptionValueError(
                "option %s: File %s does not match required filetypes: %s" % (opt, value, ', '.join([ "%s (%s)"%(nm,ft) for nm,ft in option.filetypes])))
    return value

def check_savedir(option, opt, value):
    value = value.strip('"')
    if not option.notNone and not value:
        return value
    if os.path.exists(value) and not os.path.isdir(value):
        raise optparse.OptionValueError(
            "option %s: Can't remove path %s" % (opt, value))
    return value

def check_dir(option, opt, value):
    value = value.strip('"')
    if not option.notNone and not value:
        return value
    if not os.path.exists(value):
        raise optparse.OptionValueError(
            "option %s: Does not exist %s" % (opt, value))
    if not os.path.isdir(value):
        raise optparse.OptionValueError(
            "option %s: Not a directory %s" % (opt, value))
    return value

class Option(optparse.Option):
    ATTRS = optparse.Option.ATTRS + ['notNone','filetypes','name','text','multichoices','remember']
    TYPES = optparse.Option.TYPES + ("password","file","savefile", "dir", "savedir", "files","multichoice")
    TYPE_CHECKER = copy.copy(optparse.Option.TYPE_CHECKER)
    TYPE_CHECKER["file"] = check_file
    TYPE_CHECKER["files"] = check_files
    TYPE_CHECKER["savefile"] = check_savefile
    TYPE_CHECKER["savedir"] = check_savedir
    TYPE_CHECKER["dir"] = check_dir
    TYPE_CHECKER["multichoice"] = check_multichoice

class OptionParser( optparse.OptionParser ):
    def __init__(self, *args, **kwargs ):
        kwargs['option_class'] = Option
	if 'dotfilename' in kwargs:
	    self.dotfilename = kwargs['dotfilename']
	    del kwargs['dotfilename']
        optparse.OptionParser.__init__( self, *args, **kwargs )
        
    def check_values (self, values, args):
        for option in self.option_list:
            if (isinstance(option, Option) and
                option.notNone and
                (getattr(values,option.dest) == "" or
                 getattr(values,option.dest) == None)):
		self.error("%s is empty" % option)
        return (values, args)

    def get_defaults(self):
        values = {}
	for (g,o) in self.iteropts():
	    if o.dest != None:
                if o.default == optparse.NO_DEFAULT or \
               	   o.default == None:
                    values[o.dest] = ''
                else:
                    values[o.dest] = o.default
        values['-args-'] = ''
	return values

    def iteropts(self):
	for o in self.option_list:
	    yield (None,o)
	for og in self.option_groups:
	    for o in og.option_list:
		yield (og,o)

    def grpopts(self):
	from collections import defaultdict
	d = defaultdict(list)
	for (g,o) in self.iteropts():
	    d[g].append(o)
	return d

class UserCancelledError( Exception ):
    pass

class Progress(object):

    def __init__(self,quiet=0):
	self._quiet = 0
	self.quiet(quiet)

    def quiet(self,q):
	oldq = self._quiet
	if isinstance(q,bool):
	    self._quiet = 2*q;
	else:
	    assert isinstance(q,int)
	    self._quiet = q
	return oldq

    def message(self,message):
	if self._quiet >= 2:
	    return
	self.initbar(message,nl=True)

    def stage(self,message,max=None,min=None,elapsed=True):
	self.elapsed = elapsed
	self.max = None
	self.min = 0
	if max != None:
	    self.max = float(max)
	if min != None:
	    self.min = float(min)
        self.value = 0
	if self._quiet >= 2:
	    return
	self.start = time.time()
        if self.max:
            self.initprogressbar(message)
        else:
            self.initbar(message)
    def update(self,increment=1,newvalue=None):
	if self._quiet >= 1:
	    return
        if self.max != None:
            if newvalue != None:
                self.value = newvalue
            else:
                self.value += increment
            self.updateprogressbar(math.floor(1000*(self.value-self.min)/(self.max-self.min)))
        else:
            self.updatebar()

    def done(self):
	if self._quiet >= 1:
	    return
	if self.max != None:
	    self.doneprogressbar()
	else:
	    self.donebar()

class ProgressText(Progress):
    def __init__(self,*args,**kwargs):
	super(ProgressText,self).__init__(*args,**kwargs)
        self.handle = sys.stdout
        self.barwidth = 10
        self.maxwidth = 60
        self.symbol = "*"
        self.bs = chr(8)
	self.neednl = False
    
    def initbar(self,message,nl=False):
	if self.neednl:
	    self.handle.write('\n')
	self.neednl = False
        print >>self.handle, message,
	self.handle.flush()
        self.barpos = 0
	self.toright = True
	if nl:
	    self.handle.write('\n')
	else:
	    self.neednl = True

    @staticmethod
    def deltaformat(delta):
	sd = map(float,str(delta).split(',',1)[-1].split(':'))
	hours,minutes,seconds = sd
	days = delta.days
	if days > 0:
	    return "%d days, %d:%02d"%(days,hours,minutes)
	if hours > 0:
	    return "%d:%02d:%02d hrs"%(hours,minutes,int(seconds))
	if minutes > 0:
	    return "%d:%02d min"%(minutes,int(seconds))
	return "%.2f sec"%(seconds,)

    def donebar(self):
	if self.elapsed:
	    d = timedelta(seconds=(time.time()-self.start))
	    print >>self.handle, "(%s)"%(self.deltaformat(d),)
	else:
	    print >>self.handle, ""
	self.neednl = False
	self.handle.flush()

    def updatebar(self):
	if self.neednl:
	    self.handle.write('\n')
	    self.neednl = False
	extrabs = False
        if self.barpos + self.barwidth >= self.maxwidth:
            self.toright = False
	    extrabs = True
	elif self.barpos == 0:
	    self.toright = True
	    extrabs = True
        if self.toright:
            self.barpos += 1
            self.handle.write("%s%s%s"%(self.bs*(self.barwidth+1*extrabs)," " if not extrabs else "",self.symbol*self.barwidth))
        else:
            self.barpos -= 1
            self.handle.write("%s%s%s"%(self.bs*(self.barwidth+2),self.symbol*self.barwidth,"  " if extrabs else " "))
	self.handle.flush()

    def initprogressbar(self,message):
	if self.neednl:
	    self.handle.write('\n')
	self.neednl = False
        # print >>self.handle, message
        print >>self.handle, "%-*s->|"%(self.maxwidth-3,
                                        message[:self.maxwidth-3])
	self.handle.flush()
        self.barpos = 0

    def doneprogressbar(self):
	print >>self.handle, (self.maxwidth-self.barpos)*self.symbol,
	if self.elapsed:
	    d = timedelta(seconds=(time.time()-self.start))
	    print >>self.handle, "(%s)"%(self.deltaformat(d),)
	else:
	    print >>self.handle, ""
	self.handle.flush()

    def updateprogressbar(self,value):
        newpos = int(round(self.maxwidth*float(value)/1000))
        if newpos > self.barpos:
            self.handle.write("%s"%self.symbol*(newpos-self.barpos))
	    self.handle.flush()
            self.barpos = newpos

try:
    from needswx import *
    __gui__ = True
except ImportError:
    __gui__ = False

def GUI():
    return __gui__

################################################################################

def sample_parse_args():
    usage = "usage: %prog [options] args"
    if 1 == len( sys.argv ):
        option_parser_class = OptionParserGUI
    else:
        option_parser_class = OptionParser
        
    parser = option_parser_class( usage = usage, version='0.1' )
    parser.add_option("-f", "--file", dest="filename", default = r'c:\1.txt',
                      help="read data from FILENAME")
    parser.add_option("-a", "--action", dest="action",
                      choices = ['delete', 'copy', 'move'],
                      help="Which action do you wish to take?!")
    parser.add_option("-n", "--number", dest="number", default = 23,
                      type = 'int',
                      help="Just a number")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose")
    
    (options, args) = parser.parse_args()
    return options, args

def main():
    options, args = sample_parse_args()
    print 'args: %s' % repr( args )
    print 'options: %s' % repr( options )
    
if '__main__' == __name__:
    main()
