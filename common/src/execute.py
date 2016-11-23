import os, os.path, sys

class Execute(object):
    def __init__(self,*args,**kw):
	self.path = args
	self.extn = kw.get('extn')
    def addtopath(self,dir):
        self.path.append(dir)
    def setextn(self,extn):
	self.extn = extn
    def execute(self, prog, *args, **kw):
	progpath = None
	for p in self.path:
	    progpath = os.path.join(p, prog + self.extn)
	    if os.path.exists(progpath):
		break
    	assert progpath, "Exec: %s not found" % (prog + self.extn,)
    	if kw.get('verbose', False):
        	argstr = " ".join(
            	map(lambda a: a if " " not in a else '"%s"' % a, args))
        	print >>sys.stderr, "Executing:\n  %s %s" % (prog + self.extn, argstr)
    	if progpath.endswith('.py'):
        	sys.argv = [progpath] + list(args)
        	execfile(progpath,{})
    	else:
        	status = subprocess.call([progpath] + list(args))
        	assert(status == 0)
    	return True
