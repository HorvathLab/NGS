import os, os.path, sys, subprocess

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
            progpath1 = os.path.join(p, prog + self.extn)
            # print(progpath1, file=sys.stderr)
            if os.path.exists(progpath1):
                progpath = progpath1
                break
        assert progpath, "Exec: %s not found" % (prog + self.extn,)
        if kw.get('verbose', False):
            argstr = " ".join([a if " " not in a else '"%s"' % a for a in args])
            print("Executing:\n  %s %s" % (prog + self.extn, argstr), file=sys.stderr)
        if progpath.endswith('.py'):
            sys.argv = [progpath] + list(args)
            exec(compile(open(progpath, "rb").read(), progpath, 'exec'),{})
        else:
            status = subprocess.call([progpath] + list(args))
            assert(status == 0)
        return True
