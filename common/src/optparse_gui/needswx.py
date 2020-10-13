
import sys, os, os.path, fnmatch, types, threading, time
import re, copy, io, csv, math, json
from optparse_gui import OptionParser, UserCancelledError, Progress
import optparse
from configparser import ConfigParser

import wx
from wx.lib.filebrowsebutton import FileBrowseButton

def quotedifnec(f):
    if ' ' in f:
        return '"%s"'%f
    return f

def quotedlistifnec(lf):
    retval = []
    for f in lf:
        retval.append(quotedifnec(f))
    return " ".join(retval)

class MyFileBrowseButton( FileBrowseButton ):
    def __init__(self,*args,**kw):
        if 'dotfile' in kw:
            self.dotfile = kw['dotfile']
            del kw['dotfile']
        if 'key' in kw:
            self.key = kw['key']
            del kw['key']
        self.isdir = False
        if 'isdir' in kw:
            self.isdir = kw['isdir']
            del kw['isdir']
        super(MyFileBrowseButton,self).__init__(*args,**kw)
    
    def createDialog( self, parent, id, pos, size, style, name="" ):
        """Setup the graphic representation of the dialog"""
        wx.Panel.__init__ (self, parent, id, pos, size, style, name )
        self.SetMinSize(size) # play nice with sizers

        box = wx.BoxSizer(wx.HORIZONTAL)

        # self.label = self.createLabel( )
        # box.Add( self.label, 0, wx.CENTER )

        self.textControl = self.createTextControl()
        box.Add( self.textControl, 1, wx.CENTER, 5)

        self.browseButton = self.createBrowseButton()
        box.Add( self.browseButton, 0, wx.LEFT|wx.CENTER, 5)

        self.SetAutoLayout(True)
        self.SetSizer( box )
        self.Layout()
        if type( size ) == tuple:
            size = wx.Size(*size)
        self.SetSize(-1, -1, size.width, size.height, wx.SIZE_USE_EXISTING)

    def OnBrowse (self, event = None):
        """ Going to browse for file... """
        current = self.GetValue()
        s = io.StringIO(current)
        rr = csv.reader(s,delimiter=' ',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        try:
            row = next(rr)
        except StopIteration:
            row = []
        if len(row) > 0 and os.path.exists(row[0]):
            directory,current = os.path.split(row[0])
            if len(row) > 1:
                current = []
                for r in row:
                    current.append(os.path.split(r)[1])
                current = ' '.join(map(quotedifnec,current))
                current = ""
        elif hasattr(self,'dotfile') and self.dotfile:
            config=ConfigParser()
            if os.path.exists(self.dotfile):
                config.read(self.dotfile)
            if config.has_section("LastFolder") and config.has_option("LastFolder",self.key):
                directory = config.get("LastFolder",self.key)
            else:
                directory = self.startDirectory
            current = ""
        else:
            directory = self.startDirectory
            current = ""


        if self.isdir:
            dlg = wx.DirDialog(self, self.dialogTitle, directory,
                               self.fileMode)
        else:
            dlg = wx.FileDialog(self, self.dialogTitle, directory, current,
                                self.fileMask, self.fileMode)

        if dlg.ShowModal() == wx.ID_OK:
            s = io.StringIO()
            wr = csv.writer(s,delimiter=' ',quotechar='"',quoting=csv.QUOTE_MINIMAL)
            if self.fileMode&wx.FD_MULTIPLE:
                wr.writerow(dlg.GetPaths())
                dir = os.path.split(dlg.GetPaths()[0])[0]
            else:
                wr.writerow([dlg.GetPath()])
                dir = os.path.split(dlg.GetPath())[0]
            self.SetValue(s.getvalue().strip())
            s.close()
            if hasattr(self,'dotfile') and self.dotfile:
                config=ConfigParser()        
                if os.path.exists(self.dotfile):
                    config.read(self.dotfile)
                if not config.has_section("LastFolder"):
                    config.add_section("LastFolder")
                config.set("LastFolder",self.key,dir)
                try:
                    wh = open(self.dotfile,'w')
                    config.write(wh)
                    wh.close()
                except IOError:
                    pass
        dlg.Destroy()


class OptparseDialog( wx.Dialog ):
    '''The dialog presented to the user with dynamically generated controls,
    to fill in the required options.
    Based on the wx.Dialog sample from wx Docs & Demos'''
    def __init__(
            self,
            option_parser, #The OptionParser object
            parent = None, 
            ID = 0, 
            title = 'Program Options', 
            pos=wx.DefaultPosition, 
            size=wx.DefaultSize, 
            style=wx.DEFAULT_DIALOG_STYLE,
            name = 'OptparseDialog',
            values = None,
            args = False
        ):

        self.option_parser = option_parser
        if values == None:
            values = option_parser.get_defaults()

        provider = wx.SimpleHelpProvider()
        wx.HelpProvider.Set(provider)
        
        wx.Dialog.__init__(self)
        self.SetExtraStyle(wx.FRAME_EX_CONTEXTHELP)
        self.Create(parent, ID, title, pos, size, style)

        sizer = wx.BoxSizer(wx.VERTICAL)

        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        top_label_text = '%s %s' % ( option_parser.get_prog_name(), 
                                     option_parser.get_version() )
        label = wx.StaticText(self, -1, top_label_text)
        sizer2.Add(label, 0, wx.GROW|wx.ALIGN_LEFT|wx.ALL, 5)
        if wx.Platform != "__WXMSW__":
            sizer2.AddStretchSpacer(-1)
            btn = wx.ContextHelpButton(self)
            sizer2.Add(btn, 0, wx.ALL)
        sizer.Add(sizer2,0, wx.GROW|wx.ALL, 5)
        
        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.RIGHT|wx.TOP|wx.LEFT, 5)

        nb = wx.Notebook(self, wx.ID_ANY)
        
        self.option_controls = {}
        self.option_controls.update(self.set_panel(nb,option_parser.option_list,values,'Options'))
        for g in option_parser.option_groups:
            self.option_controls.update(self.set_panel(nb,g.option_list,values,g.title))
        if args:
            self.args_ctrl = self.set_args_panel(nb,values,'Arguments')
        else:
            self.args_ctrl = None

        sizer.Add(nb, 0, wx.GROW|wx.RIGHT|wx.TOP|wx.LEFT, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.RIGHT|wx.TOP|wx.LEFT, 5)

        btnsizer = wx.BoxSizer(wx.HORIZONTAL)

        btn = wx.Button(self, wx.ID_CANCEL)
        # btn.SetHelpText("The OK button completes the dialog")
        btnsizer.Add(btn,0,wx.ALL,5)

        btnsizer.AddSpacer(100)
        btn = wx.Button(self, wx.ID_CLEAR, label="Reset")
        btn.Bind(wx.EVT_BUTTON, self.closeDialog)
        btnsizer.Add(btn,0,wx.ALL,5)

        btnsizer.AddSpacer(100)
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        # btn.SetHelpText("The Cancel button cancels the dialog.")
        btnsizer.Add(btn,0,wx.ALL,5)

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER,wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def closeDialog(self,event):
        self.state = 'Reset'
        self.Close()

    def set_panel(self,parent,opts,values,title):
        nopt = len(opts)
        s = wx.FlexGridSizer(2)
        p = wx.Panel(parent, -1)
        parent.AddPage(p,title)
        p.SetSizer(s)
        return self.add_opts(opts,values,p,s)

    def set_args_panel(self,parent,values,title):
        s = wx.FlexGridSizer(1,2)
        p = wx.Panel(parent, -1)
        parent.AddPage(p,title)
        p.SetSizer(s)
        label = wx.StaticText(p, -1, 'Arguments' )
        label.SetHelpText( 'Free-form arguments.' )
            
        ctrl = wx.TextCtrl( p, -1, '', size = ( 300, 100 ), 
                                   style=wx.TE_MULTILINE | wx.TE_PROCESS_ENTER )
        ctrl.SetHelpText( 
'''Args can either be separated by a space or a newline
Args that contain spaces must be entered like so: "arg with sapce" 
'''
)
        ctrl.Value = values['-args-']
        s.Add( label, 0, wx.ALIGN_RIGHT | wx.ALL, 5 )
        s.Add( ctrl, 1, wx.ALIGN_LEFT | wx.ALL, 5 )
        return ctrl

    def add_opts(self,opts,values,parent,sizer):
        controls = {}
        for option in opts:
            if option.dest is None:
                continue
            
            if option.help is None:
                option.help = ''

            if option.name is None:
                option.name = option.dest.title()

            label = wx.StaticText(parent, -1, option.name )
            label.SetHelpText( option.help )
            sizer.Add( label, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT|wx.ALL, 5 )

            if 'file' == option.type:
                if not option.filetypes:
                    fileMask = 'All Files|*.*'
                else:
                    fileMask = '|'.join([ "%s (%s)|%s"%(nm,ft,ft) for nm,ft in option.filetypes])
                ctrl = MyFileBrowseButton(parent, -1,
                                          size=(300, -1),
                                          fileMode=wx.FD_OPEN|wx.FD_FILE_MUST_EXIST,
                                          fileMask=fileMask,
                                          startDirectory=os.getcwd(),
                                          dotfile=self.option_parser.dotfile,
                                          key=option.dest,
                                          initialValue=str(values.get(option.dest,"")))
                sizer.Add( ctrl, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
            elif 'files' == option.type:
                if not option.filetypes:
                    fileMask = 'All Files|*.*'
                else:
                    fileMask = '|'.join([ "%s (%s)|%s"%(nm,ft,ft) for nm,ft in option.filetypes])
                if isinstance(values.get(option.dest,""),str):
                    initStr = values.get(option.dest,"")
                else:
                    initStr = str(' '.join(v if ' ' not in v else '"%s"'%v for v in values.get(option.dest,[])))
                ctrl = MyFileBrowseButton(parent, -1,
                                          size=(300, -1),
                                          fileMode=wx.FD_OPEN|wx.FD_MULTIPLE|wx.FD_FILE_MUST_EXIST,
                                          fileMask=fileMask,
                                          startDirectory=os.getcwd(),
                                          dotfile=self.option_parser.dotfile,
                                          key=option.dest,
                                          initialValue=initStr)
                sizer.Add( ctrl, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
            elif 'savefile' == option.type:
                if not option.filetypes:
                    fileMask = 'All Files|*.*'
                else:
                    fileMask = '|'.join([ "%s (%s)|%s"%(nm,ft,ft) for nm,ft in option.filetypes])
                ctrl = MyFileBrowseButton(parent, -1,
                                          size=(300, -1),
                                          fileMode=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT,
                                          fileMask=fileMask,
                                          startDirectory=os.getcwd(),
                                          dotfile=self.option_parser.dotfile,
                                          key=option.dest,
                                          initialValue=str(values.get(option.dest,"")))
                sizer.Add( ctrl, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
            elif 'savedir' == option.type:
                ctrl = MyFileBrowseButton(parent, -1,
                                          size=(300, -1),
                                          fileMode=wx.DD_DEFAULT_STYLE,
                                          startDirectory=os.getcwd(),
                                          dotfile=self.option_parser.dotfile,
                                          key=option.dest,isdir=True,
                                          initialValue=str(values.get(option.dest,"")))
                sizer.Add( ctrl, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
            elif 'dir' == option.type:
                ctrl = MyFileBrowseButton(parent, -1,
                                          size=(300, -1),
                                          fileMode=wx.DD_DEFAULT_STYLE|wx.DD_DIR_MUST_EXIST,
                                          startDirectory=os.getcwd(),
                                          dotfile=self.option_parser.dotfile,
                                          key=option.dest,isdir=True,
                                          initialValue=str(values.get(option.dest,"")))
                sizer.Add( ctrl, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
            elif 'store' == option.action:
                if 'choice' == option.type:
                    if optparse.NO_DEFAULT == option.default:
                        option.default = option.choices[0]
                    ctrl = wx.ComboBox( 
                        parent, -1, choices = option.choices,
                        style = wx.CB_DROPDOWN | wx.CB_READONLY, 
                        size=(300,-1)
                    )
                    try:
                        ind = option.choices.index(values.get(option.dest,None))
                    except (ValueError,KeyError):
                        ind = 0
                    ctrl.Select(ind)
                elif 'multichoice' == option.type:
                    if sys.platform == 'win32':
                        perentry = 13
                        pergap = 0
                        top = 5
                        bottom = 0
                        ysize = min(len(option.multichoices),5)*perentry + \
                                   (min(len(option.multichoices),5)-1)*pergap + top + bottom
                    else:
                        perentry = 22
                        pergap = 3
                        ysize = min(len(option.multichoices),5)*perentry + \
                                   (min(len(option.multichoices),5)+1)*pergap
                    ctrl = wx.ListBox(
                        parent, -1, choices = option.multichoices,
                        style = wx.LB_EXTENDED | \
                                wx.LB_HSCROLL | wx.LB_NEEDED_SB,
                        size = (300,ysize)
                    )
                    # print >>sys.stderr, values.get(option.dest),option.multichoices
                    selected = values.get(option.dest,[])
                    if isinstance(selected,str):
                        selected = selected.split(',')
                    for val in selected:
                        try:
                            ind = option.multichoices.index(val)
                            # print >>sys.stderr, val, ind
                            ctrl.Select(ind)
                        except ValueError:
                            continue
                else:
                    if option.text:
                        ctrl = wx.TextCtrl( parent, -1, "", size = ( 300, 100 ), 
                                            style=wx.TE_MULTILINE | wx.TE_PROCESS_ENTER )
                    elif option.type == 'password':
                        ctrl = wx.TextCtrl( parent, -1, "", size=(300,-1), 
                                            style=wx.TE_PASSWORD )
                    else:
                        ctrl = wx.TextCtrl( parent, -1, "", size=(300,-1) )
                    if option.dest in values:
                        ctrl.Value = str( values[option.dest] )
                sizer.Add( ctrl, 0, wx.ALIGN_LEFT|wx.ALL, 5 )
            elif option.action in ( 'store_true', 'store_false' ):
                ctrl = wx.CheckBox( parent, -1, "", size = ( 300, -1 ) )
                if option.dest in values:
                    ctrl.Value = values[option.dest]
                sizer.Add( ctrl, 0, wx.ALIGN_LEFT|wx.ALL, 5 )
                
            ctrl.SetHelpText( option.help )
            controls[ option ] = ctrl            
        return controls

    def _getOptions( self ):
        option_values = {}
        for option, ctrl in self.option_controls.items():
            if option.type == 'multichoice':
                option_values[option] = ','.join([option.multichoices[i] for i in ctrl.GetSelections()])
            else:
                option_values[option] = ctrl.GetValue()
        
        return option_values

    def _getArgs( self ):
        if self.args_ctrl == None:
            return []
        args_buff = self.args_ctrl.Value
        args = re.findall( r'(?:((?:(?:\w|\d)+)|".*?"))\s*', args_buff )
        return args

    def getOptionsAndArgs( self ):
        '''Returns the tuple ( options, args )
        options -  a dictionary of option names and values
        args - a sequence of args'''
        
        option_values = self._getOptions()
        args = self._getArgs()
        return option_values, args

class EmptyNotNoneOptionError (optparse.OptionError):
    """
    Raised if a notNone option has no value.
    """

class UserCheckOptionError (optparse.OptionError):
    """
    Raised if a user supplied values check fails.
    """


class OptionParserGUI( OptionParser ):
    def __init__( self, *args, **kwargs ):
        if wx.GetApp() is None:
            self.app = wx.App( False )

        self.args = False
        if 'args' in kwargs:
            self.args = kwargs['args']
            del kwargs['args']

        dotfile = None
        if 'dotfile' in kwargs:
            dotfile = kwargs['dotfile']
            del kwargs['dotfile']

        OptionParser.__init__( self, *args, **kwargs )

        self.dotfile = self.find_dotfile(dotfile)

    def find_dotfile(self,base=None):
        if not base:
            base = self.get_prog_name()
        if 'HOMEPATH' in os.environ and 'HOMEDRIVE' in os.environ:
            home = os.environ['HOMEDRIVE'] + os.environ['HOMEPATH']
        elif 'HOME' in os.environ:
            home = os.environ['HOME']
        else:
            raise RuntimeError("Can't find home directory!")
        if base.endswith('.exe'):
            dotfile = base[:-4]+'.ini'
        elif base.endswith('.py'):
            dotfile = base[:-3]+'.ini'
        else:
            dotfile = base+'.ini'
        return os.path.join(home,"."+dotfile)

    def parse_args( self, args = None, values = None, opts = None ):
        '''
        This is the heart of it all - overrides optparse.OptionParser.parse_args
        @param arg is irrelevant and thus ignored, 
               it\'s here only for interface compatibility
        '''
        if opts != None:
            initvals = {}
            for g,o in self.iteropts():
                # print repr(g),repr(o),o.dest
                if o.dest and hasattr(opts,o.dest):
                    initvals[o.dest] = getattr(opts,o.dest)
            set_values = initvals
        else:
            set_values = None

        if opts == None:
            config=ConfigParser()
            if os.path.exists(self.dotfile):
                config.read(self.dotfile)
                if not config.has_section("VERSION") or config.get("VERSION","VERSION") != self.version:
                    os.unlink(self.dotfile)
                    config=ConfigParser()
            if config.has_section("LastValue"):
                for g,o in self.iteropts():
                    if o.dest and o.remember and config.has_option("LastValue",o.dest):
                        if set_values == None:
                            set_values = {}
                        value = json.loads(config.get("LastValue",o.dest))
                        if o.type == 'multichoice':
                            set_values[o.dest] = value.split(',')
                        elif o.type in ('savefile','savedir','file'):
                            set_values[o.dest] = quotedifnec(value)
                        elif o.type in ('files',):
                            set_values[o.dest] = quotedlistifnec(value)
                        else:
                            set_values[o.dest] = value

        good = False
        while not good:
            good = True
            dlg = OptparseDialog( option_parser = self, values = set_values, args = self.args )
            dlg.state = None
            dlg_result = dlg.ShowModal()

            if dlg_result == wx.ID_CANCEL and dlg.state == None:
                raise UserCancelledError( 'User has canceled' )

            if dlg_result == wx.ID_CANCEL and dlg.state == 'Reset':
                good = False
                if os.path.exists(self.dotfile):
                    os.unlink(self.dotfile)
                set_values = None
                continue        

            assert dlg_result == wx.ID_OK

            if values is None:
                values = self.get_default_values()
            
            option_values, args = dlg.getOptionsAndArgs()

            set_values = {'-args-':''}
            for option, value in option_values.items():
                set_values[option.dest] = value
            if dlg.args_ctrl:
                set_values['-args-'] = dlg.args_ctrl.Value

            optmap = {}
            for option in self.option_list:
                optmap[str(option)] = option
            for gr in self.option_groups:
                for option in gr.option_list:
                    optmap[str(option)] = option

            for option, value in option_values.items():
                if option.action in ('store_true','store_false'):
                    setattr( values, option.dest, value )
                    continue
            
                if option.takes_value() is False:
                    value = None

                if isinstance(value,str):
                    value = str(value)

                if option.notNone and (value == None or value == ''):
                    self.error("%s: notNone option is empty"%(option,),
                               option,exit=False)
                    good = False
                    break

                try:
                    option.process( option, value, values, self )
                except optparse.OptionValueError as e:
                    self.error(e.msg,option,exit=False)
                    good = False
                    break

        config=ConfigParser()        
        if os.path.exists(self.dotfile):
            config.read(self.dotfile)
        if not config.has_section("LastValue"):
            config.add_section("LastValue")
        if not config.has_section("VERSION"):
            config.add_section("VERSION")
        config.set("VERSION","VERSION",self.version)
        for g,o in self.iteropts():
            if o.remember:
                if getattr(values,o.dest) not in (None,""):
                    config.set("LastValue",o.dest,json.dumps(getattr(values,o.dest)))
                else:
                    config.remove_option("LastValue",o.dest)
        try:
            wh = open(self.dotfile,'w')
            config.write(wh)
            wh.close()
        except IOError:
            pass
        return values, args

    def error( self, msg, option=None, exit=True):
        msg = re.sub(r"u'","'",msg)
        if ':' in msg:
            msg = msg.split(':',1)[1].strip()
            if option:
                msg = option.name+": "+msg
        dlg = wx.MessageDialog( None, msg, 'Error!', wx.ICON_ERROR )
        dlg.ShowModal()
        if exit:
            sys.exit(2)
        return 

class ProgressGUI(Progress):
    def __init__(self,title,*args,**kwargs):
        super(ProgressGUI,self).__init__(*args,**kwargs)
        self.title = title
        self.dialog = None
        self.lock = threading.Lock()

    def init(self,message):
        if wx.GetApp() is None:
            self.app = wx.App( False )
        if self.dialog:
            return
        args = (self.title,message+" "*(60-len(message)), 1001) 
        t = threading.Thread(target=self.start, args=args)
        t.start()
        while True:
            self.lock.acquire()
            if self.dialog:
                self.lock.release()
                break
            self.lock.release()
            time.sleep(1)

    def start(self,*args):
        self.lock.acquire()
        self.dialog = wx.ProgressDialog(style=0,*args)
        self.lock.release()
        
    def initprogressbar(self, message):
        self.init(message)
        self.updateprogressbar(0)
        # self.dialog.Update(0,message)

    def initbar(self, message):
        self.init(message)
        self.updatebar()
        # self.dialog.Update(0,message)

    def updateprogressbar(self,value):
        self.dialog.Update(value)

    def updatebar(self):
        self.dialog.Pulse()
