"""
:NAME:
cli

:DESC:
Command line interface (cli) utilities

:HELP:
    -Info
    Progressbar like utilitie
    -Trace
    Keep trace of a program run

"""

import sys
import os
import time
import re
import gzip


# PARAMETERS
#-------------------------------------------------------------------------------
VERBOSITY = 5
COLOR = 0
JUSTIFICATION_WIDTH = 36
UPDATE = 1
TRACE = 0
DIR = ''
#-------------------------------------------------------------------------------
# PARAMETERS

# OBJECTS
#-------------------------------------------------------------------------------
class Trace(object):
    def __init__(self, filemane):
        if not TRACE:
            return
        self.f = open(os.path.join(DIR, filemane), 'w')
        self.linecount = 0

    def writexy(self, x, y):
        self.write('%f\t%f' % (x,y))

    def write(self, data):
        if not TRACE:
            return
        self.f.write(str(data))
        self.f.flush()

    def writeln(self, data, count=False):
        if not TRACE:
            return
        if count:
            self.f.write('%d\t%s' % (self.linecount, str(data)))
            self.linecount += 1
        else:
            self.f.write(str(data))            
        self.write('\n')
        self.f.flush()
        

    def write_values(self, values, count=False):
        if not TRACE:
            return
        for val in values:
            self.writeln(val, count)


    def __del__(self):
        if not TRACE:
            return
        self.f.close()

class Info:
    """
    Display Info in text mode with progress meter
    Example:
    m = Info(10)
    for i in range(10):
        m('processing item%d' % i)

    """
    
    def __init__(self, length=None, showTime=False, showRemaining=True, verbosity=2):
        """
        length   -- int : lengh of iteration
                   if set -> autoincrement (you dont need to use ratio in call

        showTime -- show duration (default=False)

        """

        self.showTime = showTime
        self.showRemaining = showRemaining
        self.t = Timer()
        self.autoincrement = False
        self.verbosity = verbosity
        if length:
            self.length = length
            self.autoincrement = True
            self.i = 0

    def __call__(self, *args):
        """ Display info
        arguments cant be ratio (float in [0.0:1.0]), a message with optional tags(<b></b>)
        Examples:
        m(0.2, 'hello')
        m('hello', 0.2)
        m('hello')
        """
        if VERBOSITY < self.verbosity:
            return True
        msg = ''
        ratio = 0.0
        
        if self.autoincrement:
            self.i += 1
            ratio = self.i / float (self.length)
        for arg in args:
            if type(arg) == str:
                msg = arg
            elif type(arg) == float:
                ratio = arg
            else:
                raise TypeError
        if COLOR:
            msg = msg.replace('<b>', YELLOW)
            msg = msg.replace('</b>', RESET)
        s = ('> %-40s' % msg)
        if COLOR:
            s += RESET 


        if ratio < 1.0:
            s += progress(ratio)
        else:
            s += progress(1.0)
        if self.showTime:
            elapsed = self.t.get_value()
            remaining = max(0.0, (1-ratio) * (elapsed / ratio))

            s += ' %s' % (format_time(elapsed))
            if self.showRemaining:
                s += ' r%s' %  (format_time(remaining))
        runner(s)

        if ratio >= 1.0:
            pass
            info('')
            #runner('      ' * 10)
        
        return True

#-------------------------------------------------------------------------------
# OBJECTS


# CODE
#-------------------------------------------------------------------------------

################################
# C O L O R  C O N S T A N T S #
################################
BLACK = '\033[30m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[33m'
BLUE = '\033[34m'
MAGENTA = '\033[35m'
CYAN = '\033[36m'
WHITE = '\033[37m'

RESET = '\033[0;0m'
BOLD = '\033[1m'
REVERSE = '\033[2m'

BLACKBG = '\033[40m'
REDBG = '\033[41m'
GREENBG = '\033[42m'
YELLOWBG = '\033[43m'
BLUEBG = '\033[44m'
MAGENTABG = '\033[45m'
CYANBG = '\033[46m'
WHITEBG = '\033[47m'

#
# verbosity = 0 -> no info message
#           = 1 -> min message
#
#
class Timer:
    def __init__(self):
        self.start    = time.time()
        self.duration = 0.0
        self.stopped  = False
        
    def stop(self):
        if not self.stopped:
            self.duration = time.time() - self.start
            self.stopped = True

    def get_value(self):
        return time.time() - self.start    

    def __repr__(self):
        self.duration = time.time() - self.start
        gm = time.gmtime(self.duration)
        gms = time.strftime('%Hh%Mm%Ss', gm)
        return '%.3fs (%s)' % (self.duration, gms)


def begin(str):
    """Start an action
    str -- string message to display
    """
    global JUSTIFICATION_WIDTH
    #sys.stderr.write('\033[K')
    JUSTIFICATION_WIDTH = max(JUSTIFICATION_WIDTH, len(str) + 1)
    str = (str + " ").ljust(JUSTIFICATION_WIDTH)
    sys.stderr.write(str)
    sys.stderr.flush()


def end(msg='OK', info=''):
    """end of action
    msg  -- message to display (default=OK)
    info -- display this info is msg is not OK
    
    """
    if  msg != 'OK':
	if COLOR:
      	    sys.stderr.write("[" + RED + "%s" % msg + RESET + "]\n")
	else:
            sys.stderr.write("[" + "%s" % msg  + "]\n")
    else:
	if COLOR:
            sys.stderr.write("[" + GREEN + "OK" + RESET + "] " + info + "\n")
	else:
            sys.stderr.write("[" + "OK" + "] " + info + "\n")
    sys.stderr.flush()

def warning(msg):
    """Print a warning message
    msg -- string message
    """
    info("Warning: ["  + "%s" % msg  + "]", 1)

def error(msg):
    """Print an error message
    msg -- string message
    """
    info("Error: ["  + "%s" % msg  + "]", 0)


def runner(msg, level=1):
    LINE_WIDTH = 120

    if VERBOSITY >= level:	
        if len(msg) > LINE_WIDTH:
            msg = msg[:LINE_WIDTH-3] + '...'
        if UPDATE:

            sys.stderr.write('%-120s\033[' % msg + str(LINE_WIDTH) + 'D')
        else:
            sys.stderr.write('%-120s\n' % msg)

        sys.stderr.flush()  

        
def info(msg, level=1):
    if VERBOSITY >= level:	
        sys.stderr.write(str(msg) + '\n')
        sys.stderr.flush()        
        #sys.stderr.write('\033[' + str(len(msg)) + 'D')


def msg(msg, level=1):
    if VERBOSITY >= level:	
        sys.stderr.write(str(msg) + '\n')
        sys.stderr.flush()        


def progress(ratio):
    """percent value
    Keyword arguments:
    value  -- int
    max    -- int
    return string
    
    """
    if COLOR:
        return '[%s%3s%%%s]' % (BLUE, '%3d' % int( 100.0 * ratio), RESET)
    else:
        return '[%3s%%]' % ('%3d' % int( 100.0 * ratio))


def format_time(t):
    h,mod = divmod(t, 60*60)
    m,s = divmod(mod, 60)
    #strt =  '%.3fs (%dh%02dm%02s)' % (t, h,m,s)
    #strt =  '%dh%02dm%5.2fs' % (h,m,s)
    strt =  '%dh%02dm%02ds' % (h,m,s)

    if COLOR:
        return '[%s%s%s]' % (MAGENTA, strt, RESET)
    else:
        return '[%s]' % strt


def auto_script(function, usage='', help=''):
    """funtion to cli script
    """
    import optparse

    usage = usage or '%prog [OPTIONS]\n'#function.func_code.co_name
    usage = usage + help

    parser  = optparse.OptionParser(usage)

    argNames = function.func_code.co_varnames[:function.func_code.co_argcount]
    for argName in argNames:

        try:
            parser.add_option('-' + argName[0], '--' + argName, dest=argName, action='store')
        except optparse.OptionConflictError:
            parser.add_option( '--' + argName, dest=argName, action='store')

    (options, args) = parser.parse_args()

    kwargs = {}
    for argName in argNames:
        if getattr(options, argName) != None:
            kwargs[argName] = getattr(options, argName)

    outputstr = function(**kwargs)
    if outputstr != None:
        print sys.stderr, outputstr

    
        
#-------------------------------------------------------------------------------
# CODE


# TEST
#-------------------------------------------------------------------------------
def test_runner():
    for i in range(10):
	if i == 0:
		runner('%d 123' % i)
	else:
		runner('%d kkkkk' % i)
        time.sleep(1)

def test_info(x=100):
    
    m = Info(x, 1)
    for i in range(x):

        time.sleep(1)
        m('processing item%d' % i)
#-------------------------------------------------------------------------------
# TEST





        
