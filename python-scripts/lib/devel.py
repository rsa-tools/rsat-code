"""
:NAME:
devel

:DESC:
Profiling & debuging utilities (UNIX only)

:TODO:
memory_usage universal

"""
# IMPORTS
#-------------------------------------------------------------------------------
import sys
import os
import gc
import profile
import pstats
import tempfile

#-------------------------------------------------------------------------------
# IMPORTS

# CODE
#-------------------------------------------------------------------------------
DEBUG = 1
########################################
#                                      #
#           Debug                      #
#                                      #
########################################
def debug():
    if DEBUG:
        raise

def try_exec(code, msg=''):
    try:
        sys.stderr.write(msg)
        exec(code)
        sys.stderr.write('\n')
    except:
        Devel.debug()
        sys.stderr.write('ERROR')


########################################
#                                      #
#            Profiler                  #
#                                      #
########################################
"""
PROFILE
import profile
import pstats

#profile.run('foo()')

profile.run('foo()', 'fooprof')

p = pstats.Stats('fooprof')

p.sort_stats('cumulative').print_stats(10)
p.sort_stats('time').print_stats(10)

"""
def prof(s, limit=-1):
    """
    s     --  string of function to profile: ex: 'print("1")'
    limit -- restrict output to n first calls (default=10)
    """
    tmp = tempfile.mktemp(suffix='.prof')
    profile.run(s, tmp)
    p = pstats.Stats(tmp)
    p.sort_stats('cumulative').print_stats(limit)
    #p.sort_stats('time').print_stats(10)

########################################
#                                      #
#           Memory profiling           #
#                                      #
########################################
def memory_usage(format='M'):
    """
    total memory usage in MB or kB for current process

    format -- string M: MegaBytes k: kiloBytes (default=M)
    
    LIMITATION: only work on UNIX (tested on Mac OS X, Linux)
    """
    gc.collect()
    s_in, s_out, s_err = os.popen3('ps -u -p %d' % os.getpid())
    res = s_out.readlines()[1]

    kb = int(res.split()[5])

    if format == 'M':
        return '%.3fM' % (kb / 1024.0)
    elif format == 'k':
        return '%.0fk' % (kb)
    else:
        return kb

#-------------------------------------------------------------------------------
# CODE
