"""
:NAME:
script

:DESC:
utilities for scripting

"""

from os import system
from popen2 import popen3
from commands import getoutput
from os import system
#import time


# CODE
#-------------------------------------------------------------------------------
def cmd(command):
    stdout, stdin, stderr = popen3(command)
    stdin.close()
    err = stderr.read()
    out = stdout.read()
    stdout.close()
    stderr.close()
    return out.strip()
        
#-------------------------------------------------------------------------------
# CODE


# TEST
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# TEST





        
