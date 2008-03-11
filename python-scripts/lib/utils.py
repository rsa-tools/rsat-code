"""
:NAME:
utils

:DESC:
General utilities on list dict IO ...

:HELP:

"""

import time


# CODE
#-------------------------------------------------------------------------------

########################################
#                                      #
#              list
#                                      #
########################################
def unzip(l):
    """
    """
    if len(l) == 0: 
        return ()
    r = []
    for i in range(len(l[0])):
        r.append(map( lambda e: e[i], l ) )
    return tuple(r)

#python 2.3 compatibility
def sorted(l, cmpfunc=None):
    """Return sorted copy of the list
    """
    x = copy(l)
    x.sort(cmpfunc)
    return x

#python 2.3 compatibility
def reversed(l):
    """Return reversed list
    """
    x = copy(l)
    x.reverse()
    return x
    

class SortedList(list):
    """A simple sorted list
    """

    def __init__(self, max=10, growing=True):
        """
        max       -- max length
        growing   -- growing order? Boolean (True = [1,2,3] or False = [3,2,1]
        """
        list.__init__(self)
        self.max     = max
        self.growing = growing

    def inf(self, a, b):
        return a <= b

    def sup(self, a, b):
        return a >= b    

    def append(self, elt):
        # try to append element
        i = 0
        try:
            while 1:
                if self.growing:
                    if self.inf(elt,self[i]):
                        break
                else:
                    if self.sup(elt,self[i]):
                        break
                i += 1
            self.insert(i, elt)
            if len(self) > self.max:
                self.pop()

        except IndexError:
            if len(self) < self.max:
                list.append(self, elt)


    def extend(self, other):
        for i in other:
            self.append(i)

    def __iadd__(self, other):
        self.extend(other)
        return self

    def __add__(self, l):
        raise NotImplementedError


########################################
#                                      #
#             Struct
#                                      #
########################################

class Struct:
    """Struct like empty class
    """
    def __init__(self, **kargs):
        for key,value in kargs.items():
            setattr(self, key, value)    

    def values(self):
        """return values of objects in struct
        """
        return self.__dict__.values()

########################################
#                                      #
#             Struct
#                                      #
########################################
def suffix(filename):
    return filename.split('.')[-1]

def prefix(filename):
    return filename[:filename.rfind('.')]


########################################
#                                      #
#             Cache
#                                      #
########################################
class MemoryCache(object):
    """
    simple memory cache
    key value association
    """
    def __init__(self, f):
        self.data = {}
        self.f = f
        self.access = 0
        self.success = 0

    def get(self, key):
        try:
            return self.data[key]
        except KeyError:
            self.data[key] = self.f.__call__(key)
            return self.data[key]
        
    def __call__(self, *args):
        self.access += 1
        if not self.data.has_key(args):
            self.data[args] = self.f(*args)
        else:
            self.success += 1
        return self.data[args]
        
    def stats(self):
        return 'efficency: %.3f' % ( self.success / float(self.access))

########################################
#                                      #
#              IO
#                                      #
########################################
def uopen(f, mode='r'):
    """Universal open hack
    f    -- file or filename
    mode -- string mode ('r', 'w' ...)
    Return file
    """
    if type(f) is str:
        return open(f, mode)
    
    if 'U' in mode:# and 'U' not in getattr(f, 'mode', ''):
        return StringIO.StringIO('\n'.join(f.read().splitlines()))        
    
    return f
 
########################################
#                                      #
#              Timer
#                                      #
########################################
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
 
        
#-------------------------------------------------------------------------------
# CODE


# TEST
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# TEST





        
