"""
:NAME:
cache

:DESC:
cache

"""


class MemoryCache(object):
    """
    simple function call memory cache

    
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
        #self.access += 1
        try:
            return self.data[args]
            #self.success += 1
        except KeyError:
            self.data[args] = self.f(*args)
            return self.data[args]

    def stats(self):
        if self.access == 0:
            ratio = 0.0
        else:
            ratio = self.success / float(self.access)
        return 'access:%d success:%d efficency:%.3f' % (self.access, self.success, ratio)
        
    def __repr__(self):
        return self.stats()